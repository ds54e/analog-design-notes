"""MIT: frozen one-condition initial-D1 sensing-divider measurement."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np
ROOT=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def columns(recipe, array):
    ac = recipe['kind'] == 'ac'
    assert array.ndim == 2 and array.shape[1] == 1+len(recipe['vectors'])*(2 if ac else 1)
    assert np.isfinite(array).all()
    return array[:, 0], {n: array[:, 1+2*i]+1j*array[:, 2+2*i] if ac else array[:, 1+i]
                        for i, n in enumerate(recipe['vectors'])}


def bandwidth(f, h):
    """First-edition bridge definition: -3.000 dB, linear dB/log-f crossing."""
    db = 20*np.log10(abs(h)/abs(h[0]))
    crossing = np.flatnonzero((db[:-1] > -3) & (db[1:] <= -3))
    if not len(crossing):
        return dict(status='NO_CROSSING_IN_RANGE', hz=None)
    i = int(crossing[0])
    value = np.exp(np.log(f[i])+(-3-db[i])/(db[i+1]-db[i])*np.log(f[i+1]/f[i]))
    return dict(status='CROSSING', hz=float(value), bracket_hz=[float(f[i]), float(f[i+1])],
        definition='First downward -3.000-dB crossing relative to |H(1 kHz)|; dB interpolation in log frequency, as learning.bridge.v1.')


def measure_case(name, case, arrays, plan):
    assert name == 'D1-half-divider'
    tol, fixed = plan['measurements'], plan['fixed']
    assert set(arrays) == {'analysis0', 'analysis1'}
    axes, values = [], []
    kcl, span = 0., 0.
    for i, recipe in enumerate(case['analyses']):
        axis, v = columns(recipe, arrays['analysis'+str(i)])
        assert len(axis) == (1 if recipe['kind'] == 'op' else 361)
        if recipe['kind'] == 'ac':
            assert np.all(np.diff(axis) > 0)
            assert abs(axis[0]/1e3-1) < 1e-12 and abs(axis[-1]/1e9-1) < 1e-12
        for c in case['connections'].values():
            kcl = max(kcl, float(np.max(abs(sum(v['i('+s+')'] for s in c['sensors'].values())))))
            if recipe['kind'] == 'op':
                volts = [np.zeros(len(axis)) if n == '0' else v['v('+n+')'] for n in c['nodes'].values()]
                span = max(span, float(np.ptp(np.stack(volts), axis=0).max()))
        axes.append(axis)
        values.append(v)
    assert kcl <= tol['kcl_max_a'] and span <= tol['terminal_span_max_v']
    op = {k: float(x[0]) for k, x in values[0].items()}
    expected = {'@cload[capacitance]': fixed['cload_f']}
    for name in ['load', 'bias', 'source', 'upper', 'lower']:
        expected['@r'+name+'[resistance]'] = fixed['r'+name+'_ohm']
    for key, nominal in expected.items():
        assert abs(op[key]/nominal-1) <= tol['passive_parameter_relative'], key
    source_error = {n: op['v('+n+')']-value for n, value in case['expected_source_nodes_v'].items()}
    assert max(abs(x) for x in source_error.values()) <= tol['source_voltage_v']
    available = -op['i(VtM2d)']-op['i(VtM4d)']
    load = op['v(out)']/fixed['rload_ohm']
    upper = (op['v(out)']-op['v(feedback)'])/fixed['rupper_ohm']
    lower = op['v(feedback)']/fixed['rlower_ohm']
    gate = op['i(VtM2g)']
    node_kcl = dict(output=available-load-upper,
        sensing_node=upper-lower+op['i(Vinj)'],
        negative_input=op['i(Vinj)']+gate,
        positive_input=-op['i(Vin)']-op['i(VtM1g)'])
    for node in ['tail', 'mirror', 'bias', 'vdd']:
        flow = sum(op['i('+c['sensors'][t]+')'] for c in case['connections'].values()
                   for t, n in c['nodes'].items() if n == node)
        if node == 'bias':
            flow += op['v(bias)']/fixed['rbias_ohm']
        if node == 'vdd':
            flow += op['i(Vdd)']
        node_kcl[node] = flow
    assert max(abs(x) for x in node_kcl.values()) <= tol['kcl_max_a']
    assert abs(op['v(inn)']-op['v(feedback)']) <= tol['source_voltage_v']
    assert abs(op['v(inp)']-op['v(input_access)']) <= tol['source_voltage_v']
    # Native terminal/supply current paths provide independent resistor V/I routes.
    bias_current = -op['i(VtM6d)']-op['i(VtM6g)']-op['i(VtM5g)']
    vi = {'Rbias': (op['v(bias)'], bias_current, fixed['rbias_ohm']),
          'Rsource': (op['v(source)']-op['v(inp)'], -op['i(Vin)'], fixed['rsource_ohm']),
          'Rload': (op['v(out)'], available-upper, fixed['rload_ohm']),
          'Rupper': (op['v(out)']-op['v(feedback)'], available-load, fixed['rupper_ohm']),
          'Rlower': (op['v(feedback)'], available-load-gate, fixed['rlower_ohm'])}
    readbacks = {}
    for key, (voltage, current, declared) in vi.items():
        if abs(current) < tol['vi_minimum_current_a']:
            readbacks[key] = dict(status='ILL_CONDITIONED_DC_VI', current_a=current,
                declared_ohm=declared, native_parameter_ohm=op['@'+key.lower()+'[resistance]'])
        else:
            actual = voltage/current
            difference = actual/declared-1
            assert abs(difference) <= tol['passive_vi_relative'], (key, difference)
            readbacks[key] = dict(status='MEASURED', observed_v_over_i_ohm=actual,
                declared_ohm=declared, relative_difference=difference)
    beta = fixed['rlower_ohm']/(fixed['rupper_ohm']+fixed['rlower_ohm'])
    parallel = 1/(1/fixed['rupper_ohm']+1/fixed['rlower_ohm'])
    sense_from_kcl = beta*op['v(out)']-parallel*gate
    assert abs(sense_from_kcl-op['v(inn)']) <= tol['source_voltage_v']
    # Full signed stationary power, including all four MOS terminal currents.
    device_power = {}
    for name, c in case['connections'].items():
        device_power[name] = sum((0. if n == '0' else op['v('+n+')'])*op['i('+c['sensors'][t]+')']
            for t, n in c['nodes'].items())
    resistor_power = {name: voltage**2/resistance for name, (voltage, _, resistance) in vi.items()}
    port_power = dict(Vdd=op['v(vdd)']*op['i(Vdd)'], Vin=op['v(source)']*op['i(Vin)'])
    power_residual = sum(device_power.values())+sum(resistor_power.values())+sum(port_power.values())
    assert abs(power_residual) <= tol['power_residual_w']
    f, ac = axes[1], values[1]
    h = ac['v(out)']/ac['v(source)']
    ac_upper = (ac['v(out)']-ac['v(feedback)'])/fixed['rupper_ohm']
    cap_current = -ac['i(VtM2d)']-ac['i(VtM4d)']-ac['v(out)']/fixed['rload_ohm']-ac_upper
    cmeasured = cap_current/(2j*np.pi*f*ac['v(out)'])
    mask = (f >= 1e5) & (f <= 1e7)
    cap_error = float(max(abs(cmeasured[mask]/fixed['cload_f']-1)))
    assert cap_error <= tol['capacitor_ac_relative']
    assert max(abs(ac['v(source)']-1)) <= tol['source_voltage_v']
    assert max(abs(ac['v(vdd)'])) <= tol['source_voltage_v']
    sense_ac_residual = ac_upper-ac['v(feedback)']/fixed['rlower_ohm']+ac['i(Vinj)']
    assert max(abs(sense_ac_residual)) <= tol['kcl_max_a']
    transfer = ac['v(inn)'][0]/ac['v(out)'][0]
    return dict(status='NUMERICAL_CHECKS_PASS', role=case['role'], op=op,
        input_common_mode_v=(op['v(inp)']+op['v(inn)'])/2,
        input_difference_v=op['v(inp)']-op['v(inn)'],
        output_v=op['v(out)'], ideal_output_v=op['v(source)']/beta,
        error_v=op['v(out)']-op['v(source)']/beta,
        tail_v=op['v(tail)'], tail_vsd_v=op['v(vdd)']-op['v(tail)'],
        p2_vsg_v=op['v(tail)']-op['v(inn)'], p2_vsd_v=op['v(tail)']-op['v(out)'],
        n4_vds_v=op['v(out)'], tail_delivered_a=-op['i(VtM5d)'],
        supply_delivered_a=-op['i(Vdd)'], input_delivered_a=-op['i(Vin)'],
        source_branch_a=-op['i(VtM2d)'], sink_branch_a=op['i(VtM4d)'], available_a=available,
        load_current_a=load, upper_current_a=upper, lower_current_a=lower,
        sensing_gate_current_a=gate, divider_beta=beta, divider_parallel_ohm=parallel,
        sensed_voltage_from_KCL_v=sense_from_kcl,
        divider_gate_current_shift_v=-parallel*gate,
        gain_one_khz=[float(h[0].real), float(h[0].imag)], bandwidth=bandwidth(f, h),
        sensing_transfer_one_khz=[float(transfer.real), float(transfer.imag)],
        maximum_MOS_KCL_a=kcl, maximum_terminal_span_v=span, node_KCL_a=node_kcl,
        maximum_sensing_AC_KCL_a=float(max(abs(sense_ac_residual))),
        source_voltage_errors_v=source_error, passive_parameter_readback={k: op[k] for k in expected},
        resistor_vi_readback=readbacks, capacitor_ac_relative_error=cap_error,
        stationary_power_w=dict(MOS=device_power, resistors=resistor_power, ports_absorbed=port_power,
            residual=power_residual, total_source_delivery=-sum(port_power.values())),
        added_divider_resistance_ohm=fixed['rupper_ohm']+fixed['rlower_ohm'])


def inputs(source_zip):
    plan = read(ROOT/'plan.json')
    assert sha(__file__) == plan['measurement_script_sha256']
    for name, expected in plan['selected_input_sha256'].items():
        assert sha(ROOT/name) == expected, name
    source = read(ROOT/'data/source.json')
    assert sha(source_zip) == source['package_sha256']
    with zipfile.ZipFile(source_zip) as z:
        for name, expected in source['selected_zip_members'].items():
            assert hashlib.sha256(z.read('example/'+name)).hexdigest() == expected, name
        case = json.loads(z.read('example/data/interface-inputs.json'))['cases']['D1-5p']
        assert case['physical_devices'] == plan['physical_units']
        with np.load(io.BytesIO(z.read('example/data/'+case['array_file'])), allow_pickle=False) as a:
            _, values = columns(case['analyses'][1]['recipe'], a['analysis1'])
        op = {k: float(v[0]) for k, v in values.items()}
        assert op == read(ROOT/'data/prior.json')['known']['observed_op']
    return plan


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source-zip', type=Path, required=True)
    p.add_argument('--observations', type=Path, default=ROOT/'data')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    plan = inputs(a.source_zip)
    index = read(a.observations/'index.json')
    assert index['plan_sha256'] == sha(ROOT/'plan.json') and set(index['cases']) == {'D1-half-divider'}
    name = 'D1-half-divider'
    case, entry = plan['cases'][name], index['cases'][name]
    assert sha(a.observations/entry['file']) == entry['sha256']
    with np.load(a.observations/entry['file'], allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    measured = measure_case(name, case, arrays, plan)
    prior = read(ROOT/'data/prior.json')['forecast']
    comparison = dict(prior_output_v=prior['output_v'], observed_output_v=measured['output_v'],
        output_residual_v=measured['output_v']-prior['output_v'],
        prior_gain=prior['gain'], observed_gain=measured['gain_one_khz'][0],
        gain_relative_residual=measured['gain_one_khz'][0]/prior['gain']-1,
        prior_bandwidth_hz=prior['bandwidth_minus3db_hz'], observed_bandwidth_hz=measured['bandwidth']['hz'],
        bandwidth_relative_residual=None if measured['bandwidth']['hz'] is None else measured['bandwidth']['hz']/prior['bandwidth_minus3db_hz']-1)
    result = dict(schema='learning.initial-d1-feedback-divider.measure.v1', status='NUMERICAL_CHECKS_PASS',
        plan_sha256=sha(ROOT/'plan.json'), cases={name: measured}, comparison=comparison, limits=plan['limits'])
    a.output.mkdir(parents=True)
    (a.output/'summary.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    keys = ['case', 'output_v', 'error_v', 'input_common_mode_v', 'input_difference_v', 'tail_vsd_v',
            'tail_delivered_a', 'source_branch_a', 'sink_branch_a', 'load_current_a', 'upper_current_a',
            'lower_current_a', 'sensing_gate_current_a', 'supply_delivered_a', 'input_delivered_a']
    with (a.output/'operating-point.csv').open('x', newline='') as f:
        writer = csv.DictWriter(f, keys, lineterminator='\n')
        writer.writeheader()
        writer.writerow({k: name if k == 'case' else measured[k] for k in keys})
    print(json.dumps(dict(status=result['status'], comparison=comparison), indent=2))


if __name__ == '__main__':
    main()

