"""MIT: frozen one-condition initial-D1 DC receiving-load measurement."""
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
    tol = plan['measurements']
    recipes = case['analyses']
    assert set(arrays) == {'analysis'+str(i) for i in range(len(recipes))}
    axes, values = [], []
    kcl, span = 0., 0.
    for i, recipe in enumerate(recipes):
        t, v = columns(recipe, arrays['analysis'+str(i)])
        assert len(t) == (1 if recipe['kind'] == 'op' else 361)
        if recipe['kind'] == 'ac':
            assert np.all(np.diff(t) > 0)
            assert abs(t[0]/1e3-1) < 1e-12 and abs(t[-1]/1e9-1) < 1e-12
        for c in case['connections'].values():
            kcl = max(kcl, float(np.max(abs(sum(v['i('+s+')'] for s in c['sensors'].values())))))
            if recipe['kind'] == 'op':
                volts = [np.zeros(len(t)) if n == '0' else v['v('+n+')'] for n in c['nodes'].values()]
                span = max(span, float(np.ptp(np.stack(volts), axis=0).max()))
        axes.append(t)
        values.append(v)
    assert kcl <= tol['kcl_max_a'] and span <= tol['terminal_span_max_v']
    op = {k: float(x[0]) for k, x in values[0].items()}
    fixed = plan['fixed']
    expected = {'@cload[capacitance]': fixed['cload_f'], '@rload[resistance]': fixed['rload_ohm'],
                '@rbias[resistance]': fixed['rbias_ohm']}
    if name == 'D1-100k-load':
        expected['@rsource[resistance]'] = fixed['rsource_ohm']
    for key, nominal in expected.items():
        assert abs(op[key]/nominal-1) <= tol['passive_parameter_relative']
    source_error = {k: op['v('+k+')']-v for k, v in case['expected_source_nodes_v'].items()}
    assert max(abs(x) for x in source_error.values()) <= tol['source_voltage_v']
    available = -op['i(VtM2d)']-op['i(VtM4d)']
    load = op['v(out)']/fixed['rload_ohm']
    gate_feedback = op.get('i(Vfeedback)', 0.)
    clamp = op.get('i(Vout)', 0.)
    node_kcl = {'output': available-load+gate_feedback-clamp}
    for node in ('tail', 'mirror', 'bias', 'vdd'):
        flow = sum(op['i('+c['sensors'][t]+')'] for c in case['connections'].values()
                   for t, n in c['nodes'].items() if n == node)
        if node == 'bias':
            flow += op['v(bias)']/fixed['rbias_ohm']
        if node == 'vdd':
            flow += op['i(Vdd)']
        node_kcl[node] = flow
    if name == 'D1-100k-load':
        node_kcl['input'] = -op['i(Vin)']-op['i(VtM1g)']
        node_kcl['feedback'] = op['i(Vfeedback)']+op['i(VtM2g)']
    assert max(abs(v) for v in node_kcl.values()) <= tol['kcl_max_a']
    bias_current = -op['i(VtM6d)']-op['i(VtM6g)']-op['i(VtM5g)']
    vi = {'Rbias': (op['v(bias)'], bias_current, fixed['rbias_ohm']),
          'Rload': (op['v(out)'], available+gate_feedback-clamp, fixed['rload_ohm'])}
    if name == 'D1-100k-load':
        vi['Rsource'] = (op['v(source)']-op['v(inp)'], -op['i(Vin)'], fixed['rsource_ohm'])
    readbacks = {}
    for key, (voltage, current, declared) in vi.items():
        if abs(current) < tol['vi_minimum_current_a']:
            readbacks[key] = dict(status='ILL_CONDITIONED_DC_VI', current_a=current,
                declared_ohm=declared, native_parameter_ohm=op['@'+key.lower()+'[resistance]'])
        else:
            measured = voltage/current
            relative = measured/declared-1
            assert abs(relative) <= tol['passive_vi_relative'], (key, relative)
            readbacks[key] = dict(status='MEASURED', observed_v_over_i_ohm=measured,
                declared_ohm=declared, relative_difference=relative)
    result = dict(status='NUMERICAL_CHECKS_PASS', role=case['role'], op=op,
        input_common_mode_v=(op['v(inp)']+op['v(inn)'])/2,
        input_difference_v=op['v(inp)']-op['v(inn)'],
        output_v=op['v(out)'], tail_v=op['v(tail)'], tail_vsd_v=op['v(vdd)']-op['v(tail)'],
        p2_vsg_v=op['v(tail)']-op['v(inn)'], p2_vsd_v=op['v(tail)']-op['v(out)'],
        n4_vds_v=op['v(out)'], tail_delivered_a=-op['i(VtM5d)'],
        supply_delivered_a=-op['i(Vdd)'], source_branch_a=-op['i(VtM2d)'],
        sink_branch_a=op['i(VtM4d)'], available_a=available,
        net_after_resistor_a=available-load, clamp_sink_a=clamp,
        maximum_MOS_KCL_a=kcl, maximum_terminal_span_v=span, node_KCL_a=node_kcl,
        source_voltage_errors_v=source_error, passive_parameter_readback={k: op[k] for k in expected},
        resistor_vi_readback=readbacks)
    if name == 'D1-100k-load':
        f, v = axes[1], values[1]
        h = v['v(out)']/v['v(source)']
        cap = -v['i(VtM2d)']-v['i(VtM4d)']-v['v(out)']/fixed['rload_ohm']+v['i(Vfeedback)']
        cmeasured = cap/(2j*np.pi*f*v['v(out)'])
        mask = (f >= 1e5) & (f <= 1e7)
        cap_error = float(max(abs(cmeasured[mask]/fixed['cload_f']-1)))
        assert cap_error <= tol['capacitor_ac_relative']
        assert max(abs(v['v(source)']-1)) <= tol['source_voltage_v']
        assert max(abs(v['v(vdd)'])) <= tol['source_voltage_v']
        result.update(input_source_v=op['v(source)'], error_v=op['v(out)']-op['v(source)'],
            input_delivered_a=-op['i(Vin)'], gain_one_khz=[float(h[0].real), float(h[0].imag)],
            bandwidth=bandwidth(f, h), capacitor_ac_relative_error=cap_error)
    return result


def inputs(source_zip):
    plan=read(ROOT/'plan.json')
    assert sha(__file__)==plan['measurement_script_sha256']
    for name,expected in plan['selected_input_sha256'].items():assert sha(ROOT/name)==expected,name
    source=read(ROOT/'data/source.json');assert sha(source_zip)==source['package_sha256']
    with zipfile.ZipFile(source_zip) as z:
        for name,expected in source['selected_zip_members'].items():assert hashlib.sha256(z.read('example/'+name)).hexdigest()==expected,name
        case=json.loads(z.read('example/data/interface-inputs.json'))['cases']['D1-5p']
        assert case['physical_devices']==plan['physical_units']
        with np.load(io.BytesIO(z.read('example/data/'+case['array_file'])),allow_pickle=False) as a:
            _,values=columns(case['analyses'][1]['recipe'],a['analysis1'])
        op={k:float(v[0]) for k,v in values.items()}
        assert op==read(ROOT/'data/prior.json')['known']['observed_op']
    return plan


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-zip',type=Path,required=True);p.add_argument('--observations',type=Path,default=ROOT/'data');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists();plan=inputs(a.source_zip);index=read(a.observations/'index.json')
    assert index['plan_sha256']==sha(ROOT/'plan.json') and set(index['cases'])=={'D1-100k-load'}
    case=plan['cases']['D1-100k-load'];entry=index['cases']['D1-100k-load'];assert sha(a.observations/entry['file'])==entry['sha256']
    with np.load(a.observations/entry['file'],allow_pickle=False) as z:arrays={k:z[k] for k in z.files}
    measured=measure_case('D1-100k-load',case,arrays,plan);prior=read(ROOT/'data/prior.json')
    f=prior['forecast'];comparison=dict(prior_output_v=f['output_v'],observed_output_v=measured['output_v'],output_residual_v=measured['output_v']-f['output_v'],
        prior_gain_real=f['gain_from_old_native'][0],observed_gain_real=measured['gain_one_khz'][0],gain_relative_residual=measured['gain_one_khz'][0]/f['gain_from_old_native'][0]-1,
        prior_bandwidth_hz=f['bandwidth_scaled_from_old_native_hz'],observed_bandwidth_hz=measured['bandwidth']['hz'],
        bandwidth_relative_residual=None if measured['bandwidth']['hz'] is None else measured['bandwidth']['hz']/f['bandwidth_scaled_from_old_native_hz']-1,
        known_output_v=prior['known']['observed_op']['v(out)'],known_gain=prior['known']['gain_one_khz'],known_bandwidth_hz=prior['known']['bandwidth_hz'])
    result=dict(schema='learning.initial-d1-dc-load.measure.v1',status='NUMERICAL_CHECKS_PASS',plan_sha256=sha(ROOT/'plan.json'),cases={'D1-100k-load':measured},comparison=comparison,limits=plan['limits'])
    a.output.mkdir(parents=True);(a.output/'summary.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    keys=['case','input_common_mode_v','input_difference_v','output_v','error_v','tail_v','tail_vsd_v','tail_delivered_a','source_branch_a','sink_branch_a','available_a','supply_delivered_a','input_delivered_a']
    with (a.output/'operating-point.csv').open('x',newline='') as file:
        writer=csv.DictWriter(file,keys,lineterminator='\n');writer.writeheader();writer.writerow({k:'D1-100k-load' if k=='case' else measured[k] for k in keys})
    print(json.dumps(dict(status=result['status'],comparison=comparison),indent=2))


if __name__=='__main__':main()
