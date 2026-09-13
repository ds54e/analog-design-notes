"""MIT: learning.design.v1, retrospective sizing and mirror/headroom analysis.

Saved observations only. No SPICE, optimizer, random samples or private paths.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs():
    data = ROOT/'data'
    for name, expected in read(data/'design-manifest.json')['files'].items():
        assert sha(data/name) == expected, name
    config = read(data/'design-inputs.json')
    for key, case in config['cases'].items():
        assert sha(data/case['array_file']) == case['array_sha256'], key
        assert sha(ROOT/'circuits'/case['circuit_file']) == case['circuit_sha256'], key
    headroom = read(data/'headroom-inputs.json')
    for key, case in headroom['cases'].items():
        assert sha(data/('headroom-'+case['array_file'])) == case['array_sha256'], key
    return config, headroom


def columns(case, index, filename):
    desc = next(x for x in case['analyses'] if x.get('original_analysis_index', x.get('analysis_index')) == index)
    r = desc['recipe']
    with np.load(filename, allow_pickle=False) as z:
        a = z['analysis'+str(index)]
    assert a.shape == (desc['rows'], desc['columns']) and np.isfinite(a).all()
    values = [a[:, 1+2*i]+1j*a[:, 2+2*i] if r['kind']=='ac' else a[:, 1+i]
              for i in range(len(r['vectors']))]
    return a[:, 0], dict(zip(r['vectors'], values))


def calculate():
    config, headroom = inputs()
    result = dict(schema='learning.design.v1', role=config['claim_role'], feasibility=[], resistor_designs={}, headroom={})
    for key, case in config['cases'].items():
        path = ROOT/'data'/case['array_file']
        if key.startswith('feasibility-'):
            for index in range(3):
                _, v = columns(case, index, path)
                current = float(-v['i(Vout)'][0]); gm = float(v['@m.xn.mapm045_vtg_core[gm]'][0]); gds = float(v['@m.xn.mapm045_vtg_core[gds]'][0])
                result['feasibility'].append(dict(source_case=key, original_analysis_index=index,
                    l_um=case['devices'][0]['l_m']*1e6, w_um=case['devices'][0]['w_m']*1e6,
                    gate_v=float(v['v(gate)'][0]), out_v=float(v['v(out)'][0]), source_v=0., body_v=0.,
                    id_a=current, gm_s=gm, gds_s=gds, gm_over_id_per_v=gm/current,
                    estimated_total_width_32p5u_um=2.*32.5e-6/current,
                    optimistic_gain_32p5u=-(gm/current*32.5e-6)/(85e-6+gds/current*32.5e-6)))
        elif key.startswith(('R37-', 'R75-')):
            _, values = columns(case, 0, path); _, ac = columns(case, 1, path)
            op = {k: float(v[0]) for k, v in values.items()}
            m = case['metadata']; p = m['parameters']
            current = sum(op[f'i({sensors[0]})'] for sensors in m['sensors'].values())
            derivatives = {field: sum(op[f'@m.{name}.mapm045_vtg_core[{field}]'] for name in m['sensors'])
                           for field in ['gm', 'gds', 'gmbs']}
            d = 1+p['rs']*sum(derivatives.values())
            load_g = 1/p['rd']+1/p['rload']
            kcl = {name: sum(op[f'i({sensor})'] for sensor in sensors) for name, sensors in m['sensors'].items()}
            kcl['out'] = current+op['v(out)']/p['rload']+(op['v(out)']-op['v(vdd)'])/p['rd']
            assert max(abs(x) for x in kcl.values()) < 1e-11
            result['resistor_designs'][key] = dict(parameters=p, id_a=current,
                output_v=op['v(out)'], source_v=op['v(source)'], gate_v=op['v(gate)'],
                vgs_v=op['v(gate)']-op['v(source)'], vds_v=op['v(out)']-op['v(source)'],
                total_positive_port_current_a=max(0., -op['i(Vdd)'])+max(0., -op['i(Vin)']),
                derivatives_s=derivatives, degeneration_factor=d, effective_gm_s=derivatives['gm']/d,
                retrospective_equation_gain=-derivatives['gm']/(d*load_g+derivatives['gds']),
                observed_gain_real_1khz=float((ac['v(out)']/ac['v(in)'])[0].real),
                total_width_um=p['nw']*p['nn'], drawn_mos_area_um2=m['resources']['drawn_mos_area_um2'],
                kcl_residuals_a=kcl)
        else:
            axis, v = columns(case, 0, path)
            current = -v['i(Vout)']; error = current/30e-6-1
            within = np.abs(error) <= .01
            result['mirror'] = dict(profile_tier=case['profile_tier'], original_realization_id=case['original_realization_id'],
                definition='Iout=-i(Vout); ratio error=Iout/(30 uA)-1. Retrospective 1% grid illustration; no lower-voltage knee sampled.',
                voltage_v=axis.tolist(), reference_voltage_v=v['v(ref)'].tolist(), output_current_a=current.tolist(),
                fractional_ratio_error=error.tolist(), observed_grid_within_one_percent_v=axis[within].tolist(),
                tested_voltage_interval_v=[float(axis[0]), float(axis[-1])])
    f = next(x for x in result['feasibility'] if abs(x['l_um']-.4)<1e-12 and x['original_analysis_index']==0)
    result['first_estimates'] = {}
    for name, budget in [('R37', 37.5e-6), ('R75', 75e-6)]:
        target = budget-.5/100e3; rd = (1.-.5)/budget
        result['first_estimates'][name] = dict(total_current_a=budget, external_load_current_a=.5/100e3,
            bank_current_a=target, rd_ohm=rd, load_conductance_s=1/rd+1/100e3,
            effective_gm_ignoring_gds_s=6*(1/rd+1/100e3),
            grounded_source_width_estimate_um=2*target/f['id_a'],
            minimum_units_at_4um_width=int(np.ceil(2*target/f['id_a']/4.)))
    for key, case in headroom['cases'].items():
        path = ROOT/'data'/('headroom-'+case['array_file'])
        _, v = columns(case, 0, path)
        result['headroom'][key] = dict(op={n: float(x[0]) for n, x in v.items()},
                                      original_run_id=case['run_id'])
        if key == 'B':
            axis, sweep = columns(case, 1, path)
            result['headroom'][key]['selected_common_mode'] = []
            for wanted in [.4, .5, .6, .9]:
                i = int(np.argmin(abs(axis-wanted))); assert abs(axis[i]-wanted)<1e-12
                result['headroom'][key]['selected_common_mode'].append(dict(
                    common_mode_v=float(axis[i]), tail_v=float(sweep['v(tail)'][i]),
                    bias_v=float(sweep['v(bias)'][i]), output_v=float(sweep['v(dl)'][i]),
                    tail_current_a=float(sweep['i(VtMtd)'][i]),
                    input_vgs_v=float(sweep['v(inp)'][i]-sweep['v(tail)'][i]),
                    input_vds_v=float(sweep['v(dl)'][i]-sweep['v(tail)'][i])))
    return result


def main():
    p = argparse.ArgumentParser(); p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    assert not a.output.exists(), 'Keep prior analysis; choose a new directory'
    result = calculate(); a.output.mkdir(parents=True)
    (a.output/'design-summary.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    rows = result['feasibility']
    with (a.output/'feasibility.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    m = result['mirror']
    with (a.output/'mirror.csv').open('w') as f:
        writer = csv.writer(f, lineterminator='\n'); writer.writerow(['output_v', 'reference_v', 'output_current_a', 'fractional_ratio_error'])
        writer.writerows(zip(m['voltage_v'], m['reference_voltage_v'], m['output_current_a'], m['fractional_ratio_error']))
    print(json.dumps(dict(status='SAVED_DESIGN_ANALYSIS_COMPLETE', feasibility_points=len(rows), resistor_cases=len(result['resistor_designs']), mirror_points=len(m['voltage_v']), output=str(a.output))))


if __name__ == '__main__':
    main()
