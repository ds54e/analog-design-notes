"""MIT: retrospective divider KCL/error accounting and new-OP local model."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import analyze
import prior

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-zip', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    plan = analyze.inputs(args.source_zip)
    name = 'D1-half-divider'
    case = plan['cases'][name]
    index = analyze.read(ROOT/'data/index.json')
    entry = index['cases'][name]
    assert index['plan_sha256'] == analyze.sha(ROOT/'plan.json')
    assert analyze.sha(ROOT/'data'/entry['file']) == entry['sha256']
    with np.load(ROOT/'data'/entry['file'], allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    measured = analyze.measure_case(name, case, arrays, plan)
    assert measured == analyze.read(ROOT/'data/summary.json')['cases'][name]
    op, fixed = measured['op'], plan['fixed']
    matrix, residual, drive = prior.stamp(op, case['connections'], fixed)
    h = float(np.linalg.solve(matrix, drive)[prior.NODES.index('out')])
    output_drive = [1. if n == 'out' else 0. for n in prior.NODES]
    z = float(np.linalg.solve(matrix, output_drive)[prior.NODES.index('out')])
    reduced = prior.explicit_reduction(op, fixed)
    assert np.isclose(h, reduced['gain'], rtol=1e-12)
    assert np.isclose(z, reduced['output_transimpedance_ohm'], rtol=1e-12)
    assert max(abs(residual)) <= plan['measurements']['kcl_max_a']
    pole = 1/(2*np.pi*fixed['cload_f']*z)
    bw = pole*np.sqrt(10**.3*(1+(1000/pole)**2)-1)
    beta, parallel = measured['divider_beta'], measured['divider_parallel_ohm']
    terms = dict(positive_input_source_drop=-fixed['rsource_ohm']*op['i(VtM1g)']/beta,
        finite_input_difference=-measured['input_difference_v']/beta,
        negative_input_gate_loading=parallel*op['i(VtM2g)']/beta)
    assert abs(sum(terms.values())-measured['error_v']) <= 1e-12
    f, ac = analyze.columns(case['analyses'][1], arrays['analysis1'])
    # This is a measured closed-circuit transfer, not an independent divider model.
    beta_effective = ac['v(inn)'][0]/ac['v(out)'][0]
    beta_kcl = beta-parallel*ac['i(VtM2g)'][0]/ac['v(out)'][0]
    assert abs(beta_effective-beta_kcl) < 1e-12
    input_transfer = ac['v(inp)'][0]/ac['v(source)'][0]
    input_kcl = 1-fixed['rsource_ohm']*ac['i(VtM1g)'][0]/ac['v(source)'][0]
    assert abs(input_transfer-input_kcl) < 1e-12
    # Label this explicitly: substituting observed gate transfers explains the
    # omitted gate-current path but does not produce a new independent forecast.
    open_z = z*reduced['feedback_denominator']
    h_observed_gates = reduced['a_plus']*input_transfer/(1+
        2j*np.pi*f[0]*fixed['cload_f']*open_z-reduced['a_minus']*beta_effective)
    observed_h = complex(*measured['gain_one_khz'])
    old = analyze.read(ROOT/'data/prior.json')
    result = dict(schema='learning.initial-d1-feedback-divider.network.v1',
        role='Retrospective new-OP local explanation and exact stationary error accounting; original prior remains frozen.',
        plan_sha256=analyze.sha(ROOT/'plan.json'), observation_sha256=entry['sha256'],
        old_prior_sha256=analyze.sha(ROOT/'data/prior.json'),
        nodes=prior.NODES, conductance_matrix_s=matrix.tolist(),
        full_native_operating_point_residual_a=residual.tolist(),
        reduced_network=reduced, full_stamp_gain=h, full_stamp_output_transimpedance_ohm=z,
        single_external_output_pole_hz=float(pole), single_pole_minus3db_hz=float(bw),
        gain_relative_error=h/observed_h.real-1,
        bandwidth_relative_error=float(bw/measured['bandwidth']['hz']-1),
        exact_stationary_error_accounting_v=terms,
        exact_stationary_error_sum_v=sum(terms.values()),
        accounting_limit='These correlated quantities belong to this one operating point. They are not independent random errors or isolated causal interventions.',
        observed_gate_path=dict(frequency_hz=float(f[0]),
            external_output_capacitance_f=fixed['cload_f'],
            open_output_transimpedance_ohm=open_z,
            sensing_transfer=[float(beta_effective.real), float(beta_effective.imag)],
            input_transfer=[float(input_transfer.real), float(input_transfer.imag)],
            reconstructed_gain=[float(h_observed_gates.real), float(h_observed_gates.imag)],
            native_gain=[float(observed_h.real), float(observed_h.imag)],
            real_gain_relative_error=float(h_observed_gates.real/observed_h.real-1),
            role='Uses already measured gate-current/voltage transfers at1kHz, DC MOS derivatives and the explicit external output capacitor. Explains a neglected path; not a prospective or independent transfer prediction.'),
        costs=dict(added_resistance_ohm=measured['added_divider_resistance_ohm'],
            divider_power_w=measured['stationary_power_w']['resistors']['Rupper']+measured['stationary_power_w']['resistors']['Rlower'],
            supply_delivered_a=measured['supply_delivered_a'],
            total_source_delivery_w=measured['stationary_power_w']['total_source_delivery'],
            drawn_mos_area_um2=sum(d['w_m']*d['l_m']*1e12 for d in plan['physical_units'])),
        old_and_new_derivatives={str(i):dict(old=prior.derivatives(old['known']['observed_op'],i),
            new=prior.derivatives(op,i), order=['gm','gmb','gds'], units='S') for i in range(1,7)},
        limits=['Finite gain-of-two connection changes function and adds divider loading; do not rank against unity feedback as the same job.',
            'The unfitted new-OP model still omits gate/body-current derivatives and intrinsic capacitance; retain its residuals.',
            'No new transient, noise, return-ratio, current-limited supply or LDO result.'])
    args.output.mkdir(parents=True)
    (args.output/'network.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    with (args.output/'error-accounting.csv').open('x', newline='') as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(['term','output_error_v'])
        writer.writerows(terms.items())
        writer.writerow(['sum',sum(terms.values())])
    print(json.dumps({k:result[k] for k in ['gain_relative_error','bandwidth_relative_error',
        'exact_stationary_error_accounting_v','observed_gate_path','costs']}, indent=2))


if __name__ == '__main__':
    main()
