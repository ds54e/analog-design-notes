"""MIT: finite, retrospective supply-path and fixed-target loss calculations.

Reads fifteen selected members of one immutable public ZIP. No simulator,
private repository, historical runner or automatic acquisition is imported.
"""
import argparse
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CANDIDATES = ['A75', 'R75', 'F75S', 'F75L']
CONDITIONS = ['nominal', 'supply_097']


def read_source(archive):
    source = json.loads((ROOT / 'data/source.json').read_text())
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == source['archive_sha256']
    with zipfile.ZipFile(archive) as z:
        selected = {name: z.read('example/' + name) for name in source['files']}
    for name, contents in selected.items():
        assert hashlib.sha256(contents).hexdigest() == source['files'][name], name
    manifest = json.loads(selected['data/manifest.json'])
    for name in selected:
        if name.startswith('data/') and name != 'data/manifest.json':
            assert source['files'][name] == manifest['files'][name[5:]], name
    return source, selected


def csv_rows(selected, name):
    return list(csv.DictReader(io.StringIO(selected['data/' + name].decode())))


def local_paths(parameters, op):
    """Separate a four-node linear KCL solve from the reduced expression."""
    p = parameters

    def total(prefix, field):
        count = p.get({'n': 'nn', 'p': 'np', 'r': 'nr'}[prefix], 0)
        return sum(op[f'@m.x{prefix}{i}.mapm045_vtg_core[{field}]'] for i in range(count))

    gm, gb, gd = (total('n', x) for x in ['gm', 'gmbs', 'gds'])
    gp, gdp = total('p', 'gm'), total('p', 'gds')
    gr = total('r', 'gm') + total('r', 'gds')
    rs = p['rs']; rt = p['rsrc'] + p.get('rin', 0.)
    gf = 1 / p['rf'] if 'rf' in p else 0.
    gload = 1 / p['rload']
    drain_g = 1 / p['rd'] if p['load'] == 'resistor' else 0.
    bias_g = 1 / p['rbias'] if gr else 0.

    # Unknowns: source, gate, out, pbias (the last is a dummy zero for R75).
    # Excitations: ground-referred Vin, VDD, output current injected into out,
    # and positive Vbiasdiag = pbias - pdrive (lower PMOS gate).
    matrix = np.array([
        [1 / rs + gm + gb + gd, -gm, -gd, 0.],
        [0., 1 / rt + gf, -gf, 0.],
        [-gm - gb - gd, gm - gf, gd + gload + gf + gdp + drain_g, gp],
        [0., 0., 0., gr + bias_g if gr else 1.],
    ])
    rhs = np.array([
        [0., 0., 0., 0.],
        [1 / rt, 0., 0., 0.],
        [0., gp + gdp + drain_g, 1., gp],
        [0., gr, 0., 0.],
    ])
    solved = np.linalg.solve(matrix, rhs)
    np.testing.assert_allclose(matrix @ solved, rhs, rtol=2e-12, atol=1e-15)

    degeneration = 1 + rs * (gm + gb + gd)
    effective_gm = gm / degeneration
    go = gload + gd / degeneration + gdp + drain_g
    if gf:
        denominator = go + (1 + effective_gm * rt) / (p['rf'] + rt)
        gain = -(effective_gm * p['rf'] - 1) / (
            go * (p['rf'] + rt) + 1 + effective_gm * rt)
    else:
        denominator = go
        gain = -effective_gm / go
    zo = 1 / denominator
    vsg_fraction = bias_g / (gr + bias_g) if gr else 0.
    reference_injection = gp * vsg_fraction
    direct_injection = gdp + drain_g
    line = zo * (reference_injection + direct_injection)
    reduced = np.array([gain, line, zo, gp * zo])
    np.testing.assert_allclose(solved[2], reduced, rtol=2e-12, atol=1e-12)
    if gr:
        np.testing.assert_allclose(solved[3, 1], 1 - vsg_fraction, rtol=2e-12)
    return dict(reference_units=p.get('nr', 0), rbias_ohm=p.get('rbias'),
        reference_conductance_s=gr, rbias_times_reference_conductance=p.get('rbias', 0) * gr,
        pbias_tracking_v_per_v=1 - vsg_fraction if gr else None,
        load_vsg_supply_fraction=vsg_fraction if gr else None,
        load_gm_s=gp, load_gds_s=gdp, effective_input_gm_s=effective_gm,
        injection_ohm_estimate=zo, input_v_per_v_estimate=gain,
        bias_gate_drop_v_per_v_estimate=gp * zo if gr else None,
        direct_supply_injection_a_per_v=direct_injection,
        reference_supply_injection_a_per_v=reference_injection,
        total_supply_injection_a_per_v=reference_injection + direct_injection,
        direct_supply_path_v_per_v=zo * direct_injection,
        reference_supply_path_v_per_v=zo * reference_injection,
        supply_v_per_v_estimate=line,
        matrix_reduced_max_relative_difference=float(max(abs(solved[2] - reduced) /
            np.where(abs(reduced) > 0, abs(reduced), 1.))))


def calculate(archive):
    source, selected = read_source(archive)
    config = json.loads(selected['data/configuration.json'])
    contract = json.loads(selected['contract-v1.json'])
    plan = json.loads(selected['confirmation-plan.json'])
    assert contract['candidates'] == CANDIDATES
    assert plan['indices'] == list(range(2000, 2521))
    assert [c['id'] for c in plan['conditions']] == CONDITIONS
    original = json.loads(selected['data/statistics.json'])
    nominal = {r['candidate']: r for r in csv_rows(selected, 'nominal.csv')}
    sensitivity = {r['candidate']: r for r in csv_rows(selected, 'sensitivities.csv')}
    op = {c: {} for c in CANDIDATES}
    for r in csv_rows(selected, 'nominal-op.csv'):
        assert r['vector'] not in op[r['candidate']]
        op[r['candidate']][r['vector']] = float(r['value'])
    condition_rows = {(r['candidate'], r['condition']): r for r in csv_rows(selected, 'conditions.csv')}
    paths = []
    for c in CANDIDATES:
        model = config['candidates'][c]
        assert hashlib.sha256(selected['circuits/' + c + '.cir']).hexdigest() == model['circuit_sha256']
        row = dict(candidate=c, **local_paths(model['parameters'], op[c]))
        for name in ['injection_ohm', 'input_v_per_v', 'supply_v_per_v', 'bias_gate_drop_v_per_v']:
            measured = float(sensitivity[c][name]) if sensitivity[c][name] else None
            row[name + '_observed'] = measured
            row[name + '_relative_discrepancy'] = row[name + '_estimate'] / measured - 1 if measured else None
        row['nominal_out_v'] = float(nominal[c]['out_v'])
        for key in ['area_um2', 'current_a', 'internal_resistance_ohm', 'input_noise_v']:
            row[key] = float(nominal[c][key])
        row['supply_097_local_out_estimate_v'] = row['nominal_out_v'] - .03 * row['supply_v_per_v_estimate']
        observed = condition_rows.get((c, 'supply_097'))
        row['supply_097_nominal_device_observed_v'] = float(observed['out_v']) if observed else None
        row['supply_097_extrapolation_residual_v'] = (
            float(observed['out_v']) - row['supply_097_local_out_estimate_v'] if observed else None)
        paths.append(row)

    confirmation = csv_rows(selected, 'confirmation.csv')
    assert len(confirmation) == len(CANDIDATES) * len(CONDITIONS) * len(plan['indices'])
    keys = [(r['candidate'], r['condition'], int(r['index'])) for r in confirmation]
    assert len(set(keys)) == len(keys)
    expected = {(c, condition, i) for c in CANDIDATES for condition in CONDITIONS for i in plan['indices']}
    assert set(keys) == expected
    limits = contract['limits']; job = contract['function']
    summaries = []; errors = []
    for c in CANDIDATES:
        for condition in CONDITIONS:
            rows = sorted([r for r in confirmation if r['candidate'] == c and r['condition'] == condition], key=lambda r: int(r['index']))
            assert len(rows) == len(plan['indices'])
            assert all(r['numerical_usable'] == 'True' and r['terminal_supported'] == 'True' for r in rows)
            e = np.array([float(r['error_v']) for r in rows])
            np.testing.assert_allclose(e, [float(r['out_v']) - job['output_center_v'] for r in rows], rtol=0, atol=1e-15)
            counts = dict(center_only_pass=0, five_checks_pass=0, below_center_band=0,
                above_center_band=0, gain_fail=0, band_fail=0, bandwidth_fail=0, current_fail=0)
            for r, error in zip(rows, e):
                tests = dict(center=abs(error) <= limits['center_abs_error_v'],
                    gain=abs(float(r['gain']) / job['signed_gain'] - 1) <= limits['gain_relative_error'],
                    band=float(r['band_min_magnitude']) >= abs(job['signed_gain']) * (1 - limits['gain_relative_error']) and
                         float(r['band_max_magnitude']) <= abs(job['signed_gain']) * (1 + limits['gain_relative_error']),
                    bandwidth=bool(r['bandwidth_hz']) and float(r['bandwidth_hz']) >= limits['bandwidth_min_hz'],
                    current=float(r['current_a']) <= limits['current_max_a'])
                passed = all(tests.values())
                assert passed == (r['dc_ac_current_pass'] == 'True')
                counts['center_only_pass'] += int(tests['center'])
                counts['five_checks_pass'] += int(passed)
                counts['below_center_band'] += int(error < -limits['center_abs_error_v'])
                counts['above_center_band'] += int(error > limits['center_abs_error_v'])
                for name in ['gain', 'band', 'bandwidth', 'current']:
                    counts[name + '_fail'] += int(not tests[name])
                errors.append(dict(candidate=c, condition=condition, index=int(r['index']),
                    error_v=float(error), center_pass=bool(tests['center']), five_checks_pass=passed,
                    realization_id=r['realization_id'], original_attempt_sha256=r['attempt_sha256']))
            mean = float(np.mean(e)); sd = float(np.std(e, ddof=1)); mse = float(np.mean(e * e))
            mean_square = mean * mean; spread_square = (len(e) - 1) / len(e) * sd * sd
            np.testing.assert_allclose(mse, mean_square + spread_square, rtol=2e-13, atol=1e-18)
            row = dict(candidate=c, condition=condition, n=len(e), mean_error_v=mean, sd_error_v=sd,
                rms_error_v=float(np.sqrt(mse)), mae_v=float(np.mean(abs(e))),
                mean_square_v2=mean_square, spread_square_v2=spread_square,
                mean_square_fraction_of_mse=mean_square / mse,
                rms_identity_residual_v2=mse - mean_square - spread_square, **counts)
            old = original['summaries'][c + '/' + condition]
            for key in ['mean_error_v', 'sd_error_v', 'rms_error_v', 'mae_v']:
                np.testing.assert_allclose(row[key], old[key], rtol=2e-13, atol=1e-15)
            assert row['five_checks_pass'] == old['dc_ac_current_pass']
            summaries.append(row)
    paired = []
    for c in CANDIDATES:
        a = np.array([r['error_v'] for r in errors if r['candidate'] == c and r['condition'] == 'nominal'])
        b = np.array([r['error_v'] for r in errors if r['candidate'] == c and r['condition'] == 'supply_097'])
        shift = b - a
        paired.append(dict(candidate=c, paired_coordinates=len(a), mean_supply_error_change_v=float(np.mean(shift)),
            sd_supply_error_change_v=float(np.std(shift, ddof=1)),
            nominal_device_local_change_estimate_v=-.03 * next(r['supply_v_per_v_estimate'] for r in paths if r['candidate'] == c),
            role='Known paired observations; the nominal-device derivative is not a complete per-realization predictor.'))
    return dict(schema='learning.supply-loss.v1', role=source['role'], source_archive_sha256=source['archive_sha256'],
        selected_source_members=len(source['files']), paths=paths, loss=summaries, paired_supply_changes=paired,
        center_band_v=limits['center_abs_error_v'], errors=errors,
        limitations=['Known OP nodal approximation neglects leakage derivatives and capacitance; it is not a new confirmation.',
            'Original local supply slopes use backward perturbations at the 1-V model boundary; the reduced calculation is a derivative.',
            'Retain original mean, sample SD, RMS, MAE, center and five-check definitions; no re-centering or statistical refit.',
            'Fixed modeled physical coordinates are not manufacturing yield, and changed inventories are not identical complete physical circuits.'])


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--feedback-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    result = calculate(a.feedback_zip)
    for key, name in [('paths', 'supply-paths.csv'), ('loss', 'loss-summary.csv'),
                      ('paired_supply_changes', 'paired-supply.csv'), ('errors', 'known-errors.csv')]:
        write_csv(a.output / name, result[key])
    summary = {k: v for k, v in result.items() if k != 'errors'}
    summary['known_condition_rows'] = len(result['errors'])
    (a.output / 'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(status='PASS', candidates=len(result['paths']), known_conditions=len(result['errors']), output=str(a.output))))


if __name__ == '__main__':
    main()
