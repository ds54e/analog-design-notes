"""MIT: regenerate the finite allocation confirmation statistics from public CSV.

Intervals describe uncertainty conditional on this model and frozen pairing.
This script neither samples a process nor runs SPICE.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import t


def read(data, name):
    with (data/(name+'.csv')).open() as f:
        return list(csv.DictReader(f))


def calculate(data):
    config = json.loads((data/'configuration.json').read_text()); plan = config['plan']
    rows = read(data, 'confirmation'); expected = {(n, i, c['id']) for n in config['candidates'] for i in config['indices'] for c in config['conditions']}
    keys = [(r['candidate'], int(r['index']), r['condition']) for r in rows]
    assert len(keys) == len(set(keys)) and set(keys) == expected
    summaries = {}; contrasts = []; pairs = plan['contrasts']; family = len(pairs)*len(plan['conditions'])
    for condition in plan['conditions']:
        c = condition['id']; groups = {}
        for name in plan['candidates']:
            assigned = [r for r in rows if r['condition'] == c and r['candidate'] == name]
            valid = [r for r in assigned if r['numerical_usable'] == 'True' and r['terminal_supported'] == 'True' and r['out_v']]
            groups[name] = {int(r['index']): float(r['error_v']) for r in valid}
            v = np.array([float(r['error_v']) for r in valid]); pred = np.array([float(r['predicted_out_v'])-float(r['out_v']) for r in valid])
            # Recompute the uniform pass independently of its saved Boolean.
            for r in valid:
                passed = abs(float(r['error_v'])) <= .025 and abs(float(r['gain'])/-6-1) <= .1 and float(r['band_min_magnitude']) >= 5.4 and float(r['band_max_magnitude']) <= 6.6 and bool(r['bandwidth_hz']) and float(r['bandwidth_hz']) >= 1e6 and float(r['current_a']) <= 80e-6
                assert bool(passed) == (r['dc_ac_current_pass'] == 'True')
            summaries[name+'/'+c] = dict(assigned=len(assigned), usable_supported=len(valid), missing_or_unsupported=len(assigned)-len(valid),
                mean_error_v=float(v.mean()), sd_error_v=float(v.std(ddof=1)), mae_v=float(np.mean(abs(v))), rms_error_v=float(np.sqrt(np.mean(v*v))),
                max_absolute_error_v=float(max(abs(v))), dc_ac_current_pass=sum(r['dc_ac_current_pass'] == 'True' for r in valid),
                prediction_rmse_v=float(np.sqrt(np.mean(pred*pred))), prediction_max_absolute_error_v=float(max(abs(pred))),
                sine_assessed=sum(r['sine_assessed'] == 'True' for r in valid), sine_all_assessed_pass=sum(r['sine_all_assessed_pass'] == 'True' for r in valid))
        for pair in pairs:
            left, right = pair['left'], pair['right']; ids = sorted(groups[left].keys() & groups[right].keys())
            delta = np.array([abs(groups[left][i])-abs(groups[right][i]) for i in ids]); n = len(delta)
            sd = float(delta.std(ddof=1)); mean = float(delta.mean()); half = float(t.ppf(1-.05/(2*family), n-1)*sd/np.sqrt(n))
            contrasts.append(dict(condition=c, left=left, right=right, n_pairs=n, mean_mae_reduction_v=mean,
                paired_sd_v=sd, interval_v=[mean-half, mean+half], family_size=family,
                meaningful_reduction_v=pair.get('meaningful_reduction_v'), full_sampled_model_comparison_allowed=n >= .99*len(plan['indices']),
                interpretation='Positive means right has lower MAE. Frozen Bonferroni paired t intervals across all stated contrasts and conditions; model-conditional approximate uncertainty.'))
    return dict(schema='ads.allocation-statistics.v1', role=plan['role'], summaries=summaries, contrasts=contrasts)


def verify_source_budget(data):
    computed = {}
    for row in read(data, 'source-transforms'):
        j = np.array([[float(row['j00']), float(row['j01'])], [float(row['j10']), float(row['j11'])]])
        s = np.array([float(row['sensitivity_delvto']), float(row['sensitivity_ln_mulu0'])])
        sigma = np.diag([float(row['sigma_vth']), float(row['sigma_beta'])]); u = s @ np.linalg.inv(j) @ sigma
        key = (row['candidate'], row['uid'].split('/')[-2]); computed[key] = computed.get(key, 0.) + float(sum(u*u))
    for row in read(data, 'source-budget'):
        assert np.isclose(computed[row['candidate'], row['role']], float(row['variance_v2']), rtol=1e-13, atol=1e-17)
    return computed


def main():
    p = argparse.ArgumentParser(); p.add_argument('--data', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    result = calculate(a.data); expected = json.loads((a.data/'statistics.json').read_text())
    assert result == expected, 'Public CSV statistics differ from audited private calculation'
    sources = verify_source_budget(a.data)
    with a.output.open('x') as f:
        f.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps(dict(status='PASS', candidate_conditions=len(result['summaries']), paired_contrasts=len(result['contrasts']), source_role_terms=len(sources))))


if __name__ == '__main__':
    main()
