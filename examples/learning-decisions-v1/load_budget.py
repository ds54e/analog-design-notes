"""MIT: two predeclared initial-D1 capacitor/deadline observations.

Reuses the verbatim earlier measurement functions. The added deadline is an
illustrative design requirement, not a change to historical D1 qualification.
No native acquisition or parameter selection occurs in this module.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from measurement_v1 import measure_case

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def deadline_result(entry, deadline):
    if entry['status'] == 'RIGHT_CENSORED':
        return 'NOT_DEMONSTRATED'
    bracket = entry['bracket_s']
    if bracket is None:
        return 'NOT_DEMONSTRATED'
    if bracket[1] <= deadline:
        return 'MEETS'
    if bracket[0] > deadline:
        return 'MISSES'
    return 'UNRESOLVED_AT_TIME_RESOLUTION'


def calculate():
    data = ROOT / 'data'; plan = read(ROOT / 'load-budget-plan.json')
    assert sha(ROOT / 'measurement_v1.py') == plan['measurement_dependency']['selected_module_sha256']
    cfg = read(data / 'interface-inputs.json')
    assert sha(ROOT / 'load-budget-plan.json') == cfg['plan_sha256']
    assert cfg['model_hashes'] == plan['model_hashes'] and cfg['apm_commit'] == plan['apm_commit']
    assert set(cfg['cases']) == {c['id'] for c in plan['cases']} == {'D1-budget-12p', 'D1-budget-15p'}
    for name, expected in read(data / 'interface-manifest.json')['files'].items():
        assert sha(data / name) == expected, name
    training_path = data / 'training-interface-summary.json'
    assert sha(training_path) == plan['known_evidence']['training_summary_sha256']
    training = read(training_path)['cases']['D1-5p']
    result = dict(schema='learning.initial-d1-load-budget.measure.v1',
        original_measurement_version=plan['measurements']['version'], plan_identity=plan['identity'],
        role='Prospective condition check on a known nominal physical core; not fresh physical samples or population confirmation.',
        deadline_s=plan['deadline']['seconds'], cases={}, prospective_estimates=plan['prospective_estimates'])
    for case in plan['cases']:
        name = case['id']; info = cfg['cases'][name]
        assert info['case'] == case and info['physical_devices'] == plan['physical_devices']
        assert info['profile_tier'] == 'ARTIFICIAL' and info['target_started_utc'] > plan['created_utc']
        assert sha(data / info['array_file']) == info['array_sha256']
        assert sha(ROOT / 'circuits' / info['circuit_file']) == info['circuit_sha256']
        with np.load(data / info['array_file'], allow_pickle=False) as z:
            arrays = {k: z[k] for k in z.files}
        measured = measure_case(info, arrays, plan)
        delta = max(abs(a['output_v'] - b['output_v']) for a, b in zip(measured['dc'], training['dc']))
        assert delta <= plan['measurements']['dc_load_difference_tolerance_v']
        measured['maximum_dc_difference_from_known_5p_v'] = delta
        forecast = plan['prospective_estimates']['cases'][name]
        edges = {}
        for edge in ['rise', 'fall']:
            entry = measured['transient']['entries'][edge]
            prediction = forecast['scaled_entry_upper_estimate_s'][edge]
            edges[edge] = dict(status=deadline_result(entry, result['deadline_s']),
                observed_entry=entry, prior_scaled_entry_estimate_s=prediction,
                observed_upper_minus_prior_estimate_s=entry['bracket_s'][1] - prediction if entry['bracket_s'] else None)
        status = 'MEETS' if all(x['status'] == 'MEETS' for x in edges.values()) else (
            'MISSES' if any(x['status'] == 'MISSES' for x in edges.values()) else 'NOT_DEMONSTRATED_OR_UNRESOLVED')
        result['cases'][name] = dict(measured=measured, edges=edges, deadline_result=status,
            predeclared_deadline_forecast=forecast['predicted_deadline_result'],
            agrees_with_forecast=status == forecast['predicted_deadline_result'],
            target_started_utc=info['target_started_utc'], saved_realization_id=info['original_realization_id'],
            original_attempt_sha256=info['original_attempt_sha256'])
    result['limits'] = [
        'A 500-ns deadline to actual OP endpoints does not remove the finite DC tracking error.',
        'Two discrete capacitor conditions on the original trajectory do not establish a continuous allowable-C boundary.',
        'The capacitor is a fixed diagnostic load, not a transistor pass-gate model or complete LDO.',
        'No new population, arbitrary initial charge, other supply/temperature or history is assessed.',
    ]
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    result = calculate()
    (a.output / 'load-budget-summary.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    print(json.dumps(dict(status='DECLARED_CONDITIONS_ANALYZED',
        cases={k: v['deadline_result'] for k, v in result['cases'].items()}, output=str(a.output))))


if __name__ == '__main__':
    main()
