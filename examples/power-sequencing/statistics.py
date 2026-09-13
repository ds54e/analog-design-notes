"""MIT: selected deterministic summaries; no sampling or statistical inference."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def summarize(data):
    manifest = json.loads((data/'manifest.json').read_text())
    for name, digest in manifest['files'].items():
        assert hashlib.sha256((data/name).read_bytes()).hexdigest() == digest, name
    def read(name):
        with (data/(name+'.csv')).open() as f: return list(csv.DictReader(f))
    static = read('static'); histories = read('histories'); plateaus = read('plateaus'); paths = read('input-paths')
    observations = json.loads((data/'observations.json').read_text())['results']
    assert len(static) == 12 and len(histories) == 47 and len(observations) == 59
    states = {}
    for r in static:
        assert abs(float(r['input_delivered_a'])-float(r['input_series_term_a'])-float(r['input_gate_term_a'])) < 1e-10
        assert abs(float(r['signed_port_power_w'])-float(r['vdd_v'])*float(r['supply_delivered_a'])-float(r['vin_v'])*float(r['input_delivered_a'])) < 1e-15
        states[r['candidate']+'/'+r['state']] = {k: float(v) for k, v in r.items() if k not in ['candidate', 'state']}
    statuses = {s: sum(r['status'] == s for r in plateaus) for s in sorted({r['status'] for r in plateaus})}
    assert not any(r['status'] not in ['OBSERVED_FINAL_ENTRY', 'IN_BAND_AT_FIRST_OBSERVATION'] for r in plateaus)
    selected = {}
    for r in histories:
        if r['kind'] == 'startup' and r['refined'] == 'False' and float(r['capacitance_f']) == 8e-12 and float(r['ramp_time_s']) == 1e-8:
            pt = next(p for p in plateaus if p['case_id'] == r['case_id'])
            selected[r['case_id']] = dict(settling_s=float(pt['resolved_entry_time_s']),
                input_delivery_peak_a=float(r['Vin_positive_delivery_peak_a']), input_absorption_peak_a=float(r['Vin_absorption_peak_a']),
                supply_delivery_peak_a=float(r['Vdd_positive_delivery_peak_a']), output_max_v=float(r['out_max_v']))
    refined = [r['refinement'] for r in observations if 'refinement' in r]; assert len(refined) == 4 and all(r['status'] == 'PASS' for r in refined)
    return dict(schema='adn.power-sequence-summary.v1', scope='Finite nominal deterministic observations; no population or all-history inference.',
        static_states=states, fast_8pf=selected, plateau_status_counts=statuses,
        maximum_input_identity_residual_a=max(float(r['maximum_input_kcl_identity_residual_a']) for r in paths),
        maximum_refinement_waveform_difference_v=max(r['maximum_piecewise_linear_waveform_difference_v'] for r in refined),
        maximum_plateau_endpoint_error_v=max(abs(float(r['endpoint_error_v'])) for r in plateaus),
        analytic_control=observations[0]['measurements']['transient']['analytic_control'])


def main():
    p = argparse.ArgumentParser(); p.add_argument('--data', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); result = summarize(a.data)
    with a.output.open('x') as f: f.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print('Verified compact data hashes and deterministic summaries; no SPICE.')


if __name__ == '__main__': main()
