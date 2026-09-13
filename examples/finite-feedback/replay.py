"""MIT: replay saved physical OP/AC or preselected 20-mV sine observations.

Uses the pinned public APM recipe API. No sampling, remapping or trimming.
"""
import argparse
import copy
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def periodic(t, y):
    """Original eight-cycle, ten-harmonic measurement; no fitted tail."""
    import numpy as np
    assert len(t) == len(y) and np.isfinite(t).all() and np.isfinite(y).all()
    assert np.all(np.diff(t) > 0)
    stop = float(t[-1]); start = stop-8/1e4
    assert start >= t[0] and max(np.diff(t[t >= start])) <= 1.1e-7
    grid = start+np.arange(8000)/1e7; values = np.interp(grid, t, y)
    ft = np.fft.rfft(values)/len(values)
    amps = np.array([2*abs(ft[n*8]) for n in range(1, 11)])
    assert amps[0] >= 1e-9
    cycles = values.reshape(8, 1000)
    delta = float(np.max(abs(cycles[-1]-cycles[-2])))
    return dict(fundamental_peak_v=float(amps[0]), thd_fraction=float(np.linalg.norm(amps[1:])/amps[0]),
                stationary=bool(delta < max(1e-7, amps[0]*1e-4)), max_last_cycle_difference_v=delta)


def resistor_checks(values, design):
    """Infer actual resistance from measured port/terminal V/I, then compare."""
    p = design['parameters']; terminals = design['terminals']; sensors = design['sensors']
    def v(node): return 0. if node == '0' else values['v('+node+')']
    def current(name, terminal): return values['i('+sensors[name][terminal]+')']
    source_sum = sum(current(n, 2) for n, nodes in terminals.items() if nodes[2] == 'source')
    drain_sum = sum(current(n, 0) for n, nodes in terminals.items() if nodes[0] == 'out')
    incoming = -values['i(Vin)']; checks = {}
    def check(name, drop, flow, expected):
        assert abs(flow) > 1e-12, (name, 'V/I readback is ill-conditioned')
        actual = drop/flow
        assert abs(actual/expected-1) < 1e-5, (name, actual, expected)
        checks[name] = dict(measured_v_over_i_ohm=actual, declared_ohm=expected)
    check('Rs', v('source'), -source_sum, p['rs'])
    if 'rf' in p:
        gate = sum(current(n, 1) for n, nodes in terminals.items() if nodes[1] == 'gate')
        feedback = gate-incoming
        check('Rsrc', v('in')-v('drive'), incoming, p['rsrc'])
        check('Rin', v('drive')-v('gate'), incoming, p['rin'])
        check('Rf', v('out')-v('gate'), feedback, p['rf'])
    else:
        feedback = 0.
        check('Rsrc', v('in')-v('gate'), incoming, p['rsrc'])
    rd_flow = 0.
    if p.get('rd'):
        rd_flow = -values['i(Vdd)']
        check('Rd', v('vdd')-v('out'), rd_flow, p['rd'])
    if p.get('rbias'):
        bias = -values['i(Vbiasdiag)']
        for name, nodes in terminals.items():
            bias -= sum(current(name, j) for j, node in enumerate(nodes) if node == 'pbias')
        check('Rbias', v('pbias'), bias, p['rbias'])
    check('Rload', v('out'), rd_flow-drain_sum-feedback, p['rload'])
    return checks


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--apm', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--binary', type=Path, default=Path('/usr/local/bin/ngspice'))
    p.add_argument('--candidate', choices=['A75', 'R75', 'F75S', 'F75L'], default='F75L')
    p.add_argument('--index', type=int, default=2000)
    p.add_argument('--condition', choices=['nominal', 'supply_097', 'signal_20m'], default='nominal')
    a = p.parse_args(); here = Path(__file__).resolve().parent; apm = a.apm.resolve(); out = a.output
    if not out.is_absolute() or out.exists() or out.resolve().is_relative_to(apm):
        raise ValueError('Select a new absolute output outside the model checkout')
    cfg = json.loads((here/'data/configuration.json').read_text())
    for name, digest in json.loads((here/'data/manifest.json').read_text())['files'].items():
        assert sha(here/'data'/name) == digest, ('public data changed', name)
    commit = subprocess.check_output(['git', '-C', str(apm), 'rev-parse', 'HEAD'], text=True).strip()
    assert commit == cfg['apm_commit']
    for name, digest in cfg['model_hashes'].items(): assert sha(apm/name) == digest, name
    def read(name):
        with (here/'data'/name).open() as f: return list(csv.DictReader(f))
    d = cfg['candidates'][a.candidate]; signal = a.condition == 'signal_20m'
    recipe = cfg['selected_signal'][a.candidate] if signal else d
    body = here/'circuits'/(a.candidate+('-signal20m' if signal else '')+'.cir')
    assert sha(body) == recipe['circuit_sha256']
    physical = [r for r in read('realizations.csv') if r['candidate'] == a.candidate and int(r['index']) == a.index]
    assert len(physical) == len(d['devices']) and all(r['status'] == 'RESOLVED' for r in physical)
    rows = {r['uid']: r for r in physical}; assert set(rows) == {u['uid'] for u in d['devices']}
    for unit in d['devices']:
        row = rows[unit['uid']]
        assert row['polarity'] == unit['polarity'] and float(row['w_m']) == unit['w_m'] and float(row['l_m']) == unit['l_m']
    main_condition = 'nominal' if signal else a.condition
    expected = next(r for r in read('confirmation.csv') if r['candidate'] == a.candidate and int(r['index']) == a.index and r['condition'] == main_condition)
    if signal:
        assert a.index in cfg['plan']['selected_large_signal']['indices']
        expected_signal = next(r for r in read('selected-signal.csv') if r['candidate'] == a.candidate and int(r['index']) == a.index)
    os.environ['APM_REPO_ROOT'] = str(apm); sys.path.insert(0, str(apm/'src'))
    import numpy as np
    from apm import research, research_spice as spice
    request = dict(schema=research.SCHEMAS['request'], circuit=str(body), devices=d['devices'],
                   other_variation_leaves=[], analyses=copy.deepcopy(recipe['analyses']))
    for analysis in request['analyses']:
        analysis['vectors'] = [v for v in analysis['vectors'] if not v.startswith('@')]
        if not signal:
            condition = next(c for c in cfg['conditions'] if c['id'] == a.condition)
            analysis['set_sources'] = condition.get('sources', {})
        if analysis['kind'] == 'tran':
            # Pinned public API emits tran step stop. The default tmax is
            # min(step, stop/50), equal to the frozen explicit 100-ns tmax.
            assert analysis['step'] == analysis['tmax'] == 1e-7 and analysis['stop'] == .001
    binding = {}; spice.flatten(body, binding)
    for name in spice.MODELS: spice.flatten(apm/name, binding)
    real = research.seal(dict(schema=research.SCHEMAS['realization'], origin='research',
        profile_tier='SOURCE_TRANSFER_HYPOTHESIS', status='RESOLVED', sample_context_id=research.sample_context(request),
        request_id=research.canonical_hash(request), input_binding=binding,
        devices=[{**u, 'raw': [float(rows[u['uid']]['delvto_v']), float(rows[u['uid']]['ln_mulu0'])]} for u in d['devices']],
        public_replay=dict(csv_sha256=sha(here/'data/realizations.csv'), original_realization_id=physical[0]['original_realization_id'],
            seed=cfg['seed'], sample_index=a.index, rule='Exact saved UID/geometry/raw; relocated file bindings, no sampling or trimming.')))
    out.mkdir(parents=True); (out/'relocated-realization.json').write_text(json.dumps(real, indent=2)+'\n')
    report = spice.execute(apm, a.binary, out/'native', request, real)
    assert report['status'] == 'PASS', report['errors']
    folder = Path(report['directory']); loaded = []
    for i, analysis in enumerate(request['analyses']):
        data = np.loadtxt(folder/f'analysis{i}.txt', skiprows=1, ndmin=2)
        assert np.isfinite(data).all() and data.shape[1] == 1+len(analysis['vectors'])*(2 if analysis['kind'] == 'ac' else 1)
        if analysis['kind'] != 'op': assert np.all(np.diff(data[:, 0]) > 0)
        loaded.append(data)
    def column(index, name):
        analysis = request['analyses'][index]; j = analysis['vectors'].index(name); data = loaded[index]
        return data[:, 1+2*j]+1j*data[:, 2+2*j] if analysis['kind'] == 'ac' else data[:, 1+j]
    assert len(loaded[0]) == 1 and len(loaded[1]) == 401
    assert abs(loaded[1][0, 0]-1e3) < 1e-8 and abs(loaded[1][-1, 0]-1e8) < .01
    max_kcl = 0.
    for i, analysis in enumerate(request['analyses']):
        if analysis['kind'] not in ['op', 'tran']: continue
        for nodes in d['terminals'].values():
            voltages = [np.zeros(len(loaded[i])) if node == '0' else column(i, 'v('+node+')') for node in nodes]
            assert np.ptp(np.stack(voltages), axis=0).max() <= 1+1e-8
        for sensors in d['sensors'].values():
            max_kcl = max(max_kcl, float(max(abs(sum(column(i, 'i('+s+')') for s in sensors)))))
    assert max_kcl < 1e-10
    op = {v: float(column(0, v)[0]) for v in request['analyses'][0]['vectors']}
    readbacks = resistor_checks(op, d); gain = column(1, 'v(out)')/column(1, 'v(in)')
    measured = dict(out_v=op['v(out)'], gain=float(gain[0].real),
                    current_a=max(0, -op['i(Vdd)'])+max(0, -op['i(Vin)']))
    errors = {k: measured[k]-float(expected[k]) for k in measured}
    assert abs(errors['out_v']) < 1e-6 and abs(errors['gain']) < 1e-3 and abs(errors['current_a']) < 1e-9, errors
    result = dict(status='PASS', candidate=a.candidate, index=a.index, condition=a.condition, measured=measured,
        differences_from_published=errors, dc_ac_current_specification_pass=expected['dc_ac_current_pass'] == 'True',
        actual_resistor_voltage_current_checks=readbacks, max_device_kcl_residual_a=max_kcl,
        device_readback='Pinned APM verified actual model/W/L/m/nf/DELVTO/MULU0 before/applied/after every analysis.',
        apm_commit=commit, run_id=report['run_id'], reproduction_status='PASS')
    if signal:
        i = 2; t = loaded[i][:, 0]; assert t[0] == 0 and abs(t[-1]-.001) < 1e-12 and max(np.diff(t)) <= 1.00001e-7
        vo = column(i, 'v(out)'); vi = column(i, 'v(in)'); s = periodic(t, vo); si = periodic(t, vi); mask = t >= .0002
        assert abs(si['fundamental_peak_v']/.02-1) < 1e-4
        measured_signal = dict(thd_fraction=s['thd_fraction'], fundamental_gain=s['fundamental_peak_v']/si['fundamental_peak_v'],
                               output_min_v=float(min(vo[mask])), output_max_v=float(max(vo[mask])))
        differences = {k: value-float(expected_signal[k]) for k, value in measured_signal.items()}
        assert abs(differences['thd_fraction']) < 1e-5 and abs(differences['fundamental_gain']) < 1e-3
        assert max(abs(differences[k]) for k in ['output_min_v', 'output_max_v']) < 1e-5
        params = d['parameters']; cap = column(i, 'i(Vcap)'); kcl = cap+vo/params['rload']
        for name, nodes in d['terminals'].items():
            if nodes[0] == 'out': kcl += column(i, 'i('+d['sensors'][name][0]+')')
        if params.get('rd'): kcl += (vo-column(i, 'v(vdd)'))/params['rd']
        if params.get('rf'): kcl += (vo-column(i, 'v(gate)'))/params['rf']
        assert max(abs(kcl)) < 1e-10 and max(abs(vo-column(i, 'v(capnode)'))) < 1e-10
        axis = np.r_[.0009, t[(t > .0009) & (t < .000925)], .000925]
        charge = float(np.trapz(np.interp(axis, t, cap), axis))
        inferred_c = charge/float(np.interp(axis[-1], t, vo)-np.interp(axis[0], t, vo))
        assert abs(inferred_c/params['cload']-1) < .001
        passed = s['stationary'] and s['thd_fraction'] <= .02 and abs(measured_signal['fundamental_gain']/6-1) <= .1 and measured_signal['output_min_v'] >= .25 and measured_signal['output_max_v'] <= .75
        assert passed == (expected_signal['dynamic_tests_pass'] == 'True')
        result.update(signal=measured_signal, signal_differences_from_published=differences,
            dynamic_specification_pass=bool(passed), joint_dc_ac_signal_pass=expected_signal['joint_dc_ac_signal_pass'] == 'True',
            effective_capacitance_from_charge_f=inferred_c, max_output_transient_kcl_a=float(max(abs(kcl))),
            scope='Saved preselected physical case, OP/AC and 20-mV sine. Charge/voltage capacitor check; no noise or population-dynamic replay.')
    else:
        result['scope'] = 'Saved physical OP/AC case. No noise, sine or full-population native replay.'
    (out/'result.json').write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
