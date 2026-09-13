"""MIT: four selected nominal DC/PWL replays through the pinned public APM API.

Exact saved units and circuit bodies are relocated, without drawing or trimming.
This is an example, not a campaign runner. SPICE is separate from site rendering.
"""
import argparse
import copy
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

EXAMPLES = ['F75L-input_only', 'F75L-supply_first-fast-8p', 'R75-input_first-fast-8p', 'F75L-supply-cycle-8p']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def entry(t, y, target, origin, horizon):
    """Same final-entry sample rule as the frozen measurement; relative bracket."""
    import numpy as np
    mask = (t >= origin) & (t <= horizon); x = t[mask]; e = np.abs(y[mask]-target)-.001
    assert len(x) >= 3 and t[-1] >= horizon*(1-1e-12)
    outside = np.flatnonzero(e > 0)
    if len(outside) and outside[-1] == len(x)-1: return dict(status='RIGHT_CENSORED', time_s=None)
    if not len(outside):
        return dict(status='IN_BAND_AT_FIRST_OBSERVATION', time_s=float(x[0]-origin), relative_bracket_s=[0., float(x[0]-origin)])
    j = int(outside[-1])
    return dict(status='OBSERVED_FINAL_ENTRY', time_s=float(x[j+1]-origin), relative_bracket_s=[float(x[j]-origin), float(x[j+1]-origin)])


def resistor_readback(observations, cfg):
    """Infer each resistance from measured node voltage and independent current.

Use the largest observed current for conditioning. Explicitly disclose a
resistor whose current is too small; never substitute its declared value as
an actual readback. DC alone cannot measure an inactive capacitor.
"""
    import numpy as np
    p = cfg['parameters']; terminals = cfg['terminals']; sensors = cfg['sensors']; candidates = {}
    for v in observations:
        def current(name, i): return v['i('+sensors[name][i]+')']
        incoming = -v['i(Vin)']; gates = sum((current(n, 1) for n, nodes in terminals.items() if nodes[1] == 'gate'), np.zeros(len(incoming)))
        feedback = gates-incoming if 'rf' in p else np.zeros(len(incoming))
        source_sum = sum((current(n, 2) for n, nodes in terminals.items() if nodes[2] == 'source'), np.zeros(len(incoming)))
        output_sum = sum((current(n, i) for n, nodes in terminals.items() for i, node in enumerate(nodes) if node == 'out'), np.zeros(len(incoming)))
        flows = dict(Rsrc=(v['v(in)']-v['v(drive)' if 'rf' in p else 'v(gate)'], incoming, p['rsrc']),
                     Rs=(v['v(source)'], -source_sum, p['rs']),
                     Rload=(v['v(out)'], (-v['i(Vdd)'] if 'rd' in p else 0.)-output_sum-feedback-v['i(Vcap)'], p['rload']))
        if 'rf' in p:
            flows.update(Rin=(v['v(drive)']-v['v(gate)'], incoming, p['rin']), Rf=(v['v(out)']-v['v(gate)'], feedback, p['rf']))
        if 'rd' in p: flows['Rd'] = (v['v(vdd)']-v['v(out)'], -v['i(Vdd)'], p['rd'])
        if 'rbias' in p:
            into_bias = sum((current(n, i) for n, nodes in terminals.items() for i, node in enumerate(nodes) if node == 'pbias'), np.zeros(len(incoming)))
            flows['Rbias'] = (v['v(pbias)'], -v['i(Vbiasdiag)']-into_bias, p['rbias'])
        for name, (drop, flow, expected) in flows.items():
            j = int(np.argmax(abs(flow))); selected = (float(drop[j]), float(flow[j]), expected)
            if name not in candidates or abs(selected[1]) > abs(candidates[name][1]): candidates[name] = selected
    checks = {}
    for name, (drop, flow, expected) in candidates.items():
        if abs(flow) <= 1e-12:
            checks[name] = dict(status='V_I_READBACK_UNRESOLVED', reason='Observed current at most 1 pA; V/I is ill-conditioned.', maximum_observed_current_a=abs(flow), declared_ohm=expected)
        else:
            actual = drop/flow; assert abs(actual/expected-1) < 1e-5, (name, actual, expected)
            checks[name] = dict(status='PASS', measured_v_over_i_ohm=actual, declared_ohm=expected, relative_difference=actual/expected-1)
    return checks


def main():
    p = argparse.ArgumentParser(); p.add_argument('--apm', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--binary', type=Path, default=Path('/usr/local/bin/ngspice')); p.add_argument('--case', choices=EXAMPLES, default=EXAMPLES[1])
    a = p.parse_args(); here = Path(__file__).resolve().parent; apm = a.apm.resolve(); out = a.output
    if not out.is_absolute() or out.exists() or out.resolve().is_relative_to(apm): raise ValueError('Use a new absolute output outside the model checkout')
    cfg = json.loads((here/'data/configuration.json').read_text()); d = cfg['cases'][a.case]; case = d['case']
    for name, digest in json.loads((here/'data/manifest.json').read_text())['files'].items(): assert sha(here/'data'/name) == digest
    assert subprocess.check_output(['git', '-C', str(apm), 'rev-parse', 'HEAD'], text=True).strip() == cfg['apm_commit']
    for name, digest in cfg['model_hashes'].items(): assert sha(apm/name) == digest
    body = here/'circuits'/(a.case+'.cir'); assert sha(body) == d['circuit_sha256'] and not case['refined']
    expected = next(r for r in json.loads((here/'data/observations.json').read_text())['results'] if r['case']['id'] == a.case)
    os.environ['APM_REPO_ROOT'] = str(apm); sys.path.insert(0, str(apm/'src'))
    import numpy as np
    from apm import research, research_spice as spice
    request = dict(schema=research.SCHEMAS['request'], circuit=str(body), devices=d['devices'], other_variation_leaves=[], analyses=copy.deepcopy(d['analyses']))
    for analysis in request['analyses']:
        analysis['vectors'] = [v for v in analysis['vectors'] if not v.startswith('@')]
        if analysis['kind'] == 'tran':
            # Public API emits tran step stop. Its implicit tmax=min(step,stop/50)
            # equals the explicit frozen tmax for every advertised example.
            assert min(analysis['step'], analysis['stop']/50) == analysis['tmax']
    binding = {}; spice.flatten(body, binding)
    for name in spice.MODELS: spice.flatten(apm/name, binding)
    units = copy.deepcopy(d['physical_devices']); assert all(u['raw'] == [0., 0.] for u in units)
    real = research.seal(dict(schema=research.SCHEMAS['realization'], origin='research', profile_tier='SOURCE_TRANSFER_HYPOTHESIS', status='RESOLVED',
        sample_context_id=research.sample_context(request), request_id=research.canonical_hash(request), input_binding=binding, devices=units,
        public_replay=dict(original_realization_id=d['original_realization_id'], saved_nominal_id=d['saved_nominal_id'], rule='Saved nominal raw/geometry/UID; only file bindings relocated.')))
    out.mkdir(parents=True); (out/'relocated-realization.json').write_text(json.dumps(real, indent=2)+'\n')
    report = spice.execute(apm, a.binary, out/'native', request, real, temperature_c=cfg['temperature_c']); assert report['status'] == 'PASS', report['errors']
    folder = Path(report['directory']); observations = []; axes = []; kcl = 0.
    for i, analysis in enumerate(request['analyses']):
        arr = np.loadtxt(folder/f'analysis{i}.txt', skiprows=1, ndmin=2)
        assert np.isfinite(arr).all() and arr.shape[1] == len(analysis['vectors'])+1
        v = {name: arr[:, j+1] for j, name in enumerate(analysis['vectors'])}; observations.append(v); axes.append(arr[:, 0])
        for nodes in d['terminals'].values():
            voltages = [np.zeros(len(arr)) if n == '0' else v['v('+n+')'] for n in nodes]
            assert np.ptp(np.stack(voltages), axis=0).max() <= 1+1e-8
        for sensors in d['sensors'].values(): kcl = max(kcl, float(max(abs(sum(v['i('+s+')'] for s in sensors)))))
        assert max(abs(v['v(out)']-v['v(capnode)'])) < 1e-10
    assert kcl < 1e-10 and len(axes[0]) == 1
    op = observations[0]; saved_op = expected['measurements']['static']['op']
    op_differences = {name: float(v[0])-saved_op[name] for name, v in op.items()}
    assert max(abs(v) for k, v in op_differences.items() if k.startswith('v(')) < 1e-6
    assert max(abs(v) for k, v in op_differences.items() if k.startswith('i(')) < 1e-9
    result = dict(status='PASS', case=a.case, apm_commit=cfg['apm_commit'], run_id=report['run_id'],
        op_differences=op_differences, resistor_readbacks=resistor_readback(observations, d), maximum_device_kcl_a=kcl,
        actual_mos_readback='Pinned APM checks actual model, W/L/m/nf, DELVTO and MULU0 before/applied/after every analysis.',
        scope='Known saved nominal replay on the tested host. No redraw, Local off-state qualification or all-history claim.')
    if case['kind'] != 'static':
        t = axes[1]; v = observations[1]; vo = v['v(out)']; tr = expected['measurements']['transient']; par = d['parameters']
        assert t[0] == 0 and abs(t[-1]-case['stop_s']) < 1e-12 and np.all(np.diff(t) > 0) and max(np.diff(t)) <= case['tmax_s']*(1+1e-6)
        for key, vector in [('input_points', 'v(in)'), ('supply_points', 'v(vdd)')]:
            points = np.array(case[key]); assert max(abs(v[vector]-np.interp(t, points[:, 0], points[:, 1]))) < 1e-9
        difference = dict(out_min_v=float(min(vo))-tr['out_min_v'], out_max_v=float(max(vo))-tr['out_max_v'])
        assert max(abs(x) for x in difference.values()) < 1e-6
        gates = sum((v['i('+d['sensors'][name][1]+')'] for name, nodes in d['terminals'].items() if nodes[1] == 'gate'), np.zeros(len(t)))
        input_model = ((v['v(in)']-vo)+par['rf']*gates)/(par['rsrc']+par['rin']+par['rf']) if 'rf' in par else gates
        input_residual = float(max(abs(-v['i(Vin)']-input_model))); assert input_residual < 1e-10
        output_kcl = v['i(Vcap)']+vo/par['rload']
        for name, nodes in d['terminals'].items():
            for j, node in enumerate(nodes):
                if node == 'out': output_kcl += v['i('+d['sensors'][name][j]+')']
        if 'rf' in par: output_kcl += (vo-v['v(gate)'])/par['rf']
        if 'rd' in par: output_kcl += (vo-v['v(vdd)'])/par['rd']
        assert max(abs(output_kcl)) < 1e-10
        port_checks = {}
        for source, node in [('Vin', 'in'), ('Vdd', 'vdd')]:
            flow = -v['i('+source+')']; power = v['v('+node+')']*flow
            measured = dict(min_delivered_current_a=float(min(flow)), max_delivered_current_a=float(max(flow)),
                signed_energy_j=float(np.trapz(power, t)), positive_delivered_energy_j=float(np.trapz(np.maximum(power, 0), t)), absorbed_energy_j=float(np.trapz(np.maximum(-power, 0), t)))
            errors = {k: val-tr['ports'][source][k] for k, val in measured.items()}
            for k, val in errors.items(): assert abs(val) < (1e-9 if k.endswith('_a') else max(1e-17, abs(tr['ports'][source][k])*1e-4)), (source, k, val)
            port_checks[source] = dict(measured=measured, differences=errors)
        entries = []
        for plateau in tr['plateaus']:
            measured = entry(t, vo, plateau['target_v'], plateau['origin_s'], plateau['horizon_s']); saved = plateau['settling']
            assert measured['status'] == saved['status']
            if measured['time_s'] is not None: assert abs(measured['time_s']-saved['time_s']) <= max(2*case['tmax_s'], .01*saved['time_s'])
            measured['endpoint_error_v'] = float(np.interp(plateau['horizon_s'], t, vo))-plateau['target_v']
            assert abs(measured['endpoint_error_v']-plateau['endpoint_error_v']) < 1e-6
            entries.append(measured)
        charges = []
        for check in tr['charge_checks']:
            start, stop = check['interval_s']; axis = np.r_[start, t[(t > start) & (t < stop)], stop]
            actual = float(np.trapz(np.interp(axis, t, v['i(Vcap)']), axis))/(float(np.interp(stop, t, vo))-float(np.interp(start, t, vo)))
            assert abs(actual/par['cload']-1) < .001
            charges.append(dict(measured_integral_i_over_delta_v_f=actual, declared_f=par['cload'], relative_difference=actual/par['cload']-1))
        result['transient'] = dict(output_extrema_differences=difference, ports=port_checks, plateaus=entries, capacitor_readback=charges,
            maximum_input_identity_residual_a=input_residual, maximum_output_kcl_a=float(max(abs(output_kcl))))
    else:
        result['capacitor_readback'] = 'DC only: capacitance is inactive and not independently measured by this public replay; private acquisition included direct parameter readback.'
    (out/'replay-result.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], case=a.case, output=str(out/'replay-result.json'))))


if __name__ == '__main__': main()
