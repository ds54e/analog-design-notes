"""MIT: one known initial-D1 high-command replay using only public APM inputs."""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from analyze import ROOT, bandwidth, columns, inputs, read, sha


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source-zip', type=Path, required=True)
    p.add_argument('--apm', type=Path, required=True)
    p.add_argument('--binary', type=Path, default=Path('/usr/local/bin/ngspice'))
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert a.output.is_absolute() and not a.output.exists()
    apm = a.apm.resolve()
    assert not a.output.resolve().is_relative_to(apm)
    plan, _ = inputs(a.source_zip)
    replay = read(ROOT/'replay-plan.json')
    assert replay['script_sha256'] == sha(__file__)
    assert replay['scientific_plan_sha256'] == sha(ROOT/'plan.json')
    version = subprocess.check_output([str(a.binary), '--version'], text=True).strip()
    assert 'ngspice-47' in version
    assert subprocess.check_output(['git', '-C', str(apm), 'rev-parse', 'HEAD'], text=True).strip() == plan['apm_commit']
    for name, expected in plan['model_hashes'].items():
        assert sha(apm/name) == expected
    name = 'D1-high-command'
    case = plan['cases'][name]
    body = ROOT/case['circuit']
    lines = {line.split()[0]: line.split() for line in body.read_text().splitlines() if line and not line.startswith('*')}
    recipes = copy.deepcopy(case['analyses'])
    for recipe in recipes:
        recipe['vectors'] = [v for v in recipe['vectors'] if not v.startswith('@')]
        for element, value in recipe.pop('elements', {}).items():
            assert float(lines[element][-1]) == value
        for source, value in recipe.pop('ac_sources', {}).items():
            assert float(lines[source][lines[source].index('AC')+1]) == value
        for source, value in recipe['set_sources'].items():
            assert float(lines[source][lines[source].index('DC')+1]) == value
    os.environ['APM_REPO_ROOT'] = str(apm)
    sys.path.insert(0, str(apm/'src'))
    from apm import research, research_spice as spice
    request = dict(schema=research.SCHEMAS['request'], circuit=str(body),
        devices=[{k: v for k, v in u.items() if k != 'raw'} for u in plan['physical_units']],
        other_variation_leaves=[], analyses=recipes)
    binding = {}
    spice.flatten(body, binding)
    for name in spice.MODELS:
        spice.flatten(apm/name, binding)
    index = read(ROOT/'data/index.json')
    original = index['cases']['D1-high-command']
    assert sha(ROOT/'data'/original['file']) == original['sha256']
    real = research.seal(dict(schema=research.SCHEMAS['realization'], origin='research', profile_tier='ARTIFICIAL',
        status='RESOLVED', request_id=research.canonical_hash(request), sample_context_id=research.sample_context(request),
        input_binding=binding, devices=plan['physical_units'],
        public_replay=dict(original_realization_id=original['realization_id'],
            rule='Relocated paths and public terminal-vector recipes; same saved six zero-raw units and high-command fixture. No redraw, rebias or fresh confirmation.')))
    a.output.mkdir(parents=True)
    for filename, value in [('request.json', request), ('realization.json', real), ('replay-plan.json', replay)]:
        (a.output/filename).write_text(json.dumps(value, indent=2)+'\n')
    report = spice.execute(apm, a.binary, a.output/'native', request, real, temperature_c=plan['temperature_c'])
    assert report['status'] == 'PASS', report
    folder = Path(report['directory'])
    saved = dict(np.load(ROOT/'data'/original['file'], allow_pickle=False))
    differences, arrays, values = {}, {}, []
    kcl, span = 0., 0.
    for i, recipe in enumerate(recipes):
        z = np.loadtxt(folder/f'analysis{i}.txt', skiprows=1, ndmin=2)
        axis, v = columns(recipe, z)
        oldaxis, old = columns(case['analyses'][i], saved['analysis'+str(i)])
        assert len(axis) == len(oldaxis) and np.allclose(axis, oldaxis, rtol=1e-12, atol=1e-12)
        if i:
            assert np.all(np.diff(axis) > 0) and len(axis) == 361
        error = {key: float(max(abs(value-old[key]))) for key, value in v.items()}
        for key, delta in error.items():
            limit = replay['tolerances']['voltage_absolute_v'] if key.startswith('v(') else replay['tolerances']['current_absolute_a']
            assert delta <= limit, (key, delta, limit)
        differences[str(i)] = dict(kind=recipe['kind'], rows=len(axis),
            maximum_voltage_difference_v=max(delta for key, delta in error.items() if key.startswith('v(')),
            maximum_current_difference_a=max(delta for key, delta in error.items() if key.startswith('i(')))
        for c in case['connections'].values():
            kcl = max(kcl, float(max(abs(sum(v['i('+s+')'] for s in c['sensors'].values())))))
            if i == 0:
                volts = [0 if n == '0' else float(v['v('+n+')'][0]) for n in c['nodes'].values()]
                span = max(span, max(volts)-min(volts))
        arrays['analysis'+str(i)] = z
        values.append(v)
    tol = plan['measurements']
    assert kcl <= tol['kcl_max_a'] and span <= tol['terminal_span_max_v']
    op, ac = values
    o = {k: float(v[0]) for k, v in op.items()}
    fixed = plan['fixed']
    avail = -o['i(VtM2d)']-o['i(VtM4d)']
    output_kcl = avail-o['v(out)']/fixed['rload_ohm']+o['i(Vfeedback)']
    input_kcl = -o['i(Vin)']-o['i(VtM1g)']
    assert max(abs(output_kcl), abs(input_kcl)) <= tol['kcl_max_a']
    flows = dict(Rsource=(o['v(source)']-o['v(inp)'], -o['i(Vin)'], fixed['rsource_ohm']),
        Rload=(o['v(out)'], avail+o['i(Vfeedback)'], fixed['rload_ohm']),
        Rbias=(o['v(bias)'], -o['i(VtM6d)']-o['i(VtM6g)']-o['i(VtM5g)'], fixed['rbias_ohm']))
    resistors = {}
    for name, (drop, current, expected) in flows.items():
        assert abs(current) >= tol['vi_minimum_current_a'], 'Retain ill-conditioned native readback; do not substitute the expected value'
        actual = drop/current
        delta = actual/expected-1
        assert abs(delta) <= tol['passive_vi_relative']
        resistors[name] = dict(observed_v_over_i_ohm=actual, expected_ohm=expected, relative_difference=delta)
    f = arrays['analysis1'][:, 0]
    h = ac['v(out)']/ac['v(source)']
    icap = -ac['i(VtM2d)']-ac['i(VtM4d)']-ac['v(out)']/fixed['rload_ohm']+ac['i(Vfeedback)']
    measured_c = icap/(2j*np.pi*f*ac['v(out)'])
    mask = (f >= 1e5) & (f <= 1e7)
    cap_error = float(max(abs(measured_c[mask]/fixed['cload_f']-1)))
    assert cap_error <= tol['capacitor_ac_relative']
    np.savez_compressed(a.output/'observations.npz', **arrays)
    result = dict(schema='learning.initial-d1-common-mode.public-replay.v1', status='PASS',
        role='One previously exposed high-command initial-D1 reproduction; not new scientific confirmation.',
        new_run_id=report['run_id'], realization_id=real['content_id'], original_realization_id=original['realization_id'],
        apm_commit=plan['apm_commit'], ngspice_version=version, binary_sha256=sha(a.binary),
        scientific_plan_sha256=sha(ROOT/'plan.json'), replay_plan_sha256=sha(ROOT/'replay-plan.json'),
        differences=differences, output_v=o['v(out)'], error_v=o['v(out)']-o['v(source)'],
        gain_one_khz=[float(h[0].real), float(h[0].imag)], bandwidth=bandwidth(f, h),
        maximum_MOS_KCL_a=kcl, maximum_terminal_span_v=span,
        output_node_KCL_a=output_kcl, input_node_KCL_a=input_kcl,
        resistor_vi_readback=resistors, capacitor_ac_relative_error=cap_error,
        MOS_readback='Pinned public APM checks actual model/W/L/m/nf/DELVTO/MULU0 before/applied/after.',
        passive_readback='Measured DC V/I for all three resistors and complex AC current/voltage for Cload; no copied expected values as observations.',
        limits='One stationary high-command OP/AC example on the same core. No transient, common-mode range qualification, noise population, regulator or independent scientific replication.')
    (a.output/'result.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
