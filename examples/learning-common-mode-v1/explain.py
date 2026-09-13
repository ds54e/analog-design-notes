"""MIT: retrospective low-frequency node equations for the observed initial D1."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analyze import ROOT, inputs, read, sha


def derivatives(op, index):
    return [op[f'@m.xm{index}.mapm045_vtg_core[{v}]'] for v in ('gm', 'gmbs', 'gds')]


def explicit(op, rload):
    gm1, gb1, gd1 = derivatives(op, 1)
    gm2, gb2, gd2 = derivatives(op, 2)
    gm3, _, gd3 = derivatives(op, 3)
    gm4, _, gd4 = derivatives(op, 4)
    _, _, gd5 = derivatives(op, 5)
    s1, s2 = gm1+gb1+gd1, gm2+gb2+gd2
    # Unknowns: tail, mirror, output. The two gate voltages are independent.
    m = np.array([[s1+s2+gd5, -gd1, -gd2],
                  [-s1, gd1+gm3+gd3, 0.],
                  [-s2, gm4, gd2+gd4+1/rload]])
    rhs_plus = np.array([gm1, -gm1, 0.])
    rhs_minus = np.array([gm2, 0., -gm2])
    ap = float(np.linalg.solve(m, rhs_plus)[2])
    am = float(np.linalg.solve(m, rhs_minus)[2])
    # Feedback substitutes v_minus=v_out into the same linear equations.
    closed = m.copy()
    closed[:, 2] -= rhs_minus
    solved = np.linalg.solve(closed, rhs_plus)
    return dict(plus=ap, minus=am, gain=ap/(1-am), direct_gain=float(solved[2]),
                tail_per_input=float(solved[0]), mirror_per_input=float(solved[1]),
                tail_gds_s=gd5, input_source_coefficient_s=s1+s2,
                tail_conductance_fraction=gd5/(s1+s2+gd5))


def stamp(op, connections, fixed):
    # Independent construction: stamp every MOS drain/source terminal current.
    names = ['tail', 'mirror', 'out', 'bias']
    m, b = np.zeros((4, 4)), np.zeros(4)

    def node(n):
        return 'out' if n == 'inn' else n

    for name, connection in connections.items():
        index = int(name.lstrip('xm'))
        gm, gb, gd = derivatives(op, index)
        n = {k: node(v) for k, v in connection['nodes'].items()}
        coeff = [('d', gd), ('g', gm), ('s', -gd-gm-gb), ('b', gb)]
        for terminal, sign in [('d', 1.), ('s', -1.)]:
            if n[terminal] not in names:
                continue
            row = names.index(n[terminal])
            for voltage, value in coeff:
                if n[voltage] in names:
                    m[row, names.index(n[voltage])] += sign*value
                elif n[voltage] == 'inp':
                    b[row] -= sign*value
                else:
                    assert n[voltage] in ('0', 'vdd')
    m[names.index('bias'), names.index('bias')] += 1/fixed['rbias_ohm']
    m[names.index('out'), names.index('out')] += 1/fixed['rload_ohm']
    x = np.linalg.solve(m, b)
    return {n: float(v) for n, v in zip(names, x)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source-zip', type=Path, required=True)
    p.add_argument('--measurements', type=Path, default=ROOT/'data/summary.json')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    plan, prior = inputs(a.source_zip)
    measured = read(a.measurements)
    assert measured['plan_sha256'] == sha(ROOT/'plan.json')
    old = read(ROOT/'data/prior.json')['cases']['D1-5p']
    new = measured['cases']['D1-high-command']
    connections = plan['cases']['D1-high-command']['connections']
    rows = []
    for name, op, native in [('known-central-0p3', old['observed_op'], old['gain_one_khz']),
                              ('new-high-command-0p7', new['op'], new['gain_one_khz'])]:
        e = explicit(op, plan['fixed']['rload_ohm'])
        s = stamp(op, connections, plan['fixed'])
        assert np.isclose(e['gain'], e['direct_gain'], rtol=1e-12, atol=1e-12)
        assert np.isclose(e['gain'], s['out'], rtol=1e-12, atol=1e-12)
        assert np.isclose(e['tail_per_input'], s['tail'], rtol=1e-12, atol=1e-12)
        assert abs(s['bias']) < 1e-12
        rows.append(dict(case=name, **e, stamped_gain=s['out'], native_one_khz_real=native[0],
            native_one_khz_imag=native[1], relative_error_to_native_real=e['gain']/native[0]-1,
            body_derivatives_retained=True))
    result = dict(schema='learning.initial-d1-common-mode.network.v1', status='ALGEBRA_CHECKS_PASS',
        role='Retrospective explanation using already observed OP derivatives; no new native target or prior prediction.',
        cases=rows, measurements_sha256=sha(a.measurements),
        assumptions=['Low frequency; all capacitances and gate/body leakage derivatives neglected.',
            'VDD is AC ground; the finite reference makes bias AC ground when its gate leakage is neglected.',
            'Input resistance has no low-frequency drop when gate leakage derivatives are omitted. Actual DC source drop remains in the measured table.',
            'Ap/Am are model-derived partial coefficients at the retained D1 operating point, not measured gains of a newly rebiased open circuit.',
            'Use h=Ap/(1-Am). Substituting one generic differential gain into A/(1+A) also assumes equal and opposite input coefficients.'],
        limits='The comparison checks local signal gain near1kHz; it does not predict broadband poles, finite overload/recovery, common-mode rejection statistics or a full operating range.')
    a.output.mkdir(parents=True)
    (a.output/'network.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    with (a.output/'network.csv').open('w') as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
