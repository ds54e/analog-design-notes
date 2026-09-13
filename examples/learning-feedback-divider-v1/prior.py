"""MIT: old-OP affine forecast for one finite initial-D1 sensing divider."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parent
NODES = ['tail', 'mirror', 'out', 'bias', 'inn', 'inp']


def read_source(path):
    source = json.loads((ROOT/'data/source.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == source['package_sha256']
    with zipfile.ZipFile(path) as z:
        for name, h in source['selected_zip_members'].items():
            assert hashlib.sha256(z.read('example/'+name)).hexdigest() == h
        case = json.loads(z.read('example/data/interface-inputs.json'))['cases']['D1-5p']
        with np.load(io.BytesIO(z.read('example/data/'+case['array_file'])), allow_pickle=False) as a:
            op_array, ac_array = a['analysis1'], a['analysis3']
    op_recipe, ac_recipe = case['analyses'][1]['recipe'], case['analyses'][3]['recipe']
    assert op_array.shape == (1, 1+len(op_recipe['vectors']))
    assert ac_array.shape == (361, 1+2*len(ac_recipe['vectors']))
    assert np.isfinite(op_array).all() and np.isfinite(ac_array).all()
    op = {n: float(op_array[0, 1+i]) for i, n in enumerate(op_recipe['vectors'])}
    values = {n: ac_array[:, 1+2*i]+1j*ac_array[:, 2+2*i] for i, n in enumerate(ac_recipe['vectors'])}
    assert op['@rload[resistance]'] == 1e6 and op['@cload[capacitance]'] == 5e-12
    assert abs(op['v(source)']-.3) <= 1e-9 and abs(op['v(vdd)']-1.) <= 1e-9
    return case, op, ac_array[:, 0], values['v(out)']/values['v(source)']


def derivatives(op, index):
    return [op[f'@m.xm{index}.mapm045_vtg_core[{v}]'] for v in ('gm', 'gmbs', 'gds')]


def stamp(op, connections, fixed):
    """Full six-node current residual and drain/source derivative stamp."""
    m = np.zeros((len(NODES), len(NODES)))
    residual = np.zeros(len(NODES))
    drive = np.zeros(len(NODES))
    for name, c in connections.items():
        gm, gb, gd = derivatives(op, int(name.lstrip('xm')))
        for terminal, node in c['nodes'].items():
            if node in NODES:
                residual[NODES.index(node)] += op['i('+c['sensors'][terminal]+')']
        for terminal, sign in [('d', 1.), ('s', -1.)]:
            node = c['nodes'][terminal]
            if node not in NODES:
                continue
            row = NODES.index(node)
            for port, value in [('d', gd), ('g', gm), ('s', -gd-gm-gb), ('b', gb)]:
                other = c['nodes'][port]
                if other in NODES:
                    m[row, NODES.index(other)] += sign*value
                else:
                    assert other in ('0', 'vdd')
    resistors = [('out', '0', fixed['rload_ohm']), ('bias', '0', fixed['rbias_ohm']),
                 ('inp', 'source', fixed['rsource_ohm']),
                 ('out', 'inn', fixed['rupper_ohm']), ('inn', '0', fixed['rlower_ohm'])]
    voltage = lambda n: 0. if n == '0' else op['v('+n+')']
    for n1, n2, resistance in resistors:
        for a, b in [(n1, n2), (n2, n1)]:
            if a not in NODES:
                continue
            row = NODES.index(a)
            residual[row] += (voltage(a)-voltage(b))/resistance
            m[row, row] += 1/resistance
            if b in NODES:
                m[row, NODES.index(b)] -= 1/resistance
            elif b == 'source':
                drive[row] += 1/resistance
            else:
                assert b == '0'
    return m, residual, drive


def explicit_reduction(op, fixed):
    """Three independent tail/mirror/output equations; eliminate the divider."""
    gm1, gb1, gd1 = derivatives(op, 1)
    gm2, gb2, gd2 = derivatives(op, 2)
    gm3, _, gd3 = derivatives(op, 3)
    gm4, _, gd4 = derivatives(op, 4)
    _, _, gd5 = derivatives(op, 5)
    s1, s2 = gm1+gb1+gd1, gm2+gb2+gd2
    total = fixed['rupper_ohm']+fixed['rlower_ohm']
    beta = fixed['rlower_ohm']/total
    m = np.array([[s1+s2+gd5, -gd1, -gd2-beta*gm2],
                  [-s1, gd1+gm3+gd3, 0.],
                  [-s2, gm4, gd2+gd4+1/fixed['rload_ohm']+1/total+beta*gm2]])
    b = np.array([gm1, -gm1, 0.])
    h = float(np.linalg.solve(m, b)[2])
    z = float(np.linalg.solve(m, [0., 0., 1.])[2])
    # A+/A- belong to the open-input OTA with the divider's DC loading included.
    opened = m.copy()
    opened[0, 2] += beta*gm2
    opened[2, 2] -= beta*gm2
    ap = float(np.linalg.solve(opened, b)[2])
    am = float(np.linalg.solve(opened, [gm2, 0., -gm2])[2])
    assert np.isclose(h, ap/(1-beta*am), rtol=1e-12)
    return dict(gain=h, output_transimpedance_ohm=z, beta=beta,
                a_plus=ap, a_minus=am, feedback_denominator=1-beta*am)


def crossing(f, h):
    db = 20*np.log10(abs(h)/abs(h[0]))
    i = int(np.flatnonzero((db[:-1] > -3) & (db[1:] <= -3))[0])
    return float(np.exp(np.log(f[i])+(-3-db[i])/(db[i+1]-db[i])*np.log(f[i+1]/f[i])))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    case, op, f, h = read_source(a.source_zip)
    fixed = json.loads((ROOT/'plan.json').read_text())['fixed']
    m, residual, drive = stamp(op, case['connections'], fixed)
    delta = np.linalg.solve(m, -residual)
    forecast_nodes = {n: float(op['v('+n+')']+delta[i]) for i, n in enumerate(NODES)}
    linear_gain = float(np.linalg.solve(m, drive)[NODES.index('out')])
    output_drive = np.array([1. if n == 'out' else 0. for n in NODES])
    z = float(np.linalg.solve(m, output_drive)[NODES.index('out')])
    reduced = explicit_reduction(op, fixed)
    assert np.isclose(linear_gain, reduced['gain'], rtol=1e-12)
    assert np.isclose(z, reduced['output_transimpedance_ohm'], rtol=1e-12)
    assert max(abs(m@delta+residual)) < 1e-15
    # DC sensing relation independently includes the known gate-terminal current.
    gate = op['i(VtM2g)']
    parallel = 1/(1/fixed['rupper_ohm']+1/fixed['rlower_ohm'])
    sense = reduced['beta']*forecast_nodes['out']-parallel*gate
    assert abs(sense-forecast_nodes['inn']) < 1e-12
    pole = 1/(2*np.pi*fixed['cload_f']*z)
    crossing_hz = pole*np.sqrt(10**.3*(1+(1000/pole)**2)-1)
    beta = reduced['beta']
    result = dict(schema='learning.initial-d1-feedback-divider.prior.v1',
        role='Old-OP affine/derivative forecast saved before the one new condition; no new target read.',
        source_sha256=hashlib.sha256(a.source_zip.read_bytes()).hexdigest(),
        old_case='D1-5p', old_realization_id=case['original_realization_id'],
        known=dict(observed_op=op, gain_one_khz=[float(h[0].real), float(h[0].imag)], bandwidth_hz=crossing(f, h)),
        node_order=NODES, conductance_matrix_s=m.tolist(), changed_connection_current_residual_a=residual.tolist(),
        signal_drive_a_per_v=drive.tolist(), predicted_node_changes_v=delta.tolist(), independent_reduction=reduced,
        forecast=dict(node_v=forecast_nodes, output_v=forecast_nodes['out'],
            ideal_output_v=.3/beta, output_error_to_ideal_v=forecast_nodes['out']-.3/beta,
            gain=linear_gain, ideal_gain=1/beta, output_transimpedance_ohm=z,
            single_output_pole_hz=float(pole), bandwidth_minus3db_hz=float(crossing_hz),
            input_common_mode_v=(forecast_nodes['inp']+forecast_nodes['inn'])/2,
            input_difference_v=forecast_nodes['inp']-forecast_nodes['inn'],
            load_current_a=forecast_nodes['out']/fixed['rload_ohm'],
            upper_current_a=(forecast_nodes['out']-forecast_nodes['inn'])/fixed['rupper_ohm'],
            lower_current_a=forecast_nodes['inn']/fixed['rlower_ohm'],
            sensed_voltage_from_independent_KCL_v=sense, known_gate_current_a=gate),
        assumptions=['Old measured terminal currents form the affine intercept; old gm/gmb/gds remain fixed.',
            'DC includes known gate/body terminal currents; their derivatives and all intrinsic capacitances are omitted.',
            'The finite source and bias resistors remain in the six-node system. The gate/source voltage derivative is1 only under the stated gate-current approximation.',
            'Only the external output capacitor contributes to the one-pole estimate. This is neither a return ratio nor a stability/phase-margin measurement.',
            'The new output differs substantially from the old0.298456V; retain prediction errors from changed headroom/derivatives.'])
    a.output.mkdir(parents=True)
    (a.output/'prior.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result['forecast'], indent=2))


if __name__ == '__main__':
    main()
