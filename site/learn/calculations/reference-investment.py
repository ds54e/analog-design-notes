"""MIT: finite saved-data arithmetic for A37-to-A75 reference investment.

No simulator, statistical refit, parameter search or new physical sample.
The source ZIP is the immutable first study; all output destinations are new.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from pathlib import Path

SOURCE_URL = 'https://ds54e.github.io/analog-design-notes/downloads/gain-study-v1.0.zip'
SOURCE_SHA256 = '4d2dda32baa15f173c6a897e7429ba0e421b0fa1738ff722fe8b889d1cf3e55f'
SELECTED = ['data/configuration.json', 'data/nominal.csv', 'data/source-budget.csv',
            'data/sensitivities.csv', 'data/statistics.json', 'data/contrasts.csv',
            'data/manifest.json', 'confirmation-plan.json',
            'development-contract-v1.json', 'circuits/A37.cir', 'circuits/A75.cir']


def digest(data):
    return hashlib.sha256(data).hexdigest()


def close(x, y, relative=1e-12, absolute=1e-15):
    assert math.isclose(x, y, rel_tol=relative, abs_tol=absolute), (x, y)


def circuit_inventory(text):
    units = []
    resistors = {}
    for line in text.splitlines():
        tokens = line.split()
        if not tokens or tokens[0].startswith('*'):
            continue
        if tokens[0].lower().startswith('x'):
            assert len(tokens) == 8, line
            name, drain, gate, source, body, model, w, l = tokens
            assert w.startswith('w=') and l.startswith('l='), line
            units.append(dict(name=name, nodes=[drain, gate, source, body],
                              model=model, w_m=float(w[2:]), l_m=float(l[2:])))
        elif tokens[0].lower().startswith('r'):
            assert len(tokens) == 4, line
            resistors[tokens[0]] = float(tokens[3])
    return units, resistors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-zip', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    source_bytes = args.source_zip.read_bytes()
    assert digest(source_bytes) == SOURCE_SHA256, 'Wrong historical package'
    with zipfile.ZipFile(io.BytesIO(source_bytes)) as archive:
        selected = {name: archive.read('example/' + name) for name in SELECTED}
    text = lambda name: selected[name].decode('utf-8')
    js = lambda name: json.loads(text(name))
    rows = lambda name: list(csv.DictReader(io.StringIO(text(name))))
    manifest = js('data/manifest.json')
    for name, data in selected.items():
        if name.startswith('data/') and name != 'data/manifest.json':
            assert digest(data) == manifest['files'][name[5:]], name
    config = js('data/configuration.json')
    assert config['apm_commit'] == '381517fda5107fabf98af7801d5a5103f38e230c'
    contract = js('development-contract-v1.json')
    assert contract['function']['signed_gain'] == -6 and contract['function']['output_center_v'] == .5
    nominal = {r['candidate']: r for r in rows('data/nominal.csv')}
    budget = rows('data/source-budget.csv')
    paths = {r['candidate']: r for r in rows('data/sensitivities.csv')}
    designs = {}
    inventories = {}
    for name, count in [('A37', 2), ('A75', 5)]:
        c = config['candidates'][name]
        assert digest(selected['circuits/' + name + '.cir']) == c['circuit_sha256']
        units, resistors = circuit_inventory(text('circuits/' + name + '.cir'))
        assert len(units) == len(c['devices']) == count + 3
        by_name = {u['name']: u for u in units}
        roles = {'input': [], 'reference': [], 'load': []}
        for d in c['devices']:
            u = by_name[d['path'][0]]
            close(u['w_m'], d['w_m']); close(u['l_m'], d['l_m'])
            assert u['model'] == 'apm045_vtg_' + ('nmos' if d['polarity'] == 'n' else 'pmos')
            roles[d['uid'].split('/')[-2]].append(u['name'])
        assert {k: len(v) for k, v in roles.items()} == dict(input=2, reference=count, load=1)
        area = math.fsum(u['w_m'] * u['l_m'] * 1e12 for u in units)
        close(area, float(nominal[name]['area_um2']))
        close(resistors['Rbias'], c['parameters']['rbias'])
        close(resistors['Rs'], c['parameters']['rs'])
        close(resistors['Rsrc'], contract['function']['source_r_ohm'])
        close(resistors['Rload'], contract['function']['load_r_ohm'])
        variance = {r['role']: float(r['variance_mV2']) for r in budget if r['candidate'] == name}
        assert set(variance) == set(roles) and all(v > 0 for v in variance.values())
        total = math.fsum(variance.values())
        designs[name] = dict(reference_units=count, total_current_uA=float(nominal[name]['current_uA']),
                            drawn_mos_area_um2=area, Rbias_ohm=resistors['Rbias'],
                            role_variance_mV2=variance,
                            role_sd_mV={k: math.sqrt(v) for k, v in variance.items()},
                            total_sd_mV=math.sqrt(total),
                            role_variance_fraction={k: v / total for k, v in variance.items()},
                            nominal_gain=float(nominal[name]['gain']),
                            bandwidth_MHz=float(nominal[name]['bandwidth_MHz']),
                            injected_current_sensitivity_kohm=float(paths[name]['injection_ohm']) / 1000,
                            supply_sensitivity_v_per_v=float(paths[name]['supply_v_per_v']))
        inventories[name] = (by_name, resistors)
    a, b = designs['A37'], designs['A75']
    # Independent circuit text comparison: all old units and every R except Rbias survive.
    for name, unit in inventories['A37'][0].items():
        assert inventories['A75'][0][name] == unit, name
    assert set(inventories['A75'][0]) - set(inventories['A37'][0]) == {'xr2', 'xr3', 'xr4'}
    for name, value in inventories['A37'][1].items():
        if name != 'Rbias': close(value, inventories['A75'][1][name])
    added_area = math.fsum(inventories['A75'][0][name]['w_m'] * inventories['A75'][0][name]['l_m'] * 1e12
                          for name in ['xr2', 'xr3', 'xr4'])
    close(added_area, b['drawn_mos_area_um2'] - a['drawn_mos_area_um2'])
    fixed_variance = a['role_variance_mV2']['input'] + a['role_variance_mV2']['load']
    n_ratio = a['reference_units'] / b['reference_units']
    scaled_reference = a['role_variance_mV2']['reference'] * n_ratio
    expected = math.sqrt(fixed_variance + scaled_reference)
    floor = math.sqrt(fixed_variance)
    investment = dict(extra_reference_units=3, extra_current_uA=b['total_current_uA']-a['total_current_uA'],
                      current_ratio=b['total_current_uA']/a['total_current_uA'], extra_drawn_area_um2=added_area,
                      inverse_unit_count_ratio=n_ratio, Rbias_ratio=b['Rbias_ohm']/a['Rbias_ohm'],
                      reference_variance_ratio=b['role_variance_mV2']['reference']/a['role_variance_mV2']['reference'],
                      inverse_count_reference_variance_mV2=scaled_reference,
                      inverse_count_total_sd_mV=expected,
                      actual_budget_total_sd_mV=b['total_sd_mV'],
                      inverse_count_total_sd_relative_residual=b['total_sd_mV']/expected-1,
                      total_sd_reduction_mV=a['total_sd_mV']-b['total_sd_mV'],
                      total_sd_reduction_fraction=1-b['total_sd_mV']/a['total_sd_mV'])
    ideal = dict(role='Restricted mathematical limit, not a finite circuit or a new native condition',
                 held_fixed='A37 input and output-load first-order variances and their transmission paths',
                 reference_variance_mV2=0, total_sd_mV=floor,
                 maximum_sd_reduction_fraction=1-floor/a['total_sd_mV'],
                 cost='No finite current/area assigned. Infinite reference averaging would consume unbounded resources.',
                 exclusions='No resizing of load/input, feedback redesign, source correlation, nonlinear tails or calibrated foundry distribution')
    stats = js('data/statistics.json'); table = rows('data/contrasts.csv')
    plan = js('confirmation-plan.json'); assert len(plan['indices']) == 464
    contrasts = []
    for condition in ['nominal', 'supply_097']:
        old = next(r for r in stats['contrasts'] if r['left']=='A37' and r['right']=='A75' and r['condition']==condition)
        row = next(r for r in table if r['left']=='A37' and r['right']=='A75' and r['condition']==condition)
        for k in ['mean_mae_reduction_v', 'meaningful_reduction_v', 'paired_sd_v']: close(float(row[k]),old[k])
        assert int(row['n_pairs']) == old['n_pairs'] == 464
        lo, hi = old['family95_t_interval_v']; close(float(row['ci_low_v']),lo); close(float(row['ci_high_v']),hi)
        close((lo+hi)/2,old['mean_mae_reduction_v'])
        contrasts.append(dict(condition=condition,role='Original known confirmation result, no new statistical fit',
                              n_pairs=464,mae_reduction_mV=1000*old['mean_mae_reduction_v'],
                              interval_mV=[1000*lo,1000*hi],required_reduction_mV=1000*old['meaningful_reduction_v'],
                              interval_lower_exceeds_required=lo>old['meaningful_reduction_v'],
                              interval_supports_positive_reduction=lo>0,
                              interval_contains_required_reduction=lo<=old['meaningful_reduction_v']<=hi))
    result = dict(schema='adn.learning-reference-investment.v1', status='PASS',
                  role='Retrospective arithmetic and restricted first-order limit from the original public study',
                  source=dict(url=SOURCE_URL,sha256=SOURCE_SHA256,
                              selected_member_sha256={n:digest(data) for n,data in selected.items()}),
                  source_profile='Historical native-Jacobian source-transfer budget; independent observable coordinates per unit. Raw DELVTO and ln(MULU0) are mapped jointly, not assumed independent. No recalibration here.',
                  designs=designs, investment=investment, restricted_floor=ideal, original_contrasts=contrasts,
                  limits='Variance/SD are model-derived local spread estimates, not measured MAE, noise, manufacturing yield or global guarantees. All source observations and confirmation exposures remain historical.')
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'reference-investment-summary.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    table=[]
    for name,d in designs.items():
        table.append(dict(case=name,role='original finite design',reference_units=d['reference_units'],
                          total_current_uA=d['total_current_uA'],drawn_mos_area_um2=d['drawn_mos_area_um2'],
                          input_sd_mV=d['role_sd_mV']['input'],load_sd_mV=d['role_sd_mV']['load'],
                          reference_sd_mV=d['role_sd_mV']['reference'],total_sd_mV=d['total_sd_mV']))
    table.append(dict(case='A37 with reference variance removed',role='restricted mathematical floor',reference_units='',
                      total_current_uA='',drawn_mos_area_um2='',input_sd_mV=a['role_sd_mV']['input'],
                      load_sd_mV=a['role_sd_mV']['load'],reference_sd_mV=0,total_sd_mV=floor))
    with (args.output/'reference-investment.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]),lineterminator='\n');w.writeheader();w.writerows(table)
    print(json.dumps(dict(status='PASS',investment=investment,restricted_floor_sd_mV=floor,original_contrasts=contrasts),indent=2))


if __name__ == '__main__':
    main()
