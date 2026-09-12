"""MIT: recompute finite-cohort passive-scenario summaries from public CSV.

This is an exact conditional description of a known MOS cohort, not a fresh
population inference. Each MOS index contributes both signs with equal weight.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def read(folder, name):
    with (folder/(name+'.csv')).open() as f: return list(csv.DictReader(f))


def mean(values):
    return math.fsum(values)/len(values)


def moments(values):
    avg = mean(values); mse = mean([x*x for x in values])
    return dict(mean_v=avg, mae_v=mean([abs(x) for x in values]), mse_v2=mse,
                rms_v=math.sqrt(mse), population_sd_v=math.sqrt(mean([(x-avg)**2 for x in values])))


def verify_gradients(data):
    raw=read(data,'gradient-observations');slopes=read(data,'sensitivities')
    cfg=json.loads((data/'configuration.json').read_text())
    op_rows=read(data,'nominal-op');results=[]
    for n in ['R37','A37','R75','A75']:
        d=cfg['candidates'][n];p=d['parameters'];op={r['vector']:float(r['value']) for r in op_rows if r['candidate']==n}
        def total(role,k):
            return math.fsum(op[f'@m.{unit["path"][0]}.mapm045_vtg_core[{k}]'] for unit in d['devices'] if unit['uid'].split('/')[-2]==role)
        g=total('input','gm')+total('input','gmbs')+total('input','gds');D=1+p['rs']*g
        go=1/p['rload']+(1/p['rd'] if p.get('rd') else 0)+total('load','gds')+total('input','gds')/D
        predicted=dict(rs=op['v(source)']/p['rs']*(D-1)/D/go,rload=op['v(out)']/p['rload']/go,
            rsrc=total('input','gm')/D/go*(op['v(in)']-op['v(gate)']))
        if p['load']=='resistor':predicted['rd']=-(op['v(vdd)']-op['v(out)'])/p['rd']/go
        else:predicted['rbias']=-total('load','gm')/go*op['v(pbias)']/(1+p['rbias']*(total('reference','gm')+total('reference','gds')))
        for r in [r for r in slopes if r['candidate']==n]:
            parameter=r['parameter'];points={float(x['log_step']):float(x['out_v']) for x in raw if x['candidate']==n and x['parameter']==parameter}
            derivative=(points[.0005]-points[-.0005])/.001
            assert math.isclose(derivative,float(r['observed_v_per_logR']),rel_tol=1e-12)
            assert math.isclose(predicted[parameter],float(r['estimate_v_per_logR']),rel_tol=1e-8,abs_tol=2e-13)
            results.append(dict(candidate=n,parameter=parameter,native_centered_derivative_v=derivative,equation_estimate_v=predicted[parameter]))
    return results


def calculate(data):
    rows = read(data, 'transfer'); base = read(data, 'known-baselines')
    slopes = read(data, 'sensitivities')
    nominal = {n: {r['parameter']: float(r['observed_v_per_logR']) for r in slopes if r['candidate']==n}
               for n in ['R37', 'R75']}
    known = {(r['candidate'], int(r['index'])): r for r in base}
    output = dict(schema='adn.passive-statistics.v1', exposure='Known MOS cohort, synthetic two-point log-passive scenarios; no population confirmation or yield claim.',
                  weighting='Both passive signs equally weighted within each of 464 retained MOS indices.', groups=[],
                  contrasts=[], curvature=[], baseline={}, slope_moments={}, gradient_recomputation=verify_gradients(data))
    for n in ['R37', 'R75']:
        b = [r for r in base if r['candidate']==n]
        assert len(b)==464 and {int(r['index']) for r in b}==set(range(1000,1464))
        output['baseline'][n] = moments([float(r['error_v']) for r in b])
        s = nominal[n]
        output['slope_moments'][n] = dict(nominal_common_v=s['rs']+s['rd'], nominal_ratio_v=(s['rs']-s['rd'])/2,
            known_OP_common_mean_square_v2=mean([float(r['s_common_v'])**2 for r in b]),
            known_OP_ratio_mean_square_v2=mean([float(r['s_ratio_v'])**2 for r in b]))
        for kind, level in [('common', .10), ('ratio', .06), ('ratio', .10)]:
            group = [r for r in rows if r['candidate']==n and r['kind']==kind and abs(float(r['log_value']))==level]
            assert len(group)==928
            good = [r for r in group if r['numerical_usable']=='True' and r['terminal_supported']=='True']
            stats = dict(candidate=n, kind=kind, log_amplitude=level, physical_circuits=464, observations=928,
                         usable=sum(r['numerical_usable']=='True' for r in group),
                         supported=sum(r['terminal_supported']=='True' for r in group),
                         dc_ac_current_pass=sum(r['dc_ac_current_pass']=='True' for r in group),
                         comparison_complete=len(good)==928)
            if good:
                actual = [float(r['error_v']) for r in good]; stats.update(moments(actual))
                for name, column in [('nominal_prediction', 'predicted_nominal_slope_v'), ('known_OP_prediction', 'predicted_known_OP_slope_v')]:
                    pred = [float(r[column])-.5 for r in good]
                    stats[name] = dict(**moments(pred),
                        prediction_error_rms_v=math.sqrt(mean([(x-y)**2 for x,y in zip(actual,pred)])),
                        prediction_error_max_abs_v=max(abs(x-y) for x,y in zip(actual,pred)))
            output['groups'].append(stats)
            if len(good)==928:
                paired = {}
                for r in good: paired.setdefault(int(r['index']), {})[float(r['log_value'])] = float(r['error_v'])
                baseline_sq=[]; odd_sq=[]; offset_curvature=[]; curvature_sq=[]; shifts=[]
                for index, signs in paired.items():
                    e0=float(known[n,index]['error_v']); mid=(signs[level]+signs[-level])/2
                    odd=(signs[level]-signs[-level])/2; shift=mid-e0
                    baseline_sq.append(e0*e0); odd_sq.append(odd*odd)
                    offset_curvature.append(2*e0*shift); curvature_sq.append(shift*shift); shifts.append(shift)
                terms=dict(baseline_v2=mean(baseline_sq), odd_response_v2=mean(odd_sq),
                           offset_even_response_v2=mean(offset_curvature), even_response_square_v2=mean(curvature_sq))
                assert math.isclose(math.fsum(terms.values()), stats['mse_v2'], rel_tol=2e-14)
                output['curvature'].append(dict(candidate=n, kind=kind, log_amplitude=level,
                    mean_even_shift_v=mean(shifts), **terms,
                    role='Retrospective exact two-sign decomposition. The even component is not a prospective higher-order predictor.'))
    for kind, level in [('common', .10), ('ratio', .06), ('ratio', .10)]:
        pair = [next(g for g in output['groups'] if g['candidate']==n and g['kind']==kind and g['log_amplitude']==level) for n in ['R37','R75']]
        if all(g['comparison_complete'] for g in pair):
            output['contrasts'].append(dict(kind=kind, log_amplitude=level,
                R75_minus_R37_mse_v2=pair[1]['mse_v2']-pair[0]['mse_v2'],
                R75_minus_R37_rms_v=pair[1]['rms_v']-pair[0]['rms_v'],
                R75_minus_R37_mae_v=pair[1]['mae_v']-pair[0]['mae_v']))
    gap=output['baseline']['R37']['mse_v2']-output['baseline']['R75']['mse_v2']
    sm=output['slope_moments'];den=sm['R75']['nominal_ratio_v']**2-sm['R37']['nominal_ratio_v']**2
    den_op=sm['R75']['known_OP_ratio_mean_square_v2']-sm['R37']['known_OP_ratio_mean_square_v2']
    output['linear_ratio_break_even_log_amplitude']=math.sqrt(gap/den)
    output['known_OP_linear_ratio_break_even_log_amplitude']=math.sqrt(gap/den_op)
    return output


def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=calculate(a.data)
    with a.output.open('x') as f: json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status='RECOMPUTED',groups=len(result['groups']),contrasts=result['contrasts'],
        linear_ratio_break_even_log_amplitude=result['linear_ratio_break_even_log_amplitude']),indent=2))


if __name__=='__main__': main()
