"""MIT: learning.initial-d1-interface.measure.v1, saved observation analysis.

Exact finite plan: same initial six-MOS core, three functional attempts and
six C-fixture OP clamps. No native execution, sampling or optimization.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from design import sha

ROOT=Path(__file__).resolve().parent


def read(path): return json.loads(path.read_text())


def inputs():
    data=ROOT/'data'; cfg=read(data/'interface-inputs.json'); plan=read(ROOT/'receiving-port-plan.json')
    assert sha(ROOT/'receiving-port-plan.json')==cfg['plan_sha256']
    for name, expected in read(data/'interface-manifest.json')['files'].items(): assert sha(data/name)==expected, name
    arrays={}
    for name, c in cfg['cases'].items():
        assert sha(data/c['array_file'])==c['array_sha256']
        assert sha(ROOT/'circuits'/c['circuit_file'])==c['circuit_sha256']
        assert c['physical_devices']==plan['physical_devices']
        with np.load(data/c['array_file'], allow_pickle=False) as z:
            arrays[name]={k:z[k] for k in z.files}
    return cfg,plan,arrays


def columns(case, arrays, index):
    desc=case['analyses'][index]; r=desc['recipe']; a=arrays['analysis'+str(index)]
    assert a.shape==(desc['rows'], desc['columns']) and np.isfinite(a).all()
    values=[a[:,1+2*i]+1j*a[:,2+2*i] if r['kind']=='ac' else a[:,1+i] for i in range(len(r['vectors']))]
    return a[:,0],dict(zip(r['vectors'],values))


def integrate(t, current):
    return np.r_[0., np.cumsum(np.diff(t)*(current[1:]+current[:-1])/2)]


def final_entry(t,y,target,band,origin,horizon):
    assert t[0]<=origin<horizon<=t[-1]*(1+1e-12)
    mask=(t>=origin)&(t<=horizon); x=t[mask]; e=np.abs(y[mask]-target)-band
    assert len(x)>2
    outside=np.flatnonzero(e>0)
    common=dict(target_v=target,band_v=band,origin_s=origin,horizon_s=horizon,observed_until_s=float(x[-1]-origin),
                final_observed_error_v=float(y[mask][-1]-target))
    if len(outside) and outside[-1]==len(x)-1:
        return dict(**common,status='RIGHT_CENSORED',time_s=None,bracket_s=None)
    if not len(outside):
        return dict(**common,status='IN_BAND_AT_FIRST_OBSERVATION',time_s=float(x[0]-origin),bracket_s=[0.,float(x[0]-origin)])
    j=int(outside[-1]);return dict(**common,status='OBSERVED_FINAL_ENTRY',time_s=float(x[j+1]-origin),bracket_s=[float(x[j]-origin),float(x[j+1]-origin)])


def bandwidth(f,h):
    """First-edition bridge definition: -3.000 dB, linear dB/log-f crossing."""
    db=20*np.log10(abs(h)/abs(h[0])); crossing=np.flatnonzero((db[:-1]>-3)&(db[1:]<=-3))
    if not len(crossing): return dict(status='NO_CROSSING_IN_RANGE',hz=None)
    i=int(crossing[0]);value=np.exp(np.log(f[i])+(-3-db[i])/(db[i+1]-db[i])*np.log(f[i+1]/f[i]))
    return dict(status='CROSSING',hz=float(value),bracket_hz=[float(f[i]),float(f[i+1])],
                definition='First downward -3.000-dB crossing relative to |H(1 kHz)|; dB interpolation in log frequency, as learning.bridge.v1.')


def measure_case(info, arrays, plan):
    case=info['case']; tol=plan['measurements']; maximum_kcl=0.;maximum_span=0.;passive={}
    all_values=[];axes=[]
    for i,desc in enumerate(info['analyses']):
        t,v=columns(info,arrays,i); axes.append(t);all_values.append(v)
        if desc['recipe']['kind']=='ac': continue
        for conn in info['connections'].values():
            maximum_kcl=max(maximum_kcl,float(max(abs(sum(v['i('+s+')'] for s in conn['sensors'].values())))))
            voltages=[np.zeros(len(t)) if n=='0' else v['v('+n+')'] for n in conn['nodes'].values()]
            maximum_span=max(maximum_span,float(np.ptp(np.stack(voltages),axis=0).max()))
        if desc['recipe']['kind']=='op':
            expected={'@cload[capacitance]':case['cload_f'], '@rload[resistance]':plan['fixed']['rload_ohm'],
                      '@rbias[resistance]':plan['fixed']['rbias_ohm']}
            if case['fixture']=='D1':expected['@rsource[resistance]']=plan['fixed']['rsource_functional_ohm']
            for name, nominal in expected.items():
                actual=float(v[name][0]);assert abs(actual/nominal-1)<=tol['passive_readback_relative_tolerance']
                passive[name]=actual
    assert maximum_kcl<=tol['kcl_max_a'] and maximum_span<=tol['terminal_span_max_v']
    result=dict(physical_units=len(info['physical_devices']), maximum_MOS_KCL_residual_a=maximum_kcl,
                maximum_terminal_span_v=maximum_span, passive_parameter_readback=passive)
    def available(v): return -v['i(VtM2d)']-v['i(VtM4d)']
    if case['fixture']=='C':
        v=all_values[0]; op={n:float(a[0]) for n,a in v.items()};ia=float(available(v)[0])
        ir=op['v(out)']/plan['fixed']['rload_ohm'];residual=ia-ir-op['i(Vout)']
        assert abs(residual)<=tol['kcl_max_a']
        result.update(op=op,available_current_a=ia,load_current_a=ir,net_current_after_resistor_a=ia-ir,
            clamp_sink_current_a=op['i(Vout)'],clamp_output_kcl_residual_a=residual,
            source_branch_current_a=-op['i(VtM2d)'],sink_branch_current_a=op['i(VtM4d)'],
            source_PMOS_vsd_v=op['v(tail)']-op['v(out)'],sink_NMOS_vds_v=op['v(out)'])
        return result
    rload=plan['fixed']['rload_ohm']; rsource=plan['fixed']['rsource_functional_ohm'];c=case['cload_f']
    ops=[]
    for i in range(3):
        v=all_values[i];op={n:float(x[0]) for n,x in v.items()}
        kcl=float((available(v)-v['v(out)']/rload+v['i(Vfeedback)'])[0])
        input_kcl=float(((v['v(source)']-v['v(inp)'])/rsource+v['i(Vin)'])[0])
        assert max(abs(kcl),abs(input_kcl))<=tol['kcl_max_a']
        ops.append(dict(input_source_v=op['v(source)'],input_node_v=op['v(inp)'],output_v=op['v(out)'],
            error_v=op['v(out)']-op['v(source)'],VDD_delivered_a=-op['i(Vdd)'],input_delivered_a=-op['i(Vin)'],
            output_kcl_residual_a=kcl,input_resistor_residual_a=input_kcl,op=op))
    result['dc']=ops
    f=axes[3];v=all_values[3];h=v['v(out)']/v['v(source)']
    assert len(f)==361 and np.all(np.diff(f)>0)
    result['ac']=dict(gain_1khz=float(abs(h[0])),bandwidth=bandwidth(f,h),peaking_db=float(20*np.log10(max(abs(h))/abs(h[0]))))
    # Independent current/admittance route reads the load in the solved AC circuit.
    ic=available(v)-v['v(out)']/rload+v['i(Vfeedback)'];derived_c=ic/(2j*np.pi*f*v['v(out)'])
    select=(f>=1e5)&(f<=1e7)
    result['ac']['maximum_inferred_capacitance_relative_error']=float(max(abs(derived_c[select]/c-1)))
    assert result['ac']['maximum_inferred_capacitance_relative_error']<1e-5
    t=axes[4];v=all_values[4];vo=v['v(out)'];vin=v['v(source)']
    assert t[0]==0 and abs(t[-1]-case['stop_s'])<1e-12 and np.all(np.diff(t)>0)
    assert max(np.diff(t))<=case['tmax_s']*(1+1e-6)
    pwl=np.array(case['pwl']);source_error=float(max(abs(vin-np.interp(t,pwl[:,0],pwl[:,1]))))
    assert source_error<=tol['pwl_max_error_v']
    net=available(v)-vo/rload+v['i(Vfeedback)'];charge=integrate(t,net);expected_charge=c*(vo-vo[0])
    charge_error=float(max(abs(charge-expected_charge))); excursion=c*float(np.ptp(vo));relative=charge_error/excursion
    assert relative<=tol['charge_relative_residual_max']
    source_residual=float(max(abs((vin-v['v(inp)'])/rsource+v['i(Vin)'])))
    assert source_residual<=tol['kcl_max_a']
    low,high=ops[0]['output_v'],ops[2]['output_v'];band=tol['step_band_fraction_of_independent_op_output_difference']*abs(high-low)
    entries={}
    for key,start,end,target,horizon in [('rise',.5e-6,.5125e-6,high,4.5e-6),('fall',4.5e-6,4.5125e-6,low,8.5e-6)]:
        origin=(start+end)/2
        entries[key]=final_entry(t,vo,target,band,origin,horizon)
    ports={}
    for source,node in [('Vdd','vdd'),('Vin','source')]:
        current=-v['i('+source+')'];power=current*v['v('+node+')']
        ports[source]=dict(min_delivered_a=float(min(current)),max_delivered_a=float(max(current)),
            signed_energy_j=float(np.trapz(power,t)),positive_delivered_energy_j=float(np.trapz(np.maximum(power,0),t)),
            absorbed_energy_j=float(np.trapz(np.maximum(-power,0),t)))
    result['transient']=dict(rows=len(t),tmax_s=case['tmax_s'],output_range_v=[float(min(vo)),float(max(vo))],
        source_pwl_max_error_v=source_error,input_resistor_max_residual_a=source_residual,
        available_current_range_a=[float(min(available(v))),float(max(available(v)))],
        net_capacitor_current_range_a=[float(min(net)),float(max(net))],
        observed_capacitor_charge_excursion_c=excursion,charge_identity_max_residual_c=charge_error,
        charge_identity_relative_residual=relative,step_band_v=band,entries=entries,ports=ports)
    return result


def calculate():
    cfg,plan,arrays=inputs();result=dict(schema=plan['measurements']['version'],plan_identity=plan['identity'],
        status='DECLARED_OBSERVATIONS_ANALYZED',cases={},comparison={})
    for name,info in cfg['cases'].items():result['cases'][name]=measure_case(info,arrays[name],plan)
    cases=result['cases'];a=cases['D1-5p'];b=cases['D1-25p'];ref=cases['D1-25p-refined']
    delta=max(abs(x['output_v']-y['output_v']) for x,y in zip(a['dc'],b['dc']))
    assert delta<=plan['measurements']['dc_load_difference_tolerance_v']
    checks={}
    for edge in ['rise','fall']:
        primary=b['transient']['entries'][edge];refined=ref['transient']['entries'][edge]
        same=primary['status']==refined['status'];difference=None;accepted=same
        if primary['time_s'] is not None and refined['time_s'] is not None:
            difference=abs(primary['time_s']-refined['time_s'])
            t=plan['measurements']['refinement_time_tolerance'];accepted=same and difference<=max(t['absolute_s'],t['relative']*primary['time_s'])
        checks[edge]=dict(status='AGREES_WITHIN_DECLARED_TOLERANCE' if accepted else 'UNRESOLVED_NUMERICAL_DISAGREEMENT',
                         time_difference_s=difference,primary=primary,refined=refined)
    result['comparison']=dict(maximum_5p_25p_dc_output_difference_v=delta,
        bandwidth_ratio_25p_to_5p=b['ac']['bandwidth']['hz']/a['ac']['bandwidth']['hz'],
        settling_ratio_25p_to_5p={k:b['transient']['entries'][k]['time_s']/a['transient']['entries'][k]['time_s']
                                if b['transient']['entries'][k]['time_s'] is not None and a['transient']['entries'][k]['time_s'] else None for k in ['rise','fall']},
        timestep_control=checks,prospective_estimates=plan['prospective_estimates'])
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists(), 'Keep prior analysis; use a new destination'
    result=calculate();a.output.mkdir(parents=True)
    (a.output/'interface-summary.json').write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    rows=[]
    for name,value in result['cases'].items():
        if name.startswith('C-'):
            rows.append(dict(case=name,input_difference_v=value['op']['v(inp)']-value['op']['v(inn)'],output_v=value['op']['v(out)'],
                available_a=value['available_current_a'],net_after_resistor_a=value['net_current_after_resistor_a'],
                source_branch_a=value['source_branch_current_a'],sink_branch_a=value['sink_branch_current_a'],
                source_PMOS_vsd_v=value['source_PMOS_vsd_v'],sink_NMOS_vds_v=value['sink_NMOS_vds_v']))
    with (a.output/'clamps.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(status=result['status'],comparison=result['comparison']),indent=2))


if __name__=='__main__':main()
