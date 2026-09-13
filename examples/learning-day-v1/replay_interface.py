"""MIT: one public native replay of initial D1 with its new 25-pF load.

Uses the pinned public APM API and saved physical units/circuit. No private
adapter, legacy runtime, sampling, optimization or site-build acquisition.
"""
import argparse
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from interface import inputs, columns, final_entry, bandwidth, integrate
from design import sha


def main():
    p=argparse.ArgumentParser();p.add_argument('--apm',type=Path,required=True)
    p.add_argument('--binary',type=Path,default=Path('/usr/local/bin/ngspice'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    apm=a.apm.resolve();out=a.output
    assert out.is_absolute() and not out.exists() and not out.resolve().is_relative_to(apm)
    cfg,plan,saved=inputs();name='D1-25p';case=cfg['cases'][name];condition=case['case']
    assert subprocess.check_output(['git','-C',str(apm),'rev-parse','HEAD'],text=True).strip()==cfg['apm_commit']
    for path,h in cfg['model_hashes'].items():assert sha(apm/path)==h,path
    os.environ['APM_REPO_ROOT']=str(apm);sys.path.insert(0,str(apm/'src'))
    from apm import research,research_spice as spice
    here=Path(__file__).resolve().parent;body=here/'circuits'/case['circuit_file']
    text=body.read_text();recipes=[]
    for desc in case['analyses']:
        recipe=copy.deepcopy(desc['recipe']);recipe['vectors']=[v for v in recipe['vectors'] if not v.startswith('@')]
        # The public API has no generic alter/parameter-vector adapter. These
        # fixed values are already present in this exact circuit body, and
        # terminal-current V/I and AC/charge checks below read their effects.
        for element,value in recipe.pop('elements',{}).items():
            line=next(x for x in text.splitlines() if x.startswith(element+' '))
            assert float(line.split()[-1])==value
        for source,value in recipe.pop('ac_sources',{}).items():
            line=next(x for x in text.splitlines() if x.startswith(source+' ')).split()
            assert float(line[line.index('AC')+1])==value
        if recipe['kind']=='tran':
            assert min(recipe['step'],recipe['stop']/50)==recipe['tmax']
        recipes.append(recipe)
    request=dict(schema=research.SCHEMAS['request'],circuit=str(body),devices=case['devices'],
                 other_variation_leaves=[],analyses=recipes)
    binding={};spice.flatten(body,binding)
    for path in spice.MODELS:spice.flatten(apm/path,binding)
    real=research.seal(dict(schema=research.SCHEMAS['realization'],origin='research',profile_tier='ARTIFICIAL',status='RESOLVED',
        sample_context_id=research.sample_context(request),request_id=research.canonical_hash(request),input_binding=binding,
        devices=case['physical_devices'],public_replay=dict(original_realization_id=case['original_realization_id'],
            case=name,rule='Relocated paths and public terminal-vector recipes; exact saved six zero-raw units and 25-pF circuit, no redraw or trim.')))
    out.mkdir(parents=True);(out/'relocated-realization.json').write_text(json.dumps(real,indent=2)+'\n')
    report=spice.execute(apm,a.binary,out/'native',request,real,temperature_c=cfg['temperature_c'])
    assert report['status']=='PASS',report['errors']
    folder=Path(report['directory']);axes=[];values=[];differences={};kcl=0.;span=0.
    for i,r in enumerate(recipes):
        data=np.loadtxt(folder/f'analysis{i}.txt',skiprows=1,ndmin=2)
        assert np.isfinite(data).all() and data.shape[1]==1+len(r['vectors'])*(2 if r['kind']=='ac' else 1)
        t=data[:,0];v={key:(data[:,1+2*j]+1j*data[:,2+2*j] if r['kind']=='ac' else data[:,1+j]) for j,key in enumerate(r['vectors'])}
        old_t,old=columns(case,saved[name],i)
        if r['kind']=='tran':
            assert t[0]==0 and abs(t[-1]-condition['stop_s'])<1e-12 and np.all(np.diff(t)>0)
            assert max(np.diff(t))<=condition['tmax_s']*(1+1e-6)
            delta=float(max(abs(np.interp(old_t,t,v['v(out)'])-old['v(out)'])))
        else:
            assert len(t)==len(old_t) and np.allclose(t,old_t,rtol=1e-11,atol=1e-12)
            delta=float(max(abs(v['v(out)']-old['v(out)'])))
        assert delta<(1e-3 if r['kind']=='ac' else 1e-6),(i,delta)
        differences[str(i)]=dict(kind=r['kind'],rows=len(t),maximum_output_absolute_difference=delta,
            units='V/V for unit AC; V for OP/transient. Adaptive transient values linearly compared at saved time coordinates.')
        if r['kind']!='ac':
            for conn in case['connections'].values():
                kcl=max(kcl,float(max(abs(sum(v['i('+s+')'] for s in conn['sensors'].values())))))
                voltages=[np.zeros(len(t)) if n=='0' else v['v('+n+')'] for n in conn['nodes'].values()]
                span=max(span,float(np.ptp(np.stack(voltages),axis=0).max()))
        axes.append(t);values.append(v)
    assert kcl<=1e-10 and span<=1+1e-8
    resistor_readback={};v=values[1]
    flows=dict(Rsource=((v['v(source)']-v['v(inp)'])[0],-v['i(Vin)'][0],1e4),
        Rload=(v['v(out)'][0],(-v['i(VtM2d)']-v['i(VtM4d)']+v['i(Vfeedback)'])[0],1e6),
        Rbias=(v['v(bias)'][0],(-v['i(VtM6d)']-v['i(VtM6g)']-v['i(VtM5g)'])[0],plan['fixed']['rbias_ohm']))
    for key,(drop,current,expected) in flows.items():
        assert abs(current)>1e-12
        observed=float(drop/current);error=observed/expected-1;assert abs(error)<1e-5
        resistor_readback[key]=dict(observed_v_over_i_ohm=observed,expected_ohm=expected,relative_difference=error)
    f=axes[3];v=values[3];h=v['v(out)']/v['v(source)'];c=condition['cload_f']
    ic=-v['i(VtM2d)']-v['i(VtM4d)']-v['v(out)']/1e6+v['i(Vfeedback)']
    measured_c=ic/(2j*np.pi*f*v['v(out)']);mask=(f>=1e5)&(f<=1e7)
    ac_c_error=float(max(abs(measured_c[mask]/c-1)));assert ac_c_error<1e-5
    t=axes[4];v=values[4];vo=v['v(out)'];pwl=np.array(condition['pwl'])
    assert max(abs(v['v(source)']-np.interp(t,pwl[:,0],pwl[:,1])))<=1e-9
    net=-v['i(VtM2d)']-v['i(VtM4d)']-vo/1e6+v['i(Vfeedback)']
    charge_error=float(max(abs(integrate(t,net)-c*(vo-vo[0])))/(c*np.ptp(vo)));assert charge_error<=.001
    low=float(values[0]['v(out)'][0]);high=float(values[2]['v(out)'][0]);band=.01*abs(high-low);entries={}
    expected=read_summary(here)['cases'][name]
    for edge,origin,target,horizon in [('rise',.50625e-6,high,4.5e-6),('fall',4.50625e-6,low,8.5e-6)]:
        item=final_entry(t,vo,target,band,origin,horizon);old=expected['transient']['entries'][edge]
        assert item['status']==old['status']
        if item['time_s'] is not None:assert abs(item['time_s']-old['time_s'])<=max(4e-9,.01*old['time_s'])
        entries[edge]=item
    result=dict(status='PASS',case=name,apm_commit=cfg['apm_commit'],new_run_id=report['run_id'],
        original_realization_id=case['original_realization_id'],differences=differences,
        maximum_MOS_KCL_residual_a=kcl,maximum_terminal_span_v=span,resistor_readback=resistor_readback,
        capacitor_ac_current_readback_max_relative_error=ac_c_error,capacitor_charge_relative_residual=charge_error,
        bandwidth=bandwidth(f,h),settling=entries,
        MOS_readback='Pinned APM verifies actual model/W/L/m/nf/DELVTO/MULU0 before/applied/after each analysis.',
        passive_readback='Public API: observed DC V/I for three resistors; complex AC current/voltage and transient charge for Cload. Native parameter vectors are not exposed by this replay.',
        limits='One nominal 25-pF initial-D1 OP/AC/standard-step replay. No Local cohort, full pass-gate charge, LDO, old selected-D1 transfer or independent scientific replication.')
    (out/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n');print(json.dumps(result,indent=2))


def read_summary(here):return json.loads((here/'interface-summary.json').read_text())


if __name__=='__main__':main()
