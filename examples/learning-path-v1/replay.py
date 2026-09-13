"""MIT: one finite public replay of the saved initial-D1 OP, DC and central AC.

Only path bindings and analysis selection change. Exact circuit bytes and six
saved zero-raw physical units are used. No calibration, sampling or legacy
campaign module is imported. This does not replay the archived loop/noise/step
analyses; analyze.py recomputes the adopted saved observations separately.
"""
import argparse
import copy
import json
import os
import subprocess
import sys
from pathlib import Path
from analyze import inputs, column


def main():
    p=argparse.ArgumentParser();p.add_argument('--apm',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--binary',type=Path,default=Path('/usr/local/bin/ngspice'))
    a=p.parse_args();here=Path(__file__).resolve().parent;apm=a.apm.resolve();out=a.output
    assert out.is_absolute() and not out.exists() and not out.resolve().is_relative_to(apm)
    cfg,arrays=inputs(here);case=cfg['cases']['D1']
    commit=subprocess.check_output(['git','-C',str(apm),'rev-parse','HEAD'],text=True).strip()
    assert commit==cfg['apm_commit']
    from analyze import sha
    for name,h in case['model_hashes'].items():assert sha(apm/name)==h,name
    os.environ['APM_REPO_ROOT']=str(apm);sys.path.insert(0,str(apm/'src'))
    import numpy as np
    from apm import research,research_spice as spice
    selected=[0,1,3];recipes=[]
    for i in selected:
        recipe=copy.deepcopy(next(x['recipe'] for x in case['analyses'] if x['analysis_index']==i))
        # The public API exposes terminal vectors, not native device diagnostics.
        recipe['vectors']=[v for v in recipe['vectors'] if not v.startswith('@')]
        ac=recipe.pop('set_ac_sources',{})
        # These values already appear in the immutable functional SPICE body.
        assert all(v==(1 if k=='Vin' else 0) for k,v in ac.items()),ac
        recipes.append(recipe)
    devices=[{k:v for k,v in d.items() if k!='raw'} for d in case['devices']]
    body=here/'circuits'/(case['stem']+'.cir')
    request=dict(schema=research.SCHEMAS['request'],circuit=str(body),devices=devices,
                 other_variation_leaves=[],analyses=recipes)
    binding={};spice.flatten(body,binding)
    for path in spice.MODELS:spice.flatten(apm/path,binding)
    real=research.seal(dict(schema=research.SCHEMAS['realization'],origin='research',
        profile_tier='ARTIFICIAL',status='RESOLVED',sample_context_id=research.sample_context(request),
        request_id=research.canonical_hash(request),input_binding=binding,devices=case['devices'],
        public_replay=dict(original_realization_id=case['original_realization_id'],
            original_run_id=case['run_id'],original_analysis_indices=selected,
            note='Relocated path bindings and selected terminal-vector recipes; same saved six zero-raw units and exact functional circuit.')))
    out.mkdir(parents=True);(out/'relocated-realization.json').write_text(json.dumps(real,indent=2)+'\n')
    report=spice.execute(apm,a.binary,out/'native',request,real,temperature_c=case['temperature_c'])
    assert report['status']=='PASS',report['errors']
    folder=Path(report['directory']);differences={};maximum_kcl=0.;maximum_span=0.
    for ni,oi in enumerate(selected):
        obs=np.loadtxt(folder/f'analysis{ni}.txt',skiprows=1,ndmin=2)
        recipe=recipes[ni];vec=recipe['vectors'];old=arrays['D1']['analysis'+str(oi)]
        assert len(obs)==len(old) and np.isfinite(obs).all()
        assert np.allclose(obs[:,0],old[:,0],rtol=1e-11,atol=1e-12)
        def col(name):
            j=vec.index(name)
            return obs[:,1+2*j]+1j*obs[:,2+2*j] if recipe['kind']=='ac' else obs[:,1+j]
        delta=float(max(abs(col('v(out)')-column(cfg,arrays,'D1',oi,'v(out)'))))
        assert delta<(1e-3 if recipe['kind']=='ac' else 1e-6),(oi,delta)
        differences[str(oi)]=dict(kind=recipe['kind'],rows=len(obs),maximum_output_absolute_difference=delta,
            units='V/V for unit AC stimulus; V for OP/DC')
        if recipe['kind'] in ('op','dc'):
            for d in case['connections'].values():
                residual=float(max(abs(sum(col('i('+s+')') for s in d['sensors'].values()))))
                voltage=[np.zeros(len(obs)) if n=='0' else col('v('+n+')') for n in d['nodes'].values()]
                maximum_kcl=max(maximum_kcl,residual)
                maximum_span=max(maximum_span,float(max(np.max(voltage,axis=0)-np.min(voltage,axis=0))))
        if oi==0:
            center=dict(output_V=float(col('v(out)')[0]),input_node_V=float(col('v(inp)')[0]),
                VDD_delivered_A=float(-col('i(Vdd)')[0]),input_delivered_A=float(-col('i(Vin)')[0]))
    assert maximum_kcl<1e-10 and maximum_span<=1+1e-8
    result=dict(status='PASS',scope='Initial-D1 OP/DC/central AC replay through the pinned public APM API. No new physical samples, no archived loop/noise/step native replay.',
        apm_commit=commit,original_run_id=case['run_id'],new_run_id=report['run_id'],center=center,
        differences=differences,maximum_MOS_KCL_residual_A=maximum_kcl,maximum_terminal_span_V=maximum_span,
        readback='Pinned APM verified actual model/W/L/m/nf/DELVTO/MULU0 before/applied/after every analysis.')
    (out/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
