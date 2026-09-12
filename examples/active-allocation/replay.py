"""MIT: replay a published saved Local OP/AC case using the pinned public APM.

Only path bindings are new. Physical UID/geometry/raw values come directly
from the public CSV; no sampler, calibration, trim or search is called.
"""
import argparse,csv,hashlib,json,os,subprocess,sys
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--apm',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--binary',type=Path,default=Path('/usr/local/bin/ngspice'))
    p.add_argument('--candidate',choices=['A37','R37','AL37W','AR37','AI37','AB37'],default='AL37W')
    p.add_argument('--index',type=int,default=100);p.add_argument('--condition',choices=['nominal','supply_097'],default='nominal')
    a=p.parse_args();here=Path(__file__).resolve().parent;apm=a.apm.resolve();out=a.output
    if not out.is_absolute() or out.exists() or out.resolve().is_relative_to(apm):raise ValueError('Select a new absolute output outside the model checkout')
    cfg=json.loads((here/'data/configuration.json').read_text());manifest=json.loads((here/'data/manifest.json').read_text())
    for n,h in manifest['files'].items():assert sha(here/'data'/n)==h,('public data changed',n)
    commit=subprocess.check_output(['git','-C',str(apm),'rev-parse','HEAD'],text=True).strip();assert commit==cfg['apm_commit']
    for n,h in cfg['model_hashes'].items():assert sha(apm/n)==h,('model changed',n)
    d=cfg['candidates'][a.candidate];body=here/'circuits'/f'{a.candidate}.cir';assert sha(body)==d['circuit_sha256']
    def read(n):
        with (here/'data'/n).open() as f:return list(csv.DictReader(f))
    selected=[r for r in read('realizations.csv') if r['candidate']==a.candidate and int(r['index'])==a.index]
    assert len(selected)==len(d['devices']) and all(r['status']=='RESOLVED' for r in selected)
    rows={r['uid']:r for r in selected};assert set(rows)=={x['uid'] for x in d['devices']}
    for device in d['devices']:
        r=rows[device['uid']];assert r['polarity']==device['polarity'] and float(r['w_m'])==device['w_m'] and float(r['l_m'])==device['l_m']
    expected=next(r for r in read('confirmation.csv') if r['candidate']==a.candidate and int(r['index'])==a.index and r['condition']==a.condition)
    os.environ['APM_REPO_ROOT']=str(apm);sys.path.insert(0,str(apm/'src'))
    import numpy as np
    from apm import research,research_spice as spice
    condition=next(c for c in cfg['conditions'] if c['id']==a.condition)
    request=dict(schema=research.SCHEMAS['request'],circuit=str(body),devices=d['devices'],other_variation_leaves=[],analyses=d['analyses'])
    # Public APM recipes expose terminal vectors. Native device-parameter
    # diagnostics used in private cause analysis are outside that recipe API.
    request['analyses'][0]['vectors']=[v for v in request['analyses'][0]['vectors'] if not v.startswith('@')]
    for r in request['analyses']:r['set_sources']=condition.get('sources',{})
    binding={};spice.flatten(body,binding)
    for n in spice.MODELS:spice.flatten(apm/n,binding)
    real=research.seal(dict(schema=research.SCHEMAS['realization'],origin='research',profile_tier='SOURCE_TRANSFER_HYPOTHESIS',status='RESOLVED',
        sample_context_id=research.sample_context(request),request_id=research.canonical_hash(request),input_binding=binding,
        devices=[{**x,'raw':[float(rows[x['uid']]['delvto_v']),float(rows[x['uid']]['ln_mulu0'])]} for x in d['devices']],
        public_replay=dict(csv_sha256=sha(here/'data/realizations.csv'),original_realization_id=selected[0]['original_realization_id'],
          seed=cfg['seed'],sample_index=a.index,rule='Relocated bindings; exact saved raw physical units, no resampling or remapping.')))
    out.mkdir(parents=True);(out/'relocated-realization.json').write_text(json.dumps(real,indent=2)+'\n')
    report=spice.execute(apm,a.binary,out/'native',request,real);assert report['status']=='PASS',report['errors']
    folder=Path(report['directory']);op=np.loadtxt(folder/'analysis0.txt',skiprows=1,ndmin=2);ac=np.loadtxt(folder/'analysis1.txt',skiprows=1,ndmin=2)
    assert op.shape[0]==1 and ac.shape[0]==401 and abs(ac[0,0]-1e3)<1e-8 and abs(ac[-1,0]-1e8)<.01
    values=dict(zip(request['analyses'][0]['vectors'],map(float,op[0,1:])))
    for nodes in d['terminals'].values():
        v=[0. if n=='0' else values['v('+n+')'] for n in nodes]
        assert max(v)-min(v)<=1+1e-8
    kcl=max(abs(sum(values['i('+s+')'] for s in group)) for group in d['sensors'].values());assert kcl<1e-11
    h=ac[:,5]+1j*ac[:,6] # AC vectors: in, gate, out, supply current, input current
    vin=ac[:,1]+1j*ac[:,2];h=h/vin
    measured=dict(out_v=values['v(out)'],gain=float(h[0].real),current_a=max(0,-values['i(Vdd)'])+max(0,-values['i(Vin)']))
    errors={k:measured[k]-float(expected[k]) for k in measured}
    assert abs(errors['out_v'])<1e-6 and abs(errors['gain'])<1e-3 and abs(errors['current_a'])<1e-9,errors
    result=dict(status='PASS',candidate=a.candidate,index=a.index,condition=a.condition,measured=measured,differences_from_published=errors,
        reproduction_status='PASS',dc_ac_current_specification_pass=expected['dc_ac_current_pass']=='True',
        readback='Pinned APM checked actual model/W/L/m/nf/DELVTO/MULU0 before/applied/after every OP/AC analysis.',max_device_kcl_residual_a=kcl,
        scope='Selected OP/AC replay with saved Local units. No random redraw, noise, sine, temperature statistics or full-population replay claimed.',run_id=report['run_id'],apm_commit=commit)
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
