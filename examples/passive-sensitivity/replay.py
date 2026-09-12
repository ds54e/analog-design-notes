"""MIT: replay a saved passive-condition target with public APM v5.0.0.

The nominal body, fixed log-resistance condition and saved MOS raw parameters
are public inputs. No sampler, trim, rebias or recentering is called.
"""
import argparse
import copy
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for name in ['apm','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--binary',type=Path,default=Path('/usr/local/bin/ngspice'))
    p.add_argument('--candidate',choices=['R37','R75'],default='R75')
    p.add_argument('--index',type=int,default=1000)
    p.add_argument('--scenario',choices=['common-0.1','common+0.1','ratio-0.06','ratio+0.06','ratio-0.1','ratio+0.1'],default='ratio+0.1')
    a=p.parse_args();here=Path(__file__).resolve().parent;apm=a.apm.resolve();out=a.output
    assert out.is_absolute() and not out.exists() and not out.resolve().is_relative_to(apm)
    cfg=json.loads((here/'data/configuration.json').read_text());manifest=json.loads((here/'data/manifest.json').read_text())
    for name,digest in manifest['files'].items():assert sha(here/'data'/name)==digest,name
    commit=subprocess.check_output(['git','-C',str(apm),'rev-parse','HEAD'],text=True).strip();assert commit==cfg['apm_commit']
    for name,digest in cfg['model_hashes'].items():assert sha(apm/name)==digest,name
    def read(name):
        with (here/'data'/name).open() as f:return list(csv.DictReader(f))
    expected=next(r for r in read('transfer.csv') if r['candidate']==a.candidate and int(r['index'])==a.index and r['scenario']==a.scenario)
    d=cfg['candidates'][a.candidate];nominal_body=here/'circuits'/(a.candidate+'.cir')
    assert sha(nominal_body)==d['circuit_sha256']
    condition=next(c for c in cfg['scenarios'] if c['id']==a.scenario);z=condition['log_value']
    changes={'Rs':z,'Rd':z} if condition['kind']=='common' else {'Rs':z/2,'Rd':-z/2}
    resistors={name:d['parameters'][name.lower()]*math.exp(delta) for name,delta in changes.items()}
    lines=[];changed=[]
    for line in nominal_body.read_text().splitlines():
        fields=line.split()
        if fields and fields[0] in resistors:
            name=fields[0];assert len(fields)==4;fields[-1]=format(resistors[name],'.17g');line=' '.join(fields);changed.append(name)
        lines.append(line)
    assert set(changed)=={'Rs','Rd'}
    selected=[r for r in read('realizations.csv') if r['candidate']==a.candidate and int(r['index'])==a.index]
    units={r['uid']:r for r in selected};assert set(units)=={u['uid'] for u in d['devices']}
    for unit in d['devices']:
        r=units[unit['uid']]
        assert r['status']=='RESOLVED' and r['polarity']==unit['polarity'] and float(r['w_m'])==unit['w_m'] and float(r['l_m'])==unit['l_m']
        assert r['original_realization_id']==expected['original_realization_id']
    out.mkdir(parents=True);body=out/'circuit.cir';body.write_text('\n'.join(lines)+'\n')
    assert sha(body)==expected['circuit_sha256'],'Changed body must match the original target exactly'
    os.environ['APM_REPO_ROOT']=str(apm);sys.path.insert(0,str(apm/'src'))
    import numpy as np
    from apm import research,research_spice as spice
    request=dict(schema=research.SCHEMAS['request'],circuit=str(body),devices=d['devices'],other_variation_leaves=[],analyses=copy.deepcopy(d['analyses']))
    request['analyses'][0]['vectors']=[v for v in request['analyses'][0]['vectors'] if not v.startswith('@')]
    for recipe in request['analyses']:recipe['set_sources']={}
    binding={};spice.flatten(body,binding)
    for name in spice.MODELS:spice.flatten(apm/name,binding)
    real=research.seal(dict(schema=research.SCHEMAS['realization'],origin='research',profile_tier='SOURCE_TRANSFER_HYPOTHESIS',status='RESOLVED',
        sample_context_id=research.sample_context(request),request_id=research.canonical_hash(request),input_binding=binding,
        devices=[{**u,'raw':[float(units[u['uid']]['delvto_v']),float(units[u['uid']]['ln_mulu0'])]} for u in d['devices']],
        public_replay=dict(original_realization_id=expected['original_realization_id'],scenario=condition,
            rule='Same MOS raw values and geometry; fixed passive target; relocated input bindings only.')))
    (out/'relocated-realization.json').write_text(json.dumps(real,indent=2)+'\n')
    report=spice.execute(apm,a.binary,out/'native',request,real);assert report['status']=='PASS',report['errors']
    folder=Path(report['directory']);op=np.loadtxt(folder/'analysis0.txt',skiprows=1,ndmin=2);ac=np.loadtxt(folder/'analysis1.txt',skiprows=1,ndmin=2)
    assert op.shape[0]==1 and ac.shape[0]==401 and abs(ac[0,0]-1e3)<1e-8 and abs(ac[-1,0]-1e8)<.01
    v=dict(zip(request['analyses'][0]['vectors'],map(float,op[0,1:])))
    for nodes in d['terminals'].values():
        values=[0. if node=='0' else v['v('+node+')'] for node in nodes]
        assert max(values)-min(values)<=1+1e-8
    kcl=max(abs(sum(v['i('+source+')'] for source in group)) for group in d['sensors'].values());assert kcl<1e-11
    # For this resistor-load circuit only, the measured supply current flows
    # through Rd, and the sum of measured source-terminal currents flows in Rs.
    # Thus actual terminal observations provide an Ohm-law readback without
    # extending the pinned public recipe's terminal-vector API.
    rs_current=-sum(v['i('+group[2]+')'] for group in d['sensors'].values())
    readback=dict(Rs=v['v(source)']/rs_current,Rd=(v['v(vdd)']-v['v(out)'])/(-v['i(Vdd)']))
    assert all(math.isclose(readback[name],value,rel_tol=2e-7) for name,value in resistors.items()),readback
    out_kcl=-v['i(Vdd)']-sum(v['i('+group[0]+')'] for group in d['sensors'].values())-v['v(out)']/d['parameters']['rload']
    assert abs(out_kcl)<1e-11
    h=(ac[:,5]+1j*ac[:,6])/(ac[:,1]+1j*ac[:,2])
    measured=dict(out_v=v['v(out)'],gain=float(h[0].real),current_a=max(0,-v['i(Vdd)'])+max(0,-v['i(Vin)']))
    differences={key:value-float(expected[key]) for key,value in measured.items()}
    assert abs(differences['out_v'])<1e-6 and abs(differences['gain'])<1e-3 and abs(differences['current_a'])<1e-9,differences
    result=dict(status='PASS',candidate=a.candidate,index=a.index,scenario=a.scenario,measured=measured,differences_from_published=differences,
        reproduction_status='PASS',dc_ac_current_specification_pass=expected['dc_ac_current_pass']=='True',
        MOS_readback='Pinned APM actual model/W/L/m/nf/DELVTO/MULU0 before/applied/after OP/AC.',
        resistor_readback_ohm=readback,resistor_readback_method='Actual source/drain supply port currents and terminal voltage drops, using circuit KCL.',
        max_device_kcl_residual_a=kcl,out_kcl_residual_a=out_kcl,apm_commit=commit,
        scope='One known MOS realization under a saved synthetic passive condition. No new sampling, sine/noise or full campaign replay.')
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
