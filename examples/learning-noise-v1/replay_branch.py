"""MIT: one finite F75L branch-current diagnostic using pinned public APM."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from analyze import calculate, digest


def main():
    p=argparse.ArgumentParser()
    for n in ['apm','feedback-zip','output']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--binary',type=Path,default=Path('/usr/local/bin/ngspice'));a=p.parse_args()
    root=Path(__file__).resolve().parent;out=a.output;apm=a.apm.resolve()
    assert out.is_absolute() and not out.exists() and not out.resolve().is_relative_to(apm)
    calculate(a.feedback_zip) # verify the selected public inputs and prior finite calculation
    source=json.loads((root/'data/source.json').read_text());plan=json.loads((root/'branch-plan.json').read_text());tol=plan['tolerances']
    assert digest(Path(__file__).read_bytes())==plan['script_sha256']
    assert digest((root/'data/source.json').read_bytes())==plan['source_sha256']
    body=root/plan['circuit'];assert digest(body.read_bytes())==plan['circuit_sha256']
    assert digest((root/plan['forecast']).read_bytes())==plan['forecast_sha256']
    original=next(c for c in source['cases'] if c['case']=='F75L-nominal');assert original['physical_units']==plan['physical_units']
    version=subprocess.check_output([str(a.binary),'--version'],text=True).strip();assert 'ngspice-47' in version
    assert subprocess.check_output(['git','-C',str(apm),'rev-parse','HEAD'],text=True).strip()==source['apm_commit']
    for path,h in source['model_hashes'].items():assert digest((apm/path).read_bytes())==h,path
    fields=[l.split() for l in body.read_text().splitlines() if l.strip() and not l.startswith('*')]
    sensors={q[2]:(q[0],q[1]) for q in fields if q[0].startswith('V') and len(q)==4 and q[3]=='0'}
    conn={q[0]:dict(sensors=[sensors[n][0] for n in q[1:5]],nodes=[sensors[n][1] for n in q[1:5]]) for q in fields if q[0].startswith('x')}
    assert set(conn)=={u['path'][0] for u in plan['physical_units']}
    sources={q[0]:float(q[q.index('AC')+1]) for q in fields if q[0][0].lower() in ['v','i'] and 'AC' in q}
    assert sources=={'Vdd':0,'Vin':0,'Idist':0,'Vbiasdiag':0,'Inoise':plan['conditions']['current_ac_a']}
    vectors=[v for v in original['op_recipe']['vectors'] if not v.startswith('@')]
    vectors += ['v(cap_probe)','v(noise_probe)','i(Vcap)','i(Vnoise)']
    assert len(vectors)==len(set(vectors))
    recipes=[dict(kind='op',vectors=vectors),dict(kind='ac',start=1e3,stop=1e6,points_per_decade=80,vectors=vectors)]
    os.environ['APM_REPO_ROOT']=str(apm);sys.path.insert(0,str(apm/'src'))
    from apm import research,research_spice as spice
    devices=[{k:v for k,v in u.items() if k!='raw'} for u in plan['physical_units']]
    request=dict(schema=research.SCHEMAS['request'],circuit=str(body),devices=devices,other_variation_leaves=[],analyses=recipes)
    binding={};spice.flatten(body,binding)
    for model in spice.MODELS:spice.flatten(apm/model,binding)
    real=research.seal(dict(schema=research.SCHEMAS['realization'],origin='nominal',profile_tier='NOMINAL',status='RESOLVED',
        sample_context_id=research.sample_context(request),request_id=research.canonical_hash(request),input_binding=binding,
        devices=plan['physical_units'],public_replay=dict(original_realization_id=original['realization_id'],rule='Same five saved zero-raw MOS units; new explicit output-to-gate AC current and zero-volt current sensors. No sample draw or rebias.')))
    out.mkdir(parents=True)
    for name,value in [('request.json',request),('realization.json',real),('frozen-plan.json',plan)]:
        (out/name).write_text(json.dumps(value,indent=2)+'\n')
    report=spice.execute(apm,a.binary,out/'native',request,real,temperature_c=plan['conditions']['temperature_c'])
    assert report['status']=='PASS',report
    native=Path(report['directory']);data=[];axes=[]
    for i,recipe in enumerate(recipes):
        z=np.loadtxt(native/f'analysis{i}.txt',skiprows=1,ndmin=2);assert np.isfinite(z).all()
        assert z.shape==(1 if i==0 else 241,1+len(vectors)*(1 if i==0 else 2))
        v={k:z[:,1+j] if i==0 else z[:,1+2*j]+1j*z[:,2+2*j] for j,k in enumerate(vectors)}
        data.append(v);axes.append(z[:,0])
    op,ac=data;f=axes[1];assert np.all(np.diff(f)>0) and abs(f[0]/1e3-1)<1e-12 and abs(f[-1]/1e6-1)<1e-12
    dc_errors={k:float(op[k][0]-original['observed_op'][k]) for k in original['observed_op'] if k.startswith('v(')}
    assert max(abs(v) for v in dc_errors.values())<tol['op_voltage_v']
    roles={u['path'][0]:u['uid'].split('/')[-2] for u in plan['physical_units']}
    def bank(v,role,terminal):return sum(v['i('+c['sensors'][terminal]+')'] for name,c in conn.items() if roles[name]==role)
    dc_kcl=0.;span=0.;ac_kcl=0.
    for c in conn.values():
        dc_kcl=max(dc_kcl,float(abs(sum(op['i('+s+')'][0] for s in c['sensors']))))
        ac_kcl=max(ac_kcl,float(np.max(abs(sum(ac['i('+s+')'] for s in c['sensors'])))))
        volts=[0 if n=='0' else float(op['v('+n+')'][0]) for n in c['nodes']];span=max(span,max(volts)-min(volts))
    assert dc_kcl<tol['dc_kcl_a'] and span<=1+tol['terminal_span_slack_v']
    iin=-op['i(Vin)'][0];ig=float(bank(op,'input',1)[0]);ifr=ig-iin-op['i(Vnoise)'][0]
    out_d=float(bank(op,'input',0)[0]+bank(op,'load',0)[0]);ibias=-float(bank(op,'reference',0)[0]+bank(op,'reference',1)[0]+op['i(Vbiasdiag)'][0])
    params=original['parameters'];observed={'rsrc':float((op['v(in)'][0]-op['v(drive)'][0])/iin),
        'rin':float((op['v(drive)'][0]-op['v(gate)'][0])/iin),'rf':float((op['v(out)'][0]-op['v(gate)'][0])/ifr),
        'rs':float(op['v(source)'][0]/(-bank(op,'input',2)[0])),'rbias':float(op['v(pbias)'][0]/ibias),
        'rload':float(op['v(out)'][0]/(-out_d-ifr-op['i(Vnoise)'][0]-op['i(Vcap)'][0]))}
    relative={k:v/params[k]-1 for k,v in observed.items()};assert max(abs(v) for v in relative.values())<tol['dc_passive_relative']
    current=ac['i(Vnoise)'];source_error=float(max(abs(current/plan['conditions']['current_ac_a']-1)))
    assert source_error<tol['source_ac_relative'] and max(abs(ac['v(in)']))<1e-15 and max(abs(ac['v(vdd)']))<1e-15
    c_observed=ac['i(Vcap)']/(2j*np.pi*f*ac['v(out)']);c_error=float(max(abs(c_observed/params['cload']-1)));assert c_error<tol['capacitor_ac_relative']
    gate_residual=-ac['i(Vin)']+(ac['v(out)']-ac['v(gate)'])/params['rf']+current-bank(ac,'input',1)
    output_residual=bank(ac,'input',0)+bank(ac,'load',0)+ac['v(out)']/params['rload']+(ac['v(out)']-ac['v(gate)'])/params['rf']+current+ac['i(Vcap)']
    ac_kcl=max(ac_kcl,float(max(abs(gate_residual))),float(max(abs(output_residual))))
    assert ac_kcl/plan['conditions']['current_ac_a']<tol['ac_kcl_relative_to_probe']
    forecast=list(csv.DictReader((root/plan['forecast']).open()));oldf=np.array([float(r['frequency_hz']) for r in forecast]);assert np.allclose(f,oldf,rtol=1e-12,atol=0)
    expected=np.array([complex(float(r['branch_transfer_real_ohm']),float(r['branch_transfer_imag_ohm'])) for r in forecast]);measured=ac['v(out)']/current
    error=float(max(abs(measured-expected)/abs(expected)));assert error<tol['branch_complex_relative']
    np.savez_compressed(out/'observations.npz',frequency_hz=f,output_v=ac['v(out)'],gate_v=ac['v(gate)'],input_current_a=-ac['i(Vin)'],branch_current_a=current,capacitor_current_a=ac['i(Vcap)'],normalized_branch_ohm=measured)
    result=dict(schema='learning.noise-branch-measurement.v1',status='PASS',role='Direct finite nominal branch-path diagnostic; later repeats are known-input reproduction',
        new_run_id=report['run_id'],realization_id=real['content_id'],original_realization_id=original['realization_id'],plan_sha256=digest((root/'branch-plan.json').read_bytes()),
        circuit_sha256=plan['circuit_sha256'],forecast_sha256=plan['forecast_sha256'],rows=[1,241],ngspice_version=version,binary_sha256=digest(a.binary.read_bytes()),apm_commit=source['apm_commit'],
        maximum_complex_forecast_relative_error=error,measured_one_khz_ohm=[float(measured[0].real),float(measured[0].imag)],dc_node_differences_v=dc_errors,
        maximum_MOS_dc_kcl_a=dc_kcl,maximum_ac_kcl_a=ac_kcl,maximum_terminal_span_v=span,source_current_readback_relative_error=source_error,capacitor_ac_relative_error=c_error,
        resistor_readback_ohm=observed,resistor_readback_relative_error=relative,physical_readback='Pinned APM model/W/L/m/nf/DELVTO/MULU0 before/applied/after; measured DC V/I for all six resistors and AC current/voltage for capacitor.',
        observations_sha256=digest((out/'observations.npz').read_bytes()),limits='One nominal linear AC transfer diagnostic with ideal current/voltage instrumentation. It is not a new full noise analysis, noise population, physical amplitude limit, or regulator.')
    (out/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
