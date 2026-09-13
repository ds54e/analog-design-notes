"""MIT: replay known F75L static curvature and its original 20-mV sine.

The public APM API uses the selected five zero-raw units and exact circuit.
This is reproduction of exposed conditions, not a new scientific cohort.
"""
import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from analyze import calculate, sha

def periodic(t, y):
    """Original eight-cycle, ten-harmonic measurement; no fitted tail."""
    import numpy as np
    assert len(t) == len(y) and np.isfinite(t).all() and np.isfinite(y).all()
    assert np.all(np.diff(t) > 0)
    stop = float(t[-1]); start = stop-8/1e4
    assert start >= t[0] and max(np.diff(t[t >= start])) <= 1.1e-7
    grid = start+np.arange(8000)/1e7; values = np.interp(grid, t, y)
    ft = np.fft.rfft(values)/len(values)
    amps = np.array([2*abs(ft[n*8]) for n in range(1, 11)])
    assert amps[0] >= 1e-9
    cycles = values.reshape(8, 1000)
    delta = float(np.max(abs(cycles[-1]-cycles[-2])))
    return dict(fundamental_peak_v=float(amps[0]), thd_fraction=float(np.linalg.norm(amps[1:])/amps[0]),
                stationary=bool(delta < max(1e-7, amps[0]*1e-4)), max_last_cycle_difference_v=delta)


def main():
    p=argparse.ArgumentParser()
    for name in ['apm','feedback-zip','output']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--binary',type=Path,default=Path('/usr/local/bin/ngspice'))
    a=p.parse_args();here=Path(__file__).resolve().parent;apm=a.apm.resolve();out=a.output
    assert out.is_absolute() and not out.exists() and not out.resolve().is_relative_to(apm)
    expected,_,_=calculate(a.feedback_zip)
    source=json.loads((here/'data/source.json').read_text())
    static=json.loads((here/'data/static-source.json').read_text())
    plan=json.loads((here/'replay-plan.json').read_text());tol=plan['tolerances']
    assert hashlib.sha256(inspect.getsource(periodic).rstrip('\n').encode()).hexdigest()==plan['measurement_function']['verbatim_body_sha256']
    case=next(r for r in source['cases'] if r['case']=='F75L-20m');p0=case['parameters']
    version=subprocess.check_output([str(a.binary),'--version'],text=True).strip()
    assert 'ngspice-47' in version
    assert subprocess.check_output(['git','-C',str(apm),'rev-parse','HEAD'],text=True).strip()==source['apm_commit']
    for path,h in source['model_hashes'].items():assert sha(apm/path)==h,path
    body=here/case['circuit'];assert sha(body)==case['circuit_sha256']
    lines=[line.split() for line in body.read_text().splitlines() if line.split() and not line.startswith('*')]
    sensors={f[2]:(f[0],f[1]) for f in lines if f[0].startswith('V') and len(f)==4 and f[3]=='0'}
    connections={}
    for f in lines:
        if f[0].startswith('x'):
            assert all(n in sensors for n in f[1:5])
            connections[f[0]]={'sensors':[sensors[n][0] for n in f[1:5]],'nodes':[sensors[n][1] for n in f[1:5]]}
    assert set(connections)=={u['path'][0] for u in case['physical_units']}
    vectors=case['transient_recipe']['vectors'][:]
    for v in ['v(drive)','v(pbias)','v(pdrive)']:
        if v not in vectors:vectors.append(v)
    assert all(not v.startswith('@') for v in vectors)
    recipes=[dict(kind='op',set_sources={'Vin':value},vectors=vectors) for value in plan['static_input_v']]
    recipes.append(dict(kind='ac',start=1e3,stop=1e8,points_per_decade=80,set_sources={'Vin':.45},vectors=vectors))
    tran=case['transient_recipe'].copy();tran['vectors']=vectors;tran['set_sources']={'Vin':.45}
    assert tran['step']==tran['tmax']==plan['transient_tmax_s'] and tran['stop']==plan['transient_stop_s']
    recipes.append(tran)
    os.environ['APM_REPO_ROOT']=str(apm);sys.path.insert(0,str(apm/'src'))
    from apm import research,research_spice as spice
    devices=[{k:v for k,v in d.items() if k!='raw'} for d in case['physical_units']]
    request=dict(schema=research.SCHEMAS['request'],circuit=str(body),devices=devices,other_variation_leaves=[],analyses=recipes)
    binding={};spice.flatten(body,binding)
    for model in spice.MODELS:spice.flatten(apm/model,binding)
    real=research.seal(dict(schema=research.SCHEMAS['realization'],origin='nominal',profile_tier='NOMINAL',status='RESOLVED',
        sample_context_id=research.sample_context(request),request_id=research.canonical_hash(request),input_binding=binding,
        devices=case['physical_units'],public_replay=dict(original_realization_id=case['rebound_realization_id'],
        rule='Exact saved F75L20m circuit/five zero-raw units; relocated paths, terminal-vector recipes and three known DC source values. No sampling, trim or redesign.')))
    out.mkdir(parents=True);(out/'relocated-realization.json').write_text(json.dumps(real,indent=2)+'\n')
    report=spice.execute(apm,a.binary,out/'native',request,real,temperature_c=case['temperature_c'])
    assert report['status']=='PASS',report
    folder=Path(report['directory']);data=[];axes=[];kcl=0.;span=0.;rows=[]
    for i,r in enumerate(recipes):
        z=np.loadtxt(folder/f'analysis{i}.txt',skiprows=1,ndmin=2)
        assert np.isfinite(z).all() and z.shape[1]==1+len(vectors)*(2 if r['kind']=='ac' else 1)
        v={key:(z[:,1+2*j]+1j*z[:,2+2*j] if r['kind']=='ac' else z[:,1+j]) for j,key in enumerate(vectors)}
        t=z[:,0];axes.append(t);data.append(v);rows.append(len(t))
        if r['kind']=='op':assert len(t)==1
        else:assert np.all(np.diff(t)>0) and abs(t[-1]-r['stop'])<max(1e-12,r['stop']*1e-10)
        if r['kind']!='ac':
            for conn in connections.values():
                currents=[v['i('+s+')'] for s in conn['sensors']]
                volts=[np.zeros(len(t)) if n=='0' else v['v('+n+')'] for n in conn['nodes']]
                kcl=max(kcl,float(np.max(abs(sum(currents)))))
                span=max(span,float(np.max(np.ptp(np.stack(volts),axis=0))))
    assert kcl<tol['mos_kcl_a'] and span<=1+tol['terminal_span_slack_v']
    old_static={r['input_delta_v']:r['original_op']['v(out)'] for r in static['cases'] if r['candidate']=='F75L'}
    op_differences=[]
    for i,x in enumerate(plan['static_input_v']):
        delta=round(x-.45,8);v=data[i]
        assert abs(v['v(in)'][0]-x)<1e-12
        difference=float(v['v(out)'][0]-old_static[delta]);assert abs(difference)<tol['op_voltage_v']
        op_differences.append(difference)
    dc_slope=float((data[2]['v(out)'][0]-data[0]['v(out)'][0])/.0002)
    dc_curve=float((data[2]['v(out)'][0]-2*data[1]['v(out)'][0]+data[0]['v(out)'][0])/1e-8)
    old_curve=next(r for r in expected['static_finite_differences'] if r['candidate']=='F75L')
    assert abs(dc_slope/old_curve['slope_v_per_v']-1)<tol['dc_slope_relative']
    assert abs(dc_curve/old_curve['second_derivative_per_v']-1)<tol['dc_curvature_relative']
    # Resistor values come from native voltage/current observations, including
    # feedback current obtained from gate KCL without inserting Rf into it.
    v=data[1];nmos=[n for n in connections if n.startswith('xn')]
    current=lambda name,j:float(v['i('+connections[name]['sensors'][j]+')'][0])
    gate=sum(current(n,1) for n in nmos);incoming=-float(v['i(Vin)'][0]);feedback=gate-incoming
    source_current=-sum(current(n,2) for n in nmos)
    drain_current=sum(current(n,0) for n,c in connections.items() if c['nodes'][0]=='out')
    bias_current=-current('xr0',0)-current('xr0',1)-float(v['i(Vbiasdiag)'][0])
    readbacks={}
    flows={'Rsrc':(float(v['v(in)'][0]-v['v(drive)'][0]),incoming,p0['rsrc']),
        'Rin':(float(v['v(drive)'][0]-v['v(gate)'][0]),incoming,p0['rin']),
        'Rf':(float(v['v(out)'][0]-v['v(gate)'][0]),feedback,p0['rf']),
        'Rs':(float(v['v(source)'][0]),source_current,p0['rs']),
        'Rload':(float(v['v(out)'][0]),-drain_current-feedback,p0['rload']),
        'Rbias':(float(v['v(pbias)'][0]),bias_current,p0['rbias'])}
    for name,(drop,flow,declared) in flows.items():
        assert abs(flow)>1e-12;actual=drop/flow;error=actual/declared-1
        assert abs(error)<tol['resistor_relative'],(name,error)
        readbacks[name]=dict(observed_v_over_i_ohm=actual,declared_ohm=declared,relative_difference=error)
    f=axes[3];v=data[3];cap=v['i(Vcap)']/(2j*np.pi*f*v['v(capnode)'])
    mask=(f>=1e4)&(f<=1e7);cap_error=float(np.max(abs(cap[mask]/p0['cload']-1)))
    assert cap_error<tol['capacitor_ac_relative']
    t=axes[4];v=data[4];assert t[0]==0 and max(np.diff(t))<=tran['tmax']*(1+1e-6)
    with np.load(here/'data/port-cycles.npz',allow_pickle=False) as z:saved=z['F75L-20m']
    grid=saved[:,0];differences={}
    for field,vector,sign in [('input_v','v(in)',1),('output_v','v(out)',1),('gate_v','v(gate)',1),('source_v','v(source)',1),
            ('supply_delivered_current_a','i(Vdd)',-1),('input_delivered_current_a','i(Vin)',-1),('capacitor_current_a','i(Vcap)',1)]:
        observed=np.interp(grid,t,sign*v[vector]);reference=saved[:,source['columns'].index(field)]
        difference=float(np.max(abs(observed-reference)));limit=tol['cycle_voltage_v'] if field.endswith('_v') else tol['cycle_current_a']
        assert difference<limit,(field,difference);differences[field]=difference
    start=np.flatnonzero(t>=.0009)[0];tt=t[start:];ic=v['i(Vcap)'][start:];vo=v['v(out)'][start:]
    charge=np.r_[0.,np.cumsum(np.diff(tt)*(ic[1:]+ic[:-1])/2)]
    charge_error=float(np.max(abs(charge-p0['cload']*(vo-vo[0])))/(p0['cload']*np.ptp(vo)))
    assert charge_error<tol['capacitor_charge_relative']
    eight=periodic(t,v['v(out)']);old_thd=float(case['original_metrics']['thd_fraction'])
    assert eight['stationary'] and abs(eight['thd_fraction']/old_thd-1)<tol['thd_relative']
    assert eight['thd_fraction']<.02
    result=dict(status='PASS',role='Replay of known F75L conditions, not new scientific targets',
        apm_commit=source['apm_commit'],ngspice_version=version,binary_sha256=sha(a.binary),new_run_id=report['run_id'],
        original_sine_realization_id=case['rebound_realization_id'],rows=rows,op_output_differences_v=op_differences,
        dc_slope_v_per_v=dc_slope,dc_second_derivative_per_v=dc_curve,
        dc_curvature_relative_difference=dc_curve/old_curve['second_derivative_per_v']-1,
        cycle_maximum_differences=differences,eight_cycle_measurement=eight,
        old_eight_cycle_thd_fraction=old_thd,resistor_readbacks=readbacks,
        capacitor_ac_relative_error=cap_error,capacitor_charge_relative_residual=charge_error,
        maximum_MOS_terminal_kcl_a=kcl,maximum_terminal_span_v=span,
        physical_readback='Pinned APM checks actual model/W/L/m/nf/DELVTO/MULU0 before/applied/after; V/I, AC current and charge check actual passives.',
        limits='One known five-unit nominal design. Same-host public recipe reproduction does not qualify a new population or high-amplitude range.')
    (out/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
