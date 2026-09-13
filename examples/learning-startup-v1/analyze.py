"""MIT: R75 ramp current, charge and unfitted RC explanations from saved data."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parent
WINDOW_S = 2e-6


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def inputs(source_zip):
    source, index = read(ROOT/'data/source.json'), read(ROOT/'data/index.json')
    assert sha(source_zip) == source['package_sha256'] == index['source_zip_sha256']
    with zipfile.ZipFile(source_zip) as z:
        selected = {n: z.read('example/'+n) for n in source['selected_zip_members']}
    for name, h in source['selected_zip_members'].items():
        assert hashlib.sha256(selected[name]).hexdigest() == h, name
    cfg = json.loads(selected['data/configuration.json'])
    original = {r['case']['id']: r for r in json.loads(selected['data/observations.json'])['results']}
    histories = {r['case_id']: r for r in csv.DictReader(io.StringIO(selected['data/histories.csv'].decode()))}
    plateaus = {r['case_id']: r for r in csv.DictReader(io.StringIO(selected['data/plateaus.csv'].decode()))}
    assert cfg['apm_commit'] == index['apm_commit'] and cfg['temperature_c'] == 26.85
    arrays = {}
    for name, row in index['cases'].items():
        case = cfg['cases'][name]
        assert sha(ROOT/'data'/row['file']) == row['sha256']
        assert (ROOT/'circuits'/(name+'.cir')).read_bytes() == selected['circuits/'+name+'.cir']
        assert sha(ROOT/'circuits'/(name+'.cir')) == case['circuit_sha256']
        assert row['original_attempt_sha256'] == case['attempt_sha256'] == original[name]['attempt_sha256']
        assert row['original_realization_id'] == case['original_realization_id']
        assert len(case['physical_devices']) == 6 and all(u['raw'] == [0., 0.] for u in case['physical_devices'])
        with np.load(ROOT/'data'/row['file'], allow_pickle=False) as z:
            arrays[name] = []
            for j, recipe in enumerate(case['analyses']):
                a = z['analysis'+str(j)]
                assert a.shape == (row['rows'][j], 1+len(recipe['vectors'])) and np.isfinite(a).all()
                arrays[name].append((a[:,0], {n:a[:,i+1] for i,n in enumerate(recipe['vectors'])}))
    return cfg, original, histories, plateaus, arrays


def passive_ramp(t, duration, vmax, v0, alpha, tau):
    x = np.clip(t, 0., duration)
    during = alpha*vmax/duration*(x+tau*np.expm1(-x/tau))+v0*np.exp(-x/tau)
    end = alpha*vmax/duration*(duration+tau*np.expm1(-duration/tau))+v0*np.exp(-duration/tau)
    after = alpha*vmax+(end-alpha*vmax)*np.exp(-np.maximum(t-duration,0.)/tau)
    return np.where(t <= duration, during, after)


def interval(t, values, start, stop):
    axis = np.r_[start, t[(t > start) & (t < stop)], stop]
    assert np.all(np.diff(axis) > 0)
    return axis, {n:np.interp(axis,t,v) for n,v in values.items()}


def csvfile(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); assert not a.output.exists()
    cfg, original, histories, plateaus, arrays = inputs(a.source_zip)
    static_name = 'R75-both_on'; p0 = cfg['cases'][static_name]['parameters']
    op = {n:float(v[0]) for n,v in arrays[static_name][0][1].items()}
    initial = float(arrays['R75-input_only'][0][1]['v(out)'][0])
    final = op['v(out)']; rd, rl, rs = (p0[k] for k in ['rd','rload','rs'])
    cap = 8e-12; conductance = 1/rd+1/rl; tau = cap/conductance; alpha = (1/rd)/conductance
    for element, expected in [('rd',rd),('rload',rl),('rs',rs),('rsrc',p0['rsrc'])]:
        assert np.isclose(op['@'+element+'[resistance]'],expected,rtol=2e-11,atol=1e-12)
    assert np.isclose(op['@cload[capacitance]'],p0['cload'],rtol=2e-11,atol=1e-22)
    sums = {d:sum(op[f'@m.xn{i}.mapm045_vtg_core[{d}]'] for i in range(6)) for d in ['gm','gmbs','gds']}
    source_denominator = 1+rs*sum(sums.values())
    effective_gds = sums['gds']/source_denominator
    tau_on = cap/(conductance+effective_gds)
    # Independently solve output/source KCL for a unit output-current injection.
    total = sum(sums.values())
    node_matrix = np.array([[conductance+sums['gds'], -total],
                            [-sums['gds'], 1/rs+total]])
    node_impedance = float(np.linalg.solve(node_matrix,[1.,0.])[0])
    assert np.isclose(node_impedance,1/(conductance+effective_gds),rtol=1e-12)
    results, table, waves = {}, [], []
    for speed in ['fast','slow']:
        name = 'R75-input_first-'+speed+'-8p'; selected = cfg['cases'][name]
        case, params = selected['case'], selected['parameters']
        duration = {'fast':1e-8,'slow':1e-6}[speed]; start=2e-6; end=start+duration
        assert case['ramp_time_s']==duration and case['sequence']=='input_first'
        assert case['capacitance_f']==cap and case['stop_s']==20e-6
        assert selected['physical_devices']==cfg['cases'][static_name]['physical_devices']
        for key in ['rd','rload','rs','rsrc','nn','nw','nl']:
            assert params[key]==p0[key], key
        t,v=arrays[name][1]; assert t[0]==0 and abs(t[-1]-20e-6)<1e-15 and np.all(np.diff(t)>0)
        assert max(np.diff(t)) <= case['tmax_s']*(1+1e-6)
        for field, vector in [('supply_points','v(vdd)'),('input_points','v(in)')]:
            points=np.array(case[field]); assert max(abs(v[vector]-np.interp(t,points[:,0],points[:,1])))<1e-9
        dc=arrays[name][0][1]
        assert np.isclose(dc['@cload[capacitance]'][0],cap,rtol=2e-11,atol=1e-22)
        assert abs(v['v(out)'][0]-initial)<1e-6
        idd=-v['i(Vdd)']; iin=-v['i(Vin)']; ir=(v['v(vdd)']-v['v(out)'])/rd
        il=v['v(out)']/rl; ic=v['i(Vcap)']; inm=sum(v[f'i(Vn{i}d)'] for i in range(6))
        checks=dict(supply_resistor_kcl_max_a=float(max(abs(idd-ir))),
                    output_kcl_max_a=float(max(abs(ir-il-ic-inm))),
                    MOS_terminal_kcl_max_a=float(max(max(abs(sum(v[f'i(Vn{i}{port})'] for port in ['d','g','s','b']))) for i in range(6))))
        assert max(checks.values())<1e-10
        peak_index=int(np.argmax(idd)); peak=float(idd[peak_index])
        assert peak==original[name]['measurements']['transient']['ports']['Vdd']['max_delivered_current_a']
        assert peak==float(histories[name]['Vdd_max_delivered_current_a'])
        h=plateaus[name]; assert h['state']=='both_on' and abs(float(h['origin_s'])-end)<1e-18
        assert float(h['target_v'])==final
        ramp_out=float(np.interp(end,t,v['v(out)'])); ramp_dd=float(np.interp(end,t,idd))
        ramp_passive=float(passive_ramp(np.array([duration]),duration,1.,initial,alpha,tau)[0])
        tail_error=abs(ramp_out-final)
        assert tail_error>.001
        estimated_tail=tau_on*np.log(tail_error/.001)
        wa,wv=interval(t,dict(output=v['v(out)'],supply=v['v(vdd)'],input=v['v(in)'],idd=idd,iin=iin,ic=ic,il=il,inm=inm),start,start+WINDOW_S)
        charge=float(np.trapz(wv['ic'],wa)); expected=cap*float(wv['output'][-1]-wv['output'][0]); charge_error=charge/expected-1
        assert abs(charge_error)<cfg['plan']['measurements']['charge_relative_tolerance']
        window=dict(start_s=start,stop_s=start+WINDOW_S,duration_s=WINDOW_S,
            capacitor_charge_integral_c=charge,capacitor_charge_C_delta_V_c=expected,capacitor_charge_relative_difference=charge_error,
            supply_signed_charge_c=float(np.trapz(wv['idd'],wa)),input_signed_charge_c=float(np.trapz(wv['iin'],wa)),
            drain_terminal_charge_c=float(np.trapz(wv['inm'],wa)),resistive_load_charge_c=float(np.trapz(wv['il'],wa)),
            supply_signed_energy_j=float(np.trapz(wv['supply']*wv['idd'],wa)),input_signed_energy_j=float(np.trapz(wv['input']*wv['iin'],wa)),
            external_capacitor_energy_change_j=cap/2*float(wv['output'][-1]**2-wv['output'][0]**2))
        window['integrated_output_KCL_residual_c']=window['supply_signed_charge_c']-window['drain_terminal_charge_c']-window['resistive_load_charge_c']-charge
        bracket=[float(h['relative_bracket_start_s']),float(h['relative_bracket_end_s'])]
        result=dict(rows=len(t),ramp_time_s=duration,native_maximum_timestep_s=float(max(np.diff(t))),
            current_peak=dict(supply_delivered_a=peak,relative_to_ramp_start_s=float(t[peak_index]-start),
                supply_voltage_v=float(v['v(vdd)'][peak_index]),output_voltage_v=float(v['v(out)'][peak_index]),
                capacitor_a=float(ic[peak_index]),load_resistor_a=float(il[peak_index]),NMOS_drain_terminal_sum_a=float(inm[peak_index])),
            ramp_end=dict(observed_output_v=ramp_out,observed_supply_current_a=ramp_dd,passive_screen_output_v=ramp_passive,passive_screen_supply_current_a=(1-ramp_passive)/rd),
            original_arrival=dict(status=h['status'],target_v=final,band_v=.001,post_ramp_bracket_s=bracket,
                total_from_ramp_start_bracket_s=[x+duration for x in bracket],horizon_s=float(h['horizon_s'])),
            retrospective_local_tail=dict(ramp_end_absolute_error_v=tail_error,tau_s=tau_on,estimated_post_ramp_time_s=float(estimated_tail),
                estimate_minus_original_upper_s=float(estimated_tail-bracket[1]),role='Unfitted final-OP single-capacitor estimate anchored at already observed ramp-end error; not prospective startup prediction.'),
            checks=checks,common_window=window)
        results[name]=result
        table.append(dict(case=name,ramp_s=duration,supply_peak_a=peak,ramp_end_output_v=ramp_out,
            passive_ramp_end_output_v=ramp_passive,passive_ramp_end_supply_a=(1-ramp_passive)/rd,
            post_ramp_arrival_upper_s=bracket[1],total_arrival_upper_s=duration+bracket[1],local_tail_estimate_s=float(estimated_tail),
            capacitor_charge_c=charge,supply_charge_c=window['supply_signed_charge_c'],supply_energy_j=window['supply_signed_energy_j'],capacitor_energy_change_j=window['external_capacitor_energy_change_j']))
        # One compact drawing grid per history; never used for measurements.
        grid=np.unique(np.r_[0.,np.geomspace(1e-10,WINDOW_S,650),duration,t[peak_index]-start])
        assert grid.min()>=0 and grid.max()<=WINDOW_S
        actual={k:np.interp(start+grid,t,q) for k,q in dict(output_v=v['v(out)'],supply_v=v['v(vdd)'],supply_a=idd,capacitor_a=ic,load_a=il,nmos_drain_a=inm).items()}
        screen=passive_ramp(grid,duration,1.,initial,alpha,tau)
        for j,time in enumerate(grid):
            waves.append(dict(case=name,relative_time_s=float(time),**{k:float(v[j]) for k,v in actual.items()},passive_output_v=float(screen[j]),passive_supply_a=float((actual['supply_v'][j]-screen[j])/rd)))
    result=dict(schema='learning.startup-current.analysis.v1',role='Retrospective four-record calculation; original59-case study and exposure history unchanged.',
        source_zip_sha256=sha(a.source_zip),conditions=dict(candidate='R75',physical_units=6,raw='Saved six zero pairs',temperature_c=26.85,
            Rd_ohm=rd,Rload_ohm=rl,Rs_ohm=rs,Rsrc_ohm=p0['rsrc'],external_capacitance_f=cap,input_v=.45,supply_final_v=1.,ramp_start_s=2e-6,original_horizon_s=20e-6),
        reference=dict(initial_output_v=initial,final_output_v=final,final_supply_delivered_a=-op['i(Vdd)'],final_input_delivered_a=-op['i(Vin)']),
        passive_screen=dict(assumption='Omit all NMOS terminal current and intrinsic charging. Include Rd, Rload and external capacitor; same independently observed initial output.',
            tau_s=tau,alpha=alpha,final_output_v=alpha,short_ramp_current_screen_a=1/rd),
        final_OP_local_tail=dict(derivative_sums_s=sums,source_feedback_denominator=source_denominator,effective_gds_s=effective_gds,tau_s=tau_on,
            independent_output_source_matrix_s=node_matrix.tolist(),independent_output_transimpedance_ohm=node_impedance,
            assumptions='Fixed input gate in the local model; neglect gate/leakage and intrinsic capacitance derivatives. Retain source and body motion through Rs. No fit to transient time or waveform.'),
        cases=results,limits=['Two known input-first ramps, not arbitrary initial charge, all-history reachability, finite upstream-driver behavior or new qualification.',
            'Current maxima and arrival brackets retain the original finite time grids and definitions.',
            'The2us common energy/charge window contains different amounts of nominal operation after each startup; it is not equal-service efficiency or compact-model stored-energy accounting.',
            'MOS drain-terminal current includes compact-model conduction and displacement; it is not a separately measured channel-current component.'])
    a.output.mkdir(parents=True)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    csvfile(a.output/'comparison.csv',table);csvfile(a.output/'waveforms.csv',waves)
    print(json.dumps(dict(cases=len(results),plot_rows=len(waves),results=table),indent=2))


if __name__ == '__main__':
    main()
