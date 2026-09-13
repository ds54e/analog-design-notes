"""Two known circuits: coherent branch-noise paths and an explicit partial budget."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

ROOT=Path(__file__).resolve().parent
VERSION='learning.resistor-noise-path.v1'
BOLTZMANN_J_PER_K=1.380649e-23
BAND_HZ=(1e3,1e6)


def digest(data):return hashlib.sha256(data).hexdigest()


def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))


def integrate(f,psd):
    """Original ADS definition: endpoint insertion and linear-frequency trapezoids."""
    assert len(f)==len(psd) and np.isfinite(f).all() and np.isfinite(psd).all()
    assert np.all(np.diff(f)>0) and np.all(psd>=0)
    low,high=BAND_HZ;assert f[0]<=low*(1+1e-12) and f[-1]>=high*(1-1e-12)
    axis=np.r_[low,f[(f>low)&(f<high)],high]
    return float(np.trapz(np.interp(axis,f,psd),axis))


def low_frequency_network(case):
    """Independent gate/source/output/bias node equations, ignoring C/leakage derivatives."""
    op=case['observed_op'];p=case['parameters'];u=case['physical_units']
    def bank(role,key):
        return sum(op['@m.'+x['path'][0]+'.mapm045_vtg_core['+key+']'] for x in u if x['uid'].split('/')[-2]==role)
    gm=bank('input','gm');gmb=bank('input','gmbs');gds=bank('input','gds');gp=bank('load','gds');gmp=bank('load','gm')
    gr=bank('reference','gm')+bank('reference','gds');rt=p['rsrc']+p.get('rin',0);gf=1/p['rf'] if p.get('rf') else 0
    rd=1/p['rd'] if p.get('rd') else 0;gt=1/rt;gs=1/p['rs'];gl=1/p['rload'];total=gm+gmb+gds
    n=4 if gr else 3
    # Unknowns: gate, source, output, optional finite reference bias.
    y=np.zeros((n,n));y[0,:3]=[gt+gf,0,-gf];y[1,:3]=[-gm,gs+total,-gds];y[2,:3]=[gm-gf,-total,gds+gp+gl+rd+gf]
    if n==4:y[2,3]=gmp;y[3,3]=gr+1/p['rbias']
    drives=np.zeros((n,3));drives[0,0]=gt;drives[2,1]=1;drives[0,2]=1;drives[2,2]=-1
    h,z,j=np.linalg.solve(y,drives)[2,:]
    denominator=1+p['rs']*total;geff=gm/denominator;go=gds/denominator+gp+gl+rd
    if gf:
        beta=rt/(p['rf']+rt);z_reduced=1/(go+1/(p['rf']+rt)+geff*beta)
        h_reduced=-(geff*p['rf']-1)/(p['rf']+rt)*z_reduced
    else:z_reduced=1/go;h_reduced=-geff/go
    assert abs(h/h_reduced-1)<1e-13 and abs(z/z_reduced-1)<1e-13
    assert abs(j/(rt*h-z)-1)<1e-13
    return dict(gm_s=gm,gmb_s=gmb,gds_s=gds,pmos_gds_s=gp,source_feedback_denominator=denominator,
                effective_gm_s=geff,output_conductance_s=go,signal_gain=float(h),zout_ohm=float(z),
                gate_to_output_branch_ohm=float(j),matrix_vs_reduced_relative=max(abs(h/h_reduced-1),abs(z/z_reduced-1)),
                direct_branch_vs_two_paths_relative=abs(j/(rt*h-z)-1))


def calculate(feedback_zip,data_root=ROOT/'data'):
    source=json.loads((data_root/'source.json').read_text());body=Path(feedback_zip).read_bytes();assert digest(body)==source['old_zip_sha256']
    assert digest((data_root/'current-transfer.csv').read_bytes())==source['current_transfer_sha256']
    if (data_root/'manifest.json').exists():
        for n,h in json.loads((data_root/'manifest.json').read_text())['files'].items():assert digest((data_root/n).read_bytes())==h,n
    for case in source['cases']:
        if case['role']=='current':assert digest((data_root.parent/case['circuit']).read_bytes())==case['circuit_sha256']
    with zipfile.ZipFile(io.BytesIO(body)) as z:members={n:z.read('example/'+n) for n in source['selected_members_sha256']}
    for n,h in source['selected_members_sha256'].items():assert digest(members[n])==h,n
    ac=rows(members['data/ac.csv']);noise=rows(members['data/noise.csv']);nominal=rows(members['data/nominal.csv']);current=rows((data_root/'current-transfer.csv').read_bytes())
    spectrum=[];budget=[];summary=[]
    for name in ['R75','F75L']:
        case=next(c for c in source['cases'] if c['case']==name+'-nominal');p=case['parameters'];temperature=case['temperature_c']+273.15;assert temperature==300
        rt=p['rsrc']+p.get('rin',0);kt4=4*BOLTZMANN_J_PER_K*temperature
        aa=[r for r in ac if r['candidate']==name];nn=[r for r in noise if r['candidate']==name];zz=[r for r in current if r['candidate']==name]
        fa=np.array([float(r['frequency_hz']) for r in aa]);fz=np.array([float(r['frequency_hz']) for r in zz]);f=np.array([float(r['frequency_hz']) for r in nn]);assert np.array_equal(fa,fz)
        assert len(fa)==401 and len(f)==481 and np.all(np.diff(fa)>0)
        h=np.array([complex(float(r['gain_real']),float(r['gain_imag'])) for r in aa]);z=np.array([complex(float(r['zout_real_ohm']),float(r['zout_imag_ohm'])) for r in zz]);assert np.isfinite(h).all() and np.isfinite(z).all()
        at_noise=lambda values:np.interp(np.log(f),np.log(fa),values)
        gain2=at_noise(abs(h)**2);assert min(gain2)>1e-24
        known_output=np.array([float(r['output_psd_v2_per_hz']) for r in nn]);known_input=np.array([float(r['input_psd_v2_per_hz']) for r in nn]);assert np.array_equal(known_output/gain2,known_input)
        # Interpolate output PSDs on the original noise grid, then use the same
        # original |H|^2 interpolation for input referral. Do not square an ASD twice.
        component={'series_input':np.full(len(f),kt4*rt),
                   'output_load':kt4/p['rload']*at_noise(abs(z)**2)/gain2,
                   'drain_resistor':kt4/p['rd']*at_noise(abs(z)**2)/gain2 if p.get('rd') else np.zeros(len(f)),
                   'feedback_resistor':kt4/p['rf']*at_noise(abs(rt*h-z)**2)/gain2 if p.get('rf') else np.zeros(len(f))}
        selected=sum(component.values());residual=known_input-selected
        assert np.all(residual>=0),(name,float(min(residual)))
        wrong=kt4/p['rf']*at_noise(abs(rt*h)**2+abs(z)**2)/gain2 if p.get('rf') else np.zeros(len(f))
        cross=kt4/p['rf']*at_noise(-2*np.real(rt*h*np.conj(z)))/gain2 if p.get('rf') else np.zeros(len(f))
        if p.get('rf'):assert max(abs((wrong+cross)/component['feedback_resistor']-1))<2e-14
        total_var=integrate(f,known_input);public_rms=float(next(r['input_noise_v'] for r in nominal if r['candidate']==name));assert abs(np.sqrt(total_var)/public_rms-1)<2e-14
        for role,values in {**component,'selected_passives':selected,'unassigned_remainder':residual,'observed_total':known_input}.items():
            variance=integrate(f,values);budget.append(dict(candidate=name,component=role,variance_v2=variance,rms_v=float(np.sqrt(variance)),fraction_of_observed_variance=variance/total_var))
        analytic_input_variance=kt4*rt*(BAND_HZ[1]-BAND_HZ[0]);assert abs(integrate(f,component['series_input'])/analytic_input_variance-1)<1e-14
        matrix=low_frequency_network(case);matrix['matrix_signal_vs_native_1khz_relative']=matrix['signal_gain']/h[0].real-1;matrix['matrix_zout_vs_native_1khz_relative']=matrix['zout_ohm']/z[0].real-1
        assert abs(matrix['matrix_signal_vs_native_1khz_relative'])<.0005 and abs(matrix['matrix_zout_vs_native_1khz_relative'])<.0005
        for i,frequency in enumerate(f):
            spectrum.append(dict(candidate=name,frequency_hz=float(frequency),observed_input_psd_v2_per_hz=float(known_input[i]),
                **{key+'_input_psd_v2_per_hz':float(value[i]) for key,value in component.items()},
                selected_passives_input_psd_v2_per_hz=float(selected[i]),unassigned_input_psd_v2_per_hz=float(residual[i]),
                incoherent_feedback_input_psd_v2_per_hz=float(wrong[i]),feedback_cross_input_psd_v2_per_hz=float(cross[i])))
        row=dict(candidate=name,temperature_k=temperature,input_series_ohm=rt,band_hz=list(BAND_HZ),original_input_noise_rms_v=public_rms,
            partial_passive_rms_v=float(np.sqrt(integrate(f,selected))),partial_fraction_of_total_variance=integrate(f,selected)/total_var,
            minimum_unassigned_psd_v2_per_hz=float(min(residual)),matrix=matrix,
            one_khz=dict(signal_real=float(h[0].real),signal_imag=float(h[0].imag),zout_real_ohm=float(z[0].real),zout_imag_ohm=float(z[0].imag),
                        gate_path_real_ohm=float((rt*h[0]).real),opposite_output_path_real_ohm=float(-z[0].real),
                        combined_branch_real_ohm=float((rt*h[0]-z[0]).real),combined_branch_imag_ohm=float((rt*h[0]-z[0]).imag)))
        if p.get('rf'):
            correct=integrate(f,component['feedback_resistor']);uncorrelated=integrate(f,wrong)
            row['feedback_branch']=dict(correct_rms_v=float(np.sqrt(correct)),incorrect_independent_ends_rms_v=float(np.sqrt(uncorrelated)),
                 cross_fraction_of_correct_variance=(correct-uncorrelated)/correct,
                 correct_over_incorrect_variance=correct/uncorrelated,
                 one_khz_equivalent_input_resistance_ohm=float(abs(rt-z[0]/h[0])**2/p['rf']),
                 real_path_signs_add_at_1khz=bool((rt*h[0]).real<0 and -z[0].real<0))
        summary.append(row)
    result=dict(schema=VERSION,status='PASS',role='Retrospective resistor-budget calculation from known nominal observations; direct branch diagnostic separately identified when supplied',
       boltzmann_j_per_k=BOLTZMANN_J_PER_K,thermal_model='4*k*T/R one-sided current PSD for an ordinary resistor at actual 300 K. SI Boltzmann constant used for this analytic contribution.',
       source_sha256=digest((data_root/'source.json').read_bytes()),cases=summary,
       definitions={'referral':'Each output PSD divided by the original log-frequency interpolation of |H|^2 onto the noise grid.',
                    'branch_transfer':'Current from output to gate gives Rt*H-Zout; square the complete complex sum before adding independent-source powers.',
                    'integration':'Original linear-frequency trapezoid with exact 1 kHz/1 MHz endpoints inserted.',
                    'remainder':'Observed total PSD minus selected ordinary-resistor PSDs; unassigned MOS/other-passive sources, not a complete BSIM decomposition.'},
       limits=['Original signal/noise records are known nominal evidence; the separately predeclared direct branch diagnostic is identified below. No calibrated noise population or foundry qualification.',
               'Complex transfers are sampled at 80 points/decade, noise at 160; output powers are interpolated in log frequency before the original referral. This is a finite sampled approximation.',
               'Low-frequency nodal check neglects capacitance and leakage derivatives; it is compared at 1 kHz, not claimed valid over all frequencies.',
               'The two terminals of one resistor share the same noise source. Independent resistors contribute independent source powers; no independent-terminal assumption.',
               'No noise is added by the ideal external generators; a real driver/reference can add more.',
               'Same-function R75/F75L designs have different device/current/resistor allocations; this analysis does not isolate feedback as the only redesign change.'])
    if (data_root/'branch-result.json').exists():
        # This aggregation was added after the separately frozen direct branch
        # target. It does not alter that script, forecast, tolerance or native result.
        branch=json.loads((data_root/'branch-result.json').read_text())
        plan_path=data_root.parent/'branch-plan.json';plan=json.loads(plan_path.read_text())
        assert branch['status']=='PASS' and branch['plan_sha256']==digest(plan_path.read_bytes())
        assert digest((data_root/'branch-observations.npz').read_bytes())==branch['observations_sha256']
        assert digest((data_root/'prebranch-summary.json').read_bytes())==plan['known_analysis_summary_sha256']
        w=np.load(data_root/'branch-observations.npz',allow_pickle=False)
        measured=w['output_v']/w['branch_current_a']
        assert np.array_equal(measured,w['normalized_branch_ohm'])
        forecast=rows((data_root.parent/plan['forecast']).read_bytes())
        expected=np.array([complex(float(r['branch_transfer_real_ohm']),float(r['branch_transfer_imag_ohm'])) for r in forecast])
        error=float(max(abs(measured-expected)/abs(expected)))
        assert error==branch['maximum_complex_forecast_relative_error'] and error<plan['tolerances']['branch_complex_relative']
        hrows=[r for r in ac if r['candidate']=='F75L'];nrows=[r for r in noise if r['candidate']=='F75L']
        fa=np.array([float(r['frequency_hz']) for r in hrows]);h=np.array([complex(float(r['gain_real']),float(r['gain_imag'])) for r in hrows])
        fn=np.array([float(r['frequency_hz']) for r in nrows]);g2=np.interp(np.log(fn),np.log(fa),abs(h)**2)
        actual_psd=(4*BOLTZMANN_J_PER_K*300/1e5)*np.interp(np.log(fn),np.log(w['frequency_hz']),abs(measured)**2)/g2
        actual_variance=integrate(fn,actual_psd);prior=next(r['variance_v2'] for r in budget if r['candidate']=='F75L' and r['component']=='feedback_resistor')
        variance_change=actual_variance/prior-1
        assert abs(variance_change)<3e-6
        result['direct_branch_diagnostic']=dict(role='One predeclared direct-source AC diagnostic, now exposed; aggregation of its saved observations',
            rows=branch['rows'],plan_sha256=branch['plan_sha256'],new_run_id=branch['new_run_id'],
            maximum_complex_forecast_relative_error=error,feedback_rms_from_direct_native_path_v=float(np.sqrt(actual_variance)),
            feedback_variance_relative_to_prior_two_path_result=variance_change,
            physical_scope='Same five nominal MOS units and bias/load in a new ideal branch-source/sensor fixture; no new full noise run or physical sample.')
    return result,{'noise-budget.csv':budget,'noise-spectrum.csv':spectrum}


def main():
    p=argparse.ArgumentParser();p.add_argument('--feedback-zip',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert a.output.is_absolute() and not a.output.exists();result,tables=calculate(a.feedback_zip);a.output.mkdir(parents=True)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    for name,table in tables.items():
        with (a.output/name).open('x',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(table[0]),lineterminator='\n');writer.writeheader();writer.writerows(table)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
