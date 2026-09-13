"""MIT: selected learning.bridge.v1 saved-data calculations; no SPICE/imported runner.

Inputs are original float64 arrays at original points for four selected v2
cases, plus previously published ADS nominal values. Existing definitions are
not edited. New DC-endpoint settling is labeled separately from old median-
plateau settling. All output goes to a new explicit directory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs(here):
    manifest = json.loads((here/'data/manifest.json').read_text())
    for name, h in manifest['files'].items():
        assert sha(here/'data'/name) == h, name
    cfg = json.loads((here/'data/configuration.json').read_text())
    arrays = {}
    for key, case in cfg['cases'].items():
        assert sha(here/'circuits'/(case['stem']+'.cir')) == case['circuit_sha256']
        arrays[key] = np.load(here/'data'/case['array_file'], allow_pickle=False)
    return cfg, arrays


def column(cfg, arrays, key, i, name):
    recipe = next(d['recipe'] for d in cfg['cases'][key]['analyses'] if d['analysis_index'] == i)
    j = recipe['vectors'].index(name)
    a = arrays[key]['analysis'+str(i)]
    if recipe['kind'] == 'ac':
        return a[:, 1+2*j] + 1j*a[:, 2+2*j]
    return a[:, 1+j]


def crossings(x, y):
    hits = np.flatnonzero(y[:-1]*y[1:] < 0)
    return [float(x[i]-y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])) for i in hits]


def calculate(here):
    cfg, arrays = inputs(here)
    c = lambda key, i, n: column(cfg, arrays, key, i, n)
    axis = lambda key, i: arrays[key]['analysis'+str(i)][:, 0]
    result = dict(schema='learning.bridge.v1', category='Retrospective selected saved-data reanalysis; no new confirmation', cases={})
    # Physical KCL uses signed current into each of four MOS terminals.
    for key, case in cfg['cases'].items():
        kcl = 0.; span = 0.
        for desc in case['analyses']:
            i = desc['analysis_index']
            if desc['kind'] not in ('op', 'dc', 'tran'):
                continue
            for device in case['connections'].values():
                kcl = max(kcl, float(np.max(np.abs(sum(c(key, i, 'i('+s+')') for s in device['sensors'].values())))))
                voltages = [np.zeros(len(axis(key, i))) if n == '0' else c(key, i, 'v('+n+')') for n in device['nodes'].values()]
                span = max(span, float(np.max(np.max(voltages, axis=0)-np.min(voltages, axis=0))))
        assert kcl < 1e-10 and span <= 1+1e-8, (key, kcl, span)
        result['cases'][key] = dict(stem=case['stem'], run_id=case['run_id'], physical_units=len(case['devices']),
            maximum_MOS_KCL_residual_A=kcl, maximum_terminal_span_V=span,
            operating_nodes_V={n:float(c(key, 0, 'v('+n+')')[0]) for n in case['resources']['nodes'] if n != '0'})
    # A: ideal tail, differential current and two resistor loads.
    x = axis('A', 1); left = c('A', 1, 'i(VtM1d)'); right = c('A', 1, 'i(VtM2d)')
    j = int(np.argmin(abs(x))); slope = float((left[j+1]-right[j+1]-left[j-1]+right[j-1])/(x[j+1]-x[j-1]))
    gm = float(c('A', 0, '@m.xm1.mapm045_vtg_core[gm]')[0]); gds = float(c('A', 0, '@m.xm1.mapm045_vtg_core[gds]')[0])
    result['cases']['A'].update(left_current_A=float(left[j]), right_current_A=float(right[j]),
        differential_current_slope_S=slope, differential_gain_1kHz=float(abs(c('A', 2, 'v(odiff)')[0])),
        prediction_gm_times_parallel_load=gm/(1/90000+gds),
        drain_current_sum_range_A=[float(min(left+right)), float(max(left+right))])
    # B: same nominal pair current, added finite tail/reference, common-mode dependence.
    tail = c('B', 1, 'i(VtMtd)')
    result['cases']['B'].update(VDD_delivered_current_A=float(-c('B', 0, 'i(Vdd)')[0]),
        tail_drain_current_A=float(c('B', 0, 'i(VtMtd)')[0]),
        reference_network_current_A=float(-c('B', 0, 'i(Vdd)')[0]-c('B', 0, 'i(VtM1d)')[0]-c('B', 0, 'i(VtM2d)')[0]),
        common_mode_sweep_V=[float(axis('B', 1)[0]),float(axis('B', 1)[-1])],
        tail_current_sweep_range_A=[float(min(tail)),float(max(tail))])
    # C: positive available current enters output from the OTA, before Rload.
    x = axis('C', 1); available = -c('C', 1, 'i(VtM2d)')-c('C', 1, 'i(VtM4d)')
    loaded = available-c('C', 1, 'v(out)')/1e6
    j = int(np.argmin(abs(x))); gm = float((available[j+1]-available[j-1])/(x[j+1]-x[j-1]))
    xo = axis('C', 2); io = -c('C', 2, 'i(VtM2d)')-c('C', 2, 'i(VtM4d)')-c('C', 2, 'v(out)')/1e6
    jo = int(np.argmin(abs(xo-.3))); go = float(-(io[jo+1]-io[jo-1])/(xo[jo+1]-xo[jo-1]))
    result['cases']['C'].update(unloaded_nulls_V=crossings(x, available), loaded_nulls_V=crossings(x, loaded),
        output_load_A=.3/1e6, available_current_at_zero_diff_A=float(available[j]),
        source_sink_current_range_A=[float(min(available)), float(max(available))],
        differential_sweep_V=[float(x[0]),float(x[-1])], gm_S=gm, loaded_gout_S=go,
        predicted_loaded_gain=gm/go, predicted_load_shift_V=(.3/1e6)/gm)
    # D1: functional buffer, original 5-pF inventory. DC excludes two guard points.
    vin = axis('D1', 1); out = c('D1', 1, 'v(out)'); gain = (out[2:]-out[:-2])/(vin[2:]-vin[:-2])
    f = axis('D1', 3); h = c('D1', 3, 'v(out)')/c('D1', 3, 'v(source)')
    cross = crossings(np.log(f), 20*np.log10(abs(h)/abs(h[0]))+3)
    d = result['cases']['D1']
    d.update(input_range_V=[float(vin[1]),float(vin[-2])], output_range_V=[float(out[1]),float(out[-2])],
        worst_absolute_error_V=float(max(abs(out[1:-1]-vin[1:-1]))), incremental_gain_range=[float(min(gain)),float(max(gain))],
        maximum_VDD_delivered_quiescent_A=float(max(-c('D1', 1, 'i(Vdd)')[1:-1])),
        signal_source_delivered_current_range_A=[float(min(-c('D1', 1, 'i(Vin)')[1:-1])),float(max(-c('D1', 1, 'i(Vin)')[1:-1]))],
        signal_gain_1kHz=float(abs(h[0])), bandwidth_minus_3dB_Hz=float(np.exp(cross[0])))
    # Old v2 referral divided total integrated output noise by gain at 1 kHz.
    fnoise = axis('D1', 5); variance = float(np.trapz(c('D1', 5, 'onoise_spectrum'), fnoise))
    d['old_v2_noise'] = dict(band_Hz=[1000,1000000], output_rms_V=float(np.sqrt(variance)),
        input_equivalent_rms_V=float(np.sqrt(variance)/abs(h[0])),
        definition='PSD integrated in linear frequency, then divided by |H(1kHz)|; differs from ADS frequency-by-frequency input PSD referral.')
    t = axis('D1', 12); y = c('D1', 12, 'v(out)')
    lo = float(np.interp(.25, vin, out)); hi = float(np.interp(.35, vin, out)); band = .01*abs(hi-lo)
    d['dc_endpoint_steps'] = dict(definition='New learning.bridge.v1: independent saved DC endpoints, 1% of actual DC step, final observed entry from input-source 50% crossing; original median-plateau metrics not changed.', band_V=band, directions={})
    for name, origin, end, target in [('rise',.50625e-6,4.5e-6,hi),('fall',4.50625e-6,8.5e-6,lo)]:
        use=(t>=origin)&(t<end);tt=t[use];yy=y[use];outside=np.flatnonzero(abs(yy-target)>band)
        assert len(outside) and outside[-1]+1<len(tt), name
        j=int(outside[-1]);d['dc_endpoint_steps']['directions'][name]=dict(status='OBSERVED_FINAL_ENTRY',target_V=target,
            bracket_s=[float(tt[j]-origin),float(tt[j+1]-origin)],time_s=float(tt[j+1]-origin),observed_until_s=float(tt[-1]-origin))
    current=-c('D1',12,'i(VtM2d)')-c('D1',12,'i(VtM4d)')-y/1e6+c('D1',12,'i(Vfeedback)')
    charge=np.r_[0,np.cumsum((current[:-1]+current[1:])*.5*np.diff(t))]
    residual=float(max(abs(charge-5e-12*(y-y[0])))/(5e-12*np.ptp(y)))
    assert residual<.005
    d['output_charge_relative_residual']=residual
    # Independent low-frequency equation/KCL route for the resistor chapter.
    resistor=json.loads((here/'data/resistor-op.json').read_text()); rr={}
    for name, record in resistor.items():
        o=record['metrics']['op'];p=record['parameters']
        gm=sum(v for k,v in o.items() if k.endswith('[gm]'))
        gds=sum(v for k,v in o.items() if k.endswith('[gds]'))
        gmb=sum(v for k,v in o.items() if k.endswith('[gmbs]'))
        drain=sum(v for k,v in o.items() if k.startswith('i(Vn') and k.endswith('d)'))
        kcl=(p['vdd']-o['v(out)'])/p['rd']-o['v(out)']/p['rload']-drain
        D=1+p['rs']*(gm+gmb+gds);go=1/p['rd']+1/p['rload']+gds/D
        prediction=-gm/D/go;observed=record['metrics']['ac']['real_v_per_v']
        assert abs(kcl)<1e-11 and abs(prediction/observed-1)<.001
        rr[name]=dict(gm_S=gm,gds_S=gds,gmb_S=gmb,Id_A=drain,Vgs_V=o['v(gate)']-o['v(source)'],
            gm_over_Id_per_V=gm/drain,D=D,predicted_gain=prediction,measured_gain=observed,
            gain_relative_error=prediction/observed-1,output_KCL_residual_A=kcl,
            predicted_current_error_path_ohm=1/go)
    result['resistor_equations']=rr
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists(), 'Use a new analysis directory'
    result=calculate(Path(__file__).resolve().parent)
    a.output.mkdir(parents=True)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(status='PASS',output=str(a.output),cases=4,resistor_equations=2)))


if __name__=='__main__':main()
