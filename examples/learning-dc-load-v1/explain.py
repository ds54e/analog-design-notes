"""MIT: retrospective OP network and DC power balance after one load comparison."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from prior import read_source, explicit_matrix, stamp, derivatives

ROOT=Path(__file__).resolve().parent


def power_balance(op,connections,rload,rbias):
    def voltage(node):return 0. if node=='0' else op['v('+node+')']
    mos={name:sum(voltage(c['nodes'][terminal])*op['i('+sensor+')']
                 for terminal,sensor in c['sensors'].items()) for name,c in connections.items()}
    input_current=-op['i(Vin)'];supply_current=-op['i(Vdd)']
    resistor={'Rload':op['v(out)']**2/rload,'Rbias':op['v(bias)']**2/rbias,
              'Rsource':(op['v(source)']-op['v(inp)'])*input_current}
    delivered={'Vdd':op['v(vdd)']*supply_current,'Vin':op['v(source)']*input_current}
    residual=sum(delivered.values())-sum(mos.values())-sum(resistor.values())
    assert abs(residual)<=1e-12, residual
    return dict(mos_absorbed_w=mos,resistor_absorbed_w=resistor,port_delivered_signed_w=delivered,
                net_port_delivered_w=sum(delivered.values()),balance_residual_w=residual,
                supply_delivered_a=supply_current,input_delivered_a=input_current,
                interpretation='Stationary OP terminal power including body/gate currents. Capacitor DC current and zero-volt probe power are zero. This is not measured transient energy or a pass-device model.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-zip',type=Path,required=True);p.add_argument('--measurements',type=Path,default=ROOT/'data/summary.json');p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    plan=json.loads((ROOT/'plan.json').read_text());measured=json.loads(a.measurements.read_text());prior=json.loads((ROOT/'data/prior.json').read_text());assert measured['plan_sha256']==hashlib.sha256((ROOT/'plan.json').read_bytes()).hexdigest()
    case,oldop,f,h=read_source(a.source_zip);new=measured['cases']['D1-100k-load'];connections=plan['cases']['D1-100k-load']['connections'];assert connections==case['connections']
    rows=[]
    for name,op,rload,gain,native_bw in [('known1Mohm',oldop,1e6,prior['known']['gain_one_khz'],prior['known']['bandwidth_hz']),('observed100kohm',new['op'],1e5,new['gain_one_khz'],new['bandwidth']['hz'])]:
        n=explicit_matrix(op,rload);s=stamp(op,connections,rload,plan['fixed']['rbias_ohm'])
        for key in ['gain','output_transimpedance_ohm']:assert np.isclose(n[key],s[key],rtol=1e-12)
        fp=1/(2*np.pi*plan['fixed']['cload_f']*n['output_transimpedance_ohm'])
        model_bw=fp*np.sqrt(10**.3*(1+(1000/fp)**2)-1)
        rows.append(dict(case=name,role='Retrospective calculation at this observed OP',network=n,independent_stamp=s,
                         actual_rload_ohm=rload,output_v=op['v(out)'],error_v=op['v(out)']-op['v(source)'],
                         load_current_a=op['v(out)']/rload,tail_delivered_a=-op['i(VtM5d)'],
                         p2_source_a=-op['i(VtM2d)'],n4_sink_a=op['i(VtM4d)'],
                         input_common_mode_v=(op['v(inp)']+op['v(inn)'])/2,input_difference_v=op['v(inp)']-op['v(inn)'],
                         tail_v=op['v(tail)'],input_p1_gm_s=derivatives(op,1)[0],feedback_p2_gm_s=derivatives(op,2)[0],
                         native_gain_one_khz=gain,native_bandwidth_hz=native_bw,
                         gain_relative_error_to_native_real=n['gain']/gain[0]-1,
                         single_output_pole_hz=float(fp),single_output_minus3db_hz=float(model_bw),
                         single_output_bw_relative_error=model_bw/native_bw-1,
                         power=power_balance(op,connections,rload,plan['fixed']['rbias_ohm'])))
    old,new=rows
    result=dict(schema='learning.initial-d1-dc-load.network.v1',status='ALGEBRA_AND_POWER_CHECKS_PASS',
                role='Retrospective interpretation using the two already observed operating points; original prior forecast remains separate and unchanged.',cases=rows,
                changes=dict(measured_bandwidth_ratio=new['native_bandwidth_hz']/old['native_bandwidth_hz'],
                             measured_gain_ratio=new['native_gain_one_khz'][0]/old['native_gain_one_khz'][0],
                             load_power_increase_w=new['power']['resistor_absorbed_w']['Rload']-old['power']['resistor_absorbed_w']['Rload'],
                             supply_current_increase_a=new['power']['supply_delivered_a']-old['power']['supply_delivered_a'],
                             effective_conductance_increase_s=new['network']['effective_output_conductance_s']-old['network']['effective_output_conductance_s'],
                             externally_added_conductance_s=9e-6),
                assumptions=['Local MOS gm/gmb/gds are taken at each actual OP; all body derivatives and finite load/bias resistance are retained.','The gain network neglects leakage derivatives and intrinsic capacitances; the absolute pole estimate uses the external5pF alone.','Terminal-power balance includes actual leakage currents even though their derivatives are omitted from the local network.','Both new gain/pole values use already observed derivatives. They are explanations, not a replacement for the frozen old-OP forecast or measured response.'],
                limits='One fixed receiver load change; no transient/overload/stability/noise/population or LDO qualification. The network does not isolate a single transistor parameter as the sole cause.',measurements_sha256=hashlib.sha256(a.measurements.read_bytes()).hexdigest())
    a.output.mkdir(parents=True);(a.output/'network.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    table=[dict(case=r['case'],rload_ohm=r['actual_rload_ohm'],output_v=r['output_v'],error_v=r['error_v'],load_current_a=r['load_current_a'],tail_current_a=r['tail_delivered_a'],source_p2_a=r['p2_source_a'],sink_n4_a=r['n4_sink_a'],supply_current_a=r['power']['supply_delivered_a'],load_power_w=r['power']['resistor_absorbed_w']['Rload'],local_gain=r['network']['gain'],native_gain_real=r['native_gain_one_khz'][0],output_transimpedance_ohm=r['network']['output_transimpedance_ohm'],native_bandwidth_hz=r['native_bandwidth_hz']) for r in rows]
    with (a.output/'network.csv').open('x',newline='') as file:
        writer=csv.DictWriter(file,list(table[0]),lineterminator='\n');writer.writeheader();writer.writerows(table)
    print(json.dumps(dict(status=result['status'],changes=result['changes'],cases=[{k:v for k,v in r.items() if k not in ['power','independent_stamp']} for r in rows]),indent=2))


if __name__=='__main__':main()
