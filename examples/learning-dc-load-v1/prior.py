"""MIT: prior-only initial-D1 affine receiving-load forecast; no acquisition."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np

ROOT = Path(__file__).resolve().parent
NEW_RLOAD = 100000.0


def read_source(path):
    source = json.loads((ROOT/'data/source.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == source['package_sha256']
    with zipfile.ZipFile(path) as z:
        for name, expected in source['selected_zip_members'].items():
            assert hashlib.sha256(z.read('example/'+name)).hexdigest() == expected
        case = json.loads(z.read('example/data/interface-inputs.json'))['cases']['D1-5p']
        with np.load(io.BytesIO(z.read('example/data/'+case['array_file'])), allow_pickle=False) as a:
            op_array, ac_array = a['analysis1'], a['analysis3']
    op_recipe, ac_recipe = case['analyses'][1]['recipe'], case['analyses'][3]['recipe']
    assert op_array.shape == (1, 1+len(op_recipe['vectors']))
    op = {n: float(op_array[0, 1+i]) for i,n in enumerate(op_recipe['vectors'])}
    assert ac_array.shape == (361, 1+2*len(ac_recipe['vectors']))
    values = {n: ac_array[:,1+2*i]+1j*ac_array[:,2+2*i] for i,n in enumerate(ac_recipe['vectors'])}
    assert np.isfinite(op_array).all() and np.isfinite(ac_array).all()
    assert op['@rload[resistance]'] == 1e6 and op['@cload[capacitance]'] == 5e-12
    return case, op, ac_array[:,0], values['v(out)']/values['v(source)']


def derivatives(op, index):
    return [op[f'@m.xm{index}.mapm045_vtg_core[{v}]'] for v in ('gm','gmbs','gds')]


def explicit_matrix(op, rload):
    gm1, gb1, gd1 = derivatives(op,1)
    gm2, gb2, gd2 = derivatives(op,2)
    gm3, _, gd3 = derivatives(op,3)
    gm4, _, gd4 = derivatives(op,4)
    _, _, gd5 = derivatives(op,5)
    s1,s2 = gm1+gb1+gd1,gm2+gb2+gd2
    # Unknown tail, mirror, out; feedback sets the second input increment=out.
    m = np.array([[s1+s2+gd5,-gd1,-gd2-gm2],
                  [-s1,gd1+gm3+gd3,0.],
                  [-s2,gm4,gd2+gd4+1/rload+gm2]])
    b = np.array([gm1,-gm1,0.])
    h = float(np.linalg.solve(m,b)[2])
    z = float(np.linalg.solve(m,np.array([0.,0.,1.]))[2])
    # Eliminate the first two nodes to form an independent output scalar equation.
    go = m[2,2]-m[2,:2]@np.linalg.solve(m[:2,:2],m[:2,2])
    drive = b[2]-m[2,:2]@np.linalg.solve(m[:2,:2],b[:2])
    assert np.isclose(z,1/go,rtol=1e-12) and np.isclose(h,drive/go,rtol=1e-12)
    return dict(gain=h,output_transimpedance_ohm=z,effective_output_conductance_s=float(go))


def stamp(op, connections, rload, rbias):
    # Separate six-MOS/four-node construction, retaining the finite bias node.
    names=['tail','mirror','out','bias'];m=np.zeros((4,4));b=np.zeros(4)
    for name, c in connections.items():
        gm,gb,gd=derivatives(op,int(name.lstrip('xm')))
        nodes={k:('out' if v=='inn' else v) for k,v in c['nodes'].items()}
        for terminal,sign in [('d',1.),('s',-1.)]:
            if nodes[terminal] not in names:continue
            row=names.index(nodes[terminal])
            for port,value in [('d',gd),('g',gm),('s',-gd-gm-gb),('b',gb)]:
                node=nodes[port]
                if node in names:m[row,names.index(node)]+=sign*value
                elif node=='inp':b[row]-=sign*value
                else:assert node in ('0','vdd')
    m[3,3]+=1/rbias;m[2,2]+=1/rload
    return dict(gain=float(np.linalg.solve(m,b)[2]),output_transimpedance_ohm=float(np.linalg.solve(m,[0.,0.,1.,0.])[2]))


def bw(f,h):
    db=20*np.log10(abs(h)/abs(h[0]));i=int(np.flatnonzero((db[:-1]>-3)&(db[1:]<=-3))[0])
    return float(np.exp(np.log(f[i])+(-3-db[i])/(db[i+1]-db[i])*np.log(f[i+1]/f[i])))


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-zip',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists()
    case,op,f,h=read_source(a.source_zip);old_r=op['@rload[resistance]'];cap=op['@cload[capacitance]']
    e=explicit_matrix(op,old_r);s=stamp(op,case['connections'],old_r,op['@rbias[resistance]'])
    for key in ('gain','output_transimpedance_ohm'):assert np.isclose(e[key],s[key],rtol=1e-12)
    delta=1/NEW_RLOAD-1/old_r;factor=1+delta*e['output_transimpedance_ohm'];new=explicit_matrix(op,NEW_RLOAD)
    assert np.isclose(new['gain'],e['gain']/factor,rtol=1e-12)
    pole=1/(2*np.pi*cap*e['output_transimpedance_ohm'])
    # Analytic -3.000-dB crossing relative to 1kHz for the one-pole approximation.
    exact_bw=pole*np.sqrt(10**.3*(1+(1000/pole)**2)-1)
    result=dict(schema='learning.initial-d1-dc-load.prior.v1',role='Prior-only forecast from known central D1 observations; saved before the one new target.',
                source_sha256=hashlib.sha256(a.source_zip.read_bytes()).hexdigest(),old_case='D1-5p',old_realization_id=case['original_realization_id'],
                known=dict(observed_op=op,gain_one_khz=[float(h[0].real),float(h[0].imag)],bandwidth_hz=bw(f,h)),
                old_network=e,independent_stamp=s,
                baseline_model_checks=dict(gain_relative_error_to_native_real=e['gain']/h[0].real-1,single_output_pole_hz=float(pole),single_output_pole_minus3db_hz=float(exact_bw),single_pole_bw_relative_error=exact_bw/bw(f,h)-1),
                forecast=dict(rload_ohm=NEW_RLOAD,added_conductance_s=delta,load_factor=factor,output_v=op['v(out)']/factor,
                              source_to_output_error_v=op['v(out)']/factor-op['v(source)'],resistor_current_a=op['v(out)']/factor/NEW_RLOAD,
                              gain_from_old_network=e['gain']/factor,gain_from_old_native=[float(h[0].real/factor),float(h[0].imag/factor)],
                              bandwidth_scaled_from_old_native_hz=bw(f,h)*factor,single_output_pole_hz=float(pole*factor)),
                assumptions=['Old observed DC output is the affine expansion point; gm/gmb/gds are held at that old OP.','Fixed supply/bias and negligible gate/leakage derivatives in the local network; body derivatives are retained.','A single external-output capacitor gives a first pole estimate. Scaling observed old bandwidth is an approximation, not a return-ratio or stability measurement.','No native100kohm target has been read. The new biased state may change derivatives; keep any forecast discrepancy.'])
    a.output.mkdir(parents=True);(a.output/'prior.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('known','assumptions')},indent=2))


if __name__=='__main__':main()
