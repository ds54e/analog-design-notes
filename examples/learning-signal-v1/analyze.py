"""MIT: six known nominal cycles, finite harmonics and measured input ports.

Pure saved-data analysis. No native acquisition, fitting to unseen targets,
population claim, recentering or change to the old eight-cycle THD definition.
"""
import argparse
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np

VERSION = 'learning.nonlinearity.v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fourier(values):
    """A single cycle, excluding its repeated endpoint; peak complex phasors."""
    assert values.shape == (1000,) and np.isfinite(values).all()
    ft = np.fft.rfft(values)/len(values)
    theta = 2*np.pi*np.arange(len(values))/len(values)
    design = np.column_stack([np.ones(len(values))]+[
        f(n*theta) for n in range(1, 11) for f in [np.cos, np.sin]])
    fit = np.linalg.lstsq(design, values, rcond=None)[0]
    peak = 2*ft[1:11]
    fitted = fit[1::2]-1j*fit[2::2]
    discrepancy = max(abs(fit[0]-ft[0].real), float(np.max(abs(peak-fitted))))
    assert discrepancy < 1e-12
    return float(ft[0].real), peak, discrepancy


def local_network(name, cfg, op, nominal):
    """Held-supply local circuit equations; leakage/capacitance omitted."""
    d = cfg['candidates'][name]; p = d['parameters']
    n = [u for u in d['devices'] if u['polarity'] == 'n']
    def bank(units, quantity):
        return sum(op['@m.'+'.'.join(u['path'])+'.mapm045_vtg_core['+quantity+']'] for u in units)
    gm, gb, gd = (bank(n, x) for x in ['gm','gmbs','gds'])
    degeneration = 1+p['rs']*(gm+gb+gd)
    ge = gm/degeneration
    load = [u for u in d['devices'] if u['polarity'] == 'p' and u['path'][0].startswith('xp')]
    go = 1/p['rload']+gd/degeneration+bank(load, 'gds')
    if p.get('rd'):
        go += 1/p['rd']
    rt = p['rsrc']+p.get('rin', 0)
    if 'rf' in p:
        alpha = p['rf']/(p['rf']+rt); beta = rt/(p['rf']+rt)
        # f retains the actual Rf from out to an independently held gate.
        # It is a diagnostic gate-driven fixture, not a matched competitor.
        f1 = -(ge-1/p['rf'])/(go+1/p['rf'])
        denominator = 1-beta*f1
        h1 = alpha*f1/denominator
        curvature_factor = alpha**2/denominator**3
    else:
        alpha, beta, f1, denominator = 1., 0., -ge/go, 1.
        h1, curvature_factor = f1, 1.
    return dict(candidate=name,gm_bank_s=gm,gmb_bank_s=gb,gds_bank_s=gd,
        source_degeneration_factor=degeneration,effective_gm_s=ge,output_conductance_s=go,
        alpha=alpha,beta=beta,gate_driven_fixture_slope=f1,static_feedback_denominator=denominator,
        circuit_slope=h1,old_ac_real_gain=float(nominal['gain']),
        slope_relative_discrepancy=h1/float(nominal['gain'])-1,
        gate_fraction=alpha+beta*h1,curvature_factor_relative_to_fixture=curvature_factor,
        local_thd_coefficient_factor_relative_to_fixture=alpha/denominator**2)


def calculate(feedback_zip, here=None):
    here = here or Path(__file__).resolve().parent
    for name, digest in json.loads((here/'data/manifest.json').read_text())['files'].items():
        path = Path(name)
        assert not path.is_absolute() and '..' not in path.parts
        assert sha(here/'data'/path) == digest, name
    source = json.loads((here/'data/source.json').read_text())
    assert sha(feedback_zip) == source['old_zip_sha256']
    with zipfile.ZipFile(feedback_zip) as z:
        selected = {name:z.read('example/'+name) for name in source['selected_members_sha256']}
    for name, body in selected.items():
        assert hashlib.sha256(body).hexdigest() == source['selected_members_sha256'][name]
    manifest = json.loads(selected['data/manifest.json'])['files']
    for name, body in selected.items():
        if name.startswith('data/') and name != 'data/manifest.json':
            assert hashlib.sha256(body).hexdigest() == manifest[name[5:]]
    def read(name):
        return list(csv.DictReader(io.StringIO(selected['data/'+name].decode())))
    cfg = json.loads(selected['data/configuration.json'])
    nominal = {r['candidate']:r for r in read('nominal.csv')}
    op = {n:{r['vector']:float(r['value']) for r in read('nominal-op.csv') if r['candidate'] == n} for n in ['R75','F75L']}
    static_source = json.loads((here/'data/static-source.json').read_text())
    assert sha(here/'data/dc-points.csv') == static_source['dc_points_sha256']
    with (here/'data/dc-points.csv').open() as f:
        points = list(csv.DictReader(f))
    assert len(points) == len(static_source['cases']) == 8
    for case in static_source['cases']:
        assert sha(here/case['circuit']) == case['circuit_sha256']
        row = next(r for r in points if r['case'] == case['case'])
        for field, vector in [('input_v','v(in)'),('output_v','v(out)'),('gate_v','v(gate)'),('source_v','v(source)')]:
            assert float(row[field]) == case['original_op'][vector]
        assert all(u['raw'] == [0.,0.] for u in case['physical_units'])
    static = []
    for name in ['R75','F75L']:
        values = {float(r['input_delta_v']):float(r['output_v']) for r in points if r['candidate'] == name}
        assert values[0.] == float(nominal[name]['out_v'])
        for step in sorted((x for x in values if x > 0),reverse=True):
            slope = (values[step]-values[-step])/(2*step)
            curvature = (values[step]-2*values[0.]+values[-step])/step**2
            static.append(dict(candidate=name,input_step_v=step,slope_v_per_v=slope,
                second_derivative_per_v=curvature,center_output_v=values[0.],
                positive_output_v=values[step],negative_output_v=values[-step],
                quadratic_input_peak_for_two_percent_thd_v=4*abs(slope)*.02/abs(curvature)))
    assert sha(here/'data/port-cycles.npz') == source['port_cycles_sha256']
    spectra = []; summaries = []; waves = {}; all_fft_errors = []
    with np.load(here/'data/port-cycles.npz', allow_pickle=False) as arrays:
        assert set(arrays.files) == {r['case'] for r in source['cases']}
        for case in source['cases']:
            label = case['case']; name = case['candidate']; p = case['parameters']
            assert sha(here/case['circuit']) == case['circuit_sha256']
            assert all(u['raw'] == [0., 0.] for u in case['physical_units'])
            assert case['original_metrics']['sine_stationary'] == 'True'
            data = arrays[label]; assert data.shape == (1001, len(source['columns']))
            assert np.isfinite(data).all()
            full = {c:data[:,j].copy() for j,c in enumerate(source['columns'])}
            assert np.array_equal(full['time_s'], np.linspace(.0009,.001,1001))
            # A new one-cycle definition. The old eight-cycle metrics are kept.
            wave = {c:v[:-1] for c,v in full.items()}; waves[label] = wave
            coefficients = {}
            for key in ['input_v','output_v','gate_v','source_v','input_delivered_current_a']:
                mean, peak, difference = fourier(wave[key]); all_fft_errors.append(difference)
                coefficients[key] = (mean, peak)
            _, vgs, err = fourier(wave['gate_v']-wave['source_v']); all_fft_errors.append(err)
            input_peak = coefficients['input_v'][1][0]; amplitude = abs(input_peak)
            assert amplitude > 1e-3
            output = coefficients['output_v'][1]
            referenced = output*np.exp(-1j*np.arange(1,11)*np.angle(input_peak))
            amps = abs(output); thd = np.linalg.norm(amps[1:])/amps[0]
            for n in range(1,11):
                spectra.append(dict(case=label,candidate=name,amplitude_peak_v=case['amplitude_peak_v'],
                    refined=case['refined'],harmonic=n,peak_v=float(amps[n-1]),
                    fraction_of_fundamental=float(amps[n-1]/amps[0]),
                    in_phase_with_input_cosine_v=float(referenced[n-1].real),
                    negative_sine_quadrature_v=float(referenced[n-1].imag)))
            dc = float(case['original_metrics']['out_v'])
            rt = p['rsrc']+p.get('rin',0)
            exact_gate = wave['input_v']-rt*wave['input_delivered_current_a']
            exact_source = -p['rs']*wave['nmos_source_current_a']
            assert np.max(abs(exact_gate-wave['gate_v'])) < 1e-11
            assert np.max(abs(exact_source-wave['source_v'])) < 1e-11
            if 'rf' in p:
                divider = (p['rf']*wave['input_v']+rt*wave['output_v'])/(p['rf']+rt)
                corrected = divider-(rt*p['rf']/(rt+p['rf']))*wave['nmos_gate_current_a']
            else:
                divider = wave['input_v']; corrected = divider-rt*wave['nmos_gate_current_a']
            assert np.max(abs(corrected-wave['gate_v'])) < 1e-10
            divider_error = divider-wave['gate_v']
            summaries.append(dict(case=label,candidate=name,amplitude_peak_v=case['amplitude_peak_v'],refined=case['refined'],
                observed_input_fundamental_peak_v=float(amplitude),fundamental_peak_v=float(amps[0]),
                fundamental_gain_magnitude=float(amps[0]/amplitude),
                output_inversion_phase_error_degrees=float(np.degrees(np.angle(-output[0]/input_peak))),
                original_eight_cycle_thd_fraction=float(case['original_metrics']['thd_fraction']),
                one_cycle_thd_fraction=float(thd),thd_difference_from_original_fraction=float(thd)-float(case['original_metrics']['thd_fraction']),
                second_harmonic_peak_v=float(amps[1]),third_harmonic_peak_v=float(amps[2]),
                second_harmonic_distortion_power_fraction=float(amps[1]**2/np.sum(amps[1:]**2)),
                mean_output_v=coefficients['output_v'][0],independent_op_output_v=dc,
                mean_shift_from_op_v=coefficients['output_v'][0]-dc,
                second_harmonic_in_phase_v=float(referenced[1].real),
                apparent_quasistatic_second_derivative_per_v=float(4*referenced[1].real/amplitude**2),
                mean_shift_minus_in_phase_h2_v=coefficients['output_v'][0]-dc-float(referenced[1].real),
                output_endpoint_repeat_difference_v=float(full['output_v'][-1]-full['output_v'][0]),
                gate_fundamental_peak_v=float(abs(coefficients['gate_v'][1][0])),
                source_fundamental_peak_v=float(abs(coefficients['source_v'][1][0])),
                vgs_fundamental_peak_v=float(abs(vgs[0])),gate_over_input_fundamental=float(abs(coefficients['gate_v'][1][0])/amplitude),
                input_current_min_a=float(wave['input_delivered_current_a'].min()),
                input_current_max_a=float(wave['input_delivered_current_a'].max()),
                input_current_mean_a=coefficients['input_delivered_current_a'][0],
                zero_gate_current_divider_max_error_v=float(np.max(abs(divider_error))),
                zero_gate_current_divider_mean_removed_max_error_v=float(np.max(abs(divider_error-divider_error.mean()))),
                exact_port_gate_max_error_v=float(np.max(abs(exact_gate-wave['gate_v']))),
                exact_port_source_max_error_v=float(np.max(abs(exact_source-wave['source_v'])))))
    extrapolations = []
    for name in ['R75','F75L']:
        low = next(r for r in summaries if r['case'] == name+'-10m')
        high = next(r for r in summaries if r['case'] == name+'-20m')
        curvature = low['apparent_quasistatic_second_derivative_per_v']
        dc_curve = next(r for r in static if r['candidate'] == name and r['input_step_v'] == .0001)
        a = high['amplitude_peak_v']; slope = float(nominal[name]['gain'])
        mean = curvature*a*a/4
        extrapolations.append(dict(candidate=name,role='Retrospective 10-to-20-mV extrapolation; both targets already known',
            input_source_case=name+'-10m',comparison_case=name+'-20m',apparent_second_derivative_per_v=curvature,
            linear_slope_from_old_ac=slope,predicted_second_harmonic_peak_v=abs(mean),
            observed_second_harmonic_peak_v=high['second_harmonic_peak_v'],
            second_harmonic_relative_difference=abs(mean)/high['second_harmonic_peak_v']-1,
            predicted_mean_shift_v=mean,observed_mean_shift_v=high['mean_shift_from_op_v'],
            mean_shift_error_v=mean-high['mean_shift_from_op_v'],
            predicted_quadratic_thd_fraction=abs(curvature)*a/(4*abs(slope)),
            observed_one_cycle_thd_fraction=high['one_cycle_thd_fraction'],
            static_second_derivative_per_v=dc_curve['second_derivative_per_v'],
            static_quadratic_mean_shift_v=dc_curve['second_derivative_per_v']*a*a/4,
            static_quadratic_thd_fraction=abs(dc_curve['second_derivative_per_v'])*a/(4*abs(dc_curve['slope_v_per_v'])),
            sine_10m_curvature_relative_to_static=curvature/dc_curve['second_derivative_per_v']-1,
            sine_20m_curvature_relative_to_static=high['apparent_quasistatic_second_derivative_per_v']/dc_curve['second_derivative_per_v']-1))
    refinements = []
    for name in ['R75','F75L']:
        base = next(r for r in summaries if r['case'] == name+'-20m')
        fine = next(r for r in summaries if r['case'] == name+'-20m-refined')
        refinements.append(dict(candidate=name,thd_relative_change=fine['one_cycle_thd_fraction']/base['one_cycle_thd_fraction']-1,
            h2_relative_change=fine['second_harmonic_peak_v']/base['second_harmonic_peak_v']-1,
            mean_output_difference_v=fine['mean_output_v']-base['mean_output_v'],
            gate_peak_relative_change=fine['gate_fundamental_peak_v']/base['gate_fundamental_peak_v']-1))
    result = dict(schema=VERSION,status='PASS',role='Retrospective analysis of six known nominal cycles',native_acquisitions=0,
        old_zip_sha256=source['old_zip_sha256'],source_sha256=sha(here/'data/source.json'),
        measurement='One final cycle, first1000 of1001 uniformly resampled points; peak harmonics1..10; phase relative to measured input. Original eight-cycle THD/stationarity retained separately.',
        fft_vs_least_squares_max_absolute_difference=float(max(all_fft_errors)),cases=summaries,
        static_source_sha256=sha(here/'data/static-source.json'),static_finite_differences=static,
        retrospective_extrapolations=extrapolations,known_timestep_controls=refinements,
        local_network=[local_network(n,cfg,op[n],nominal[n]) for n in ['R75','F75L']],
        limits=['Sinusoid-inferred apparent curvature and the separately selected static finite differences are different estimates; both use known observations.',
                'R75 has two saved static steps; F75L has one. No new static refinement or native condition was acquired.',
                'One public cycle does not independently establish stationarity; the saved original eight-cycle assessment remains the source.',
                'The quasistatic polynomial omits memory and higher orders; discrepancies and original physical failures remain.',
                'R75 and F75L differ in complete design and resource allocation; their comparison does not isolate Rf alone.'])
    return result, spectra, waves


def main():
    p = argparse.ArgumentParser(); p.add_argument('--feedback-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    result, spectra, _ = calculate(a.feedback_zip)
    for name, rows in [('harmonics',spectra),('cycle-summary',result['cases']),
                       ('quadratic-extrapolation',result['retrospective_extrapolations']),
                       ('local-network',result['local_network']),('static-curvature',result['static_finite_differences'])]:
        with (a.output/(name+'.csv')).open('x',newline='') as f:
            writer = csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],cases=len(result['cases']),harmonics=len(spectra),output=str(a.output))))


if __name__ == '__main__':
    main()
