"""MIT: ordinary figures from compact finite-feedback data, without SPICE."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES = ['A75', 'R75', 'F75S', 'F75L']
LABELS = ['A75: active baseline', 'R75: resistor baseline', 'F75S: small PMOS', 'F75L: large PMOS']
COLORS = dict(zip(NAMES, ['#777777', '#0072B2', '#D55E00', '#009E73']))


def main():
    p = argparse.ArgumentParser(); p.add_argument('--data', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.hashsalt': 'ads-feedback-v1'})
    def read(name):
        with (a.data/(name+'.csv')).open() as f: return list(csv.DictReader(f))
    def save(fig, name):
        fig.savefig(a.output/(name+'.svg'), metadata={'Date': None}, bbox_inches='tight')
        fig.savefig(a.output/(name+'.png'), dpi=160, bbox_inches='tight'); plt.close(fig)
    def grid(ax):
        ax.grid(axis='x', alpha=.2); ax.set_axisbelow(True)
    budget = read('source-budget'); ypos = np.arange(len(NAMES))
    fig, ax = plt.subplots(figsize=(6.6, 4.6), layout='constrained'); left = np.zeros(len(NAMES))
    for role, color in [('input', '#0072B2'), ('reference', '#009E73'), ('load', '#D55E00')]:
        values = np.array([sum(float(r['variance_v2'])*1e6 for r in budget if r['candidate'] == n and r['role'] == role) for n in NAMES])
        ax.barh(ypos, values, left=left, label=role.capitalize(), color=color, height=.62); left += values
    ax.set_yticks(ypos, LABELS); ax.invert_yaxis(); ax.set(xlabel='Nominal first-order output variance (mV²)', title='Where the modeled error is generated')
    ax.legend(frameon=False, ncol=3, loc='lower right'); grid(ax); save(fig, 'source-budget')
    nominal = {r['candidate']: r for r in read('nominal')}
    fig, ax = plt.subplots(figsize=(6.5, 4.4), layout='constrained')
    ax.barh(ypos, [float(nominal[n]['input_noise_v'])*1e6 for n in NAMES], color=[COLORS[n] for n in NAMES], height=.62)
    ax.set_yticks(ypos, LABELS); ax.invert_yaxis(); ax.set(xlabel='Input-referred RMS noise, 1 kHz–1 MHz (µV)', title='Nominal noise cost at 75 µA'); grid(ax)
    save(fig, 'nominal-noise')
    rows = read('confirmation'); stats = json.loads((a.data/'statistics.json').read_text())
    for condition, title in [('nominal', 'Fresh physical coordinates at 1 V'), ('supply_097', 'Same physical coordinates at 0.97 V')]:
        fig, ax = plt.subplots(figsize=(6.4, 4.4), layout='constrained')
        for n in NAMES:
            values = sorted(float(r['error_v'])*1000 for r in rows if r['candidate'] == n and r['condition'] == condition and r['numerical_usable'] == 'True' and r['terminal_supported'] == 'True')
            ax.step(values, np.arange(1, len(values)+1)/len(values), where='post', color=COLORS[n], label=n,
                lw=1.8 if n in ['F75L', 'R75'] else 1.1, ls='--' if n in ['A75', 'F75S'] else '-')
        ax.axvspan(-25, 25, color='#999999', alpha=.13, label='±25 mV center band'); ax.axvline(0, color='#444444', lw=.7)
        ax.set(xlabel='Signed output error about fixed 0.5 V (mV)', ylabel='Empirical cumulative fraction', ylim=(0, 1.02), title=title)
        ax.legend(frameon=False, fontsize=9, ncol=2, loc='lower right'); ax.grid(alpha=.15)
        save(fig, 'confirmation-'+('nominal' if condition == 'nominal' else 'supply'))
    pair_names = [('A75', 'F75S'), ('A75', 'F75L'), ('F75S', 'F75L'), ('F75L', 'R75')]
    fig, ax = plt.subplots(figsize=(6.8, 4.8), layout='constrained')
    for c, color, offset, label in [('nominal', '#0072B2', -.12, '1 V'), ('supply_097', '#D55E00', .12, '0.97 V')]:
        selected = [next(r for r in stats['contrasts'] if r['condition'] == c and (r['left'], r['right']) == pair) for pair in pair_names]
        means = np.array([r['mean_mae_reduction_v']*1000 for r in selected]); intervals = np.array([r['interval_v'] for r in selected])*1000
        ax.errorbar(means, np.arange(4)+offset, xerr=[means-intervals[:, 0], intervals[:, 1]-means], fmt='o', capsize=3, color=color, label=label)
    for i, value in enumerate([20, 25, 4, 8]):
        ax.plot([value, value], [i-.32, i+.32], color='#666666', ls=':', lw=1.8, label='Frozen useful-effect threshold' if i == 0 else None)
    ax.axvline(0, color='#444444', lw=.8); ax.set_yticks(np.arange(4), [a+' → '+b for a, b in pair_names]); ax.invert_yaxis()
    ax.set(xlabel='Mean absolute-error reduction favoring right candidate (mV)', title='Where a useful accuracy benefit is established')
    ax.legend(frameon=False, fontsize=9, loc='upper left', bbox_to_anchor=(0, -.18), ncol=3); grid(ax); save(fig, 'paired-comparisons')
    fig, ax = plt.subplots(figsize=(6.5, 4.5), layout='constrained')
    for c, color, offset, label in [('nominal', '#0072B2', -.17, '1 V'), ('supply_097', '#D55E00', .17, '0.97 V')]:
        group = [stats['summaries'][n+'/'+c] for n in NAMES]
        ax.barh(ypos+offset, [100*r['dc_ac_current_pass']/r['assigned'] for r in group], height=.3, color=color, label=label)
    ax.set_yticks(ypos, LABELS); ax.invert_yaxis(); ax.set(xlabel='Assessed DC/AC/current pass fraction (%)', xlim=(0, 100), title='Accuracy improvement does not guarantee function')
    ax.legend(frameon=False, loc='upper right'); grid(ax); save(fig, 'assessed-function')
    fig, ax = plt.subplots(figsize=(6.5, 4.5), layout='constrained')
    for c, color, offset, label in [('nominal', '#0072B2', -.17, '1 V'), ('supply_097', '#D55E00', .17, '0.97 V')]:
        ax.barh(ypos+offset, [1000*stats['summaries'][n+'/'+c]['prediction_rmse_v'] for n in NAMES], height=.3, color=color, label=label)
    ax.set_yticks(ypos, LABELS); ax.invert_yaxis(); ax.set(xlabel='Saved pre-target prediction RMSE (mV)', title='Limit of the nominal linear voltage prediction')
    ax.legend(frameon=False, loc='lower right'); grid(ax); save(fig, 'prediction-error')
    dynamic_figures(read, save, nominal)
    print('Generated eleven SVG/PNG figures from public compact data; no SPICE.')


def dynamic_figures(read, save, nominal):
    dynamics = read('dynamics'); selected = read('selected-signal')
    fig, ax = plt.subplots(figsize=(6.5, 4.5), layout='constrained')
    for name in NAMES:
        group = sorted([r for r in dynamics if r['candidate'] == name and r['kind'] == 'sine' and r['refined'] == 'False'], key=lambda r: float(r['amplitude_peak_v']))
        x = [5]+[float(r['amplitude_peak_v'])*1000 for r in group]
        y = [float(nominal[name]['sine_thd_fraction'])*100]+[float(r['thd_fraction'])*100 for r in group]
        ax.plot(x, y, 'o-', color=COLORS[name], label=name+' nominal')
        sparse = [float(r['thd_fraction'])*100 for r in selected if r['candidate'] == name]
        offset = dict(A75=-.42, R75=-.14, F75S=.14, F75L=.42)[name]
        ax.scatter(np.full(len(sparse), 20+offset), sparse, marker='x', color=COLORS[name], s=35)
    ax.axhline(2, color='#555555', ls='--', label='Fixed 2% limit')
    ax.scatter([], [], marker='x', color='#222222', label='6 Local cases/design at 20 mV')
    ax.set(xlabel='Input sine peak amplitude (mV), 10 kHz', ylabel='Output THD (%)', title='More signal changes the useful choice', xticks=[5, 10, 20])
    ax.legend(frameon=False, fontsize=9, ncol=2); ax.grid(alpha=.15)
    save(fig, 'signal-distortion')
    fig, ax = plt.subplots(figsize=(6.5, 4.5), layout='constrained')
    for name in NAMES:
        values = [(1., float(nominal[name]['bandwidth_hz'])/1e6)]
        values += [(float(r['capacitance_f'])*1e12, float(r['bandwidth_hz'])/1e6) for r in dynamics if r['candidate'] == name and r['kind'] == 'capacitance']
        values.sort(); ax.loglog(*zip(*values), 'o-', color=COLORS[name], label=name)
    ax.axhline(1, color='#555555', ls='--', label='1-MHz requirement')
    ax.set(xlabel='External load capacitance (pF)', ylabel='First −3-dB bandwidth (MHz)', title='Feedback is faster than A75, slower than R75')
    ax.set_xticks([.2, 1, 8, 16], ['0.2', '1', '8', '16']); ax.legend(frameon=False); ax.grid(alpha=.15, which='both')
    save(fig, 'capacitive-load')
    waves = read('step-waveforms')
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.3), layout='constrained')
    for name in NAMES:
        group = [r for r in waves if r['candidate'] == name and r['refined'] == 'False' and float(r['capacitance_f']) == 8e-12]
        for ax, origin, end in [(axes[0], 2e-6, 7.9e-6), (axes[1], 8e-6, 14e-6)]:
            points = [r for r in group if origin <= float(r['time_s']) <= end]
            ax.plot([(float(r['time_s'])-origin)*1e6 for r in points], [float(r['output_v']) for r in points], color=COLORS[name], label=name)
            ax.set(xlabel='Time since ramp begins (µs)', ylabel='Output voltage (V)'); ax.grid(alpha=.15)
    axes[0].set_title('Input 0.43 → 0.47 V, 8-pF load'); axes[1].set_title('Input 0.47 → 0.43 V, 8-pF load')
    axes[0].legend(frameon=False, ncol=4, fontsize=9)
    save(fig, 'step-response')
    waves = read('sine-waveforms')
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.3), layout='constrained')
    for name in ['R75', 'F75L']:
        group = [r for r in waves if r['candidate'] == name and r['refined'] == 'False' and float(r['amplitude_peak_v']) == .02]
        for ax, column in [(axes[0], 'input_delivered_current_a'), (axes[1], 'supply_delivered_current_a')]:
            ax.plot([(float(r['time_s'])-.0009)*1e6 for r in group], [float(r[column])*1e6 for r in group], color=COLORS[name], label=name)
            ax.set(xlabel='Time within final 10-kHz cycle (µs)', ylabel='Delivered current (µA)'); ax.grid(alpha=.15)
    axes[0].axhline(0, color='#555555', lw=.7); axes[0].set_title('Input source: negative means absorption'); axes[1].set_title('Supply source: DC ceiling differs from signal peaks')
    axes[0].legend(frameon=False); save(fig, 'port-currents')


if __name__ == '__main__':
    main()
