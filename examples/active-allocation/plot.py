"""MIT: ordinary figures from compact active-allocation data, without SPICE."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES = ['A37', 'AL37W', 'AR37', 'AI37', 'AB37', 'R37']
LABELS = ['A37 baseline', 'AL37W: load', 'AR37: reference', 'AI37: input', 'AB37: balanced', 'R37 baseline']
COLORS = dict(zip(NAMES, ['#777777', '#D55E00', '#CC79A7', '#E69F00', '#009E73', '#0072B2']))


def main():
    p = argparse.ArgumentParser(); p.add_argument('--data', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.hashsalt': 'ads-allocation-v1'})
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
    ax.set_yticks(ypos, LABELS); ax.invert_yaxis(); ax.set(xlabel='Input-referred RMS noise, 1 kHz–1 MHz (µV)', title='Nominal noise cost at 37.5 µA'); grid(ax)
    save(fig, 'nominal-noise')
    rows = read('confirmation'); stats = json.loads((a.data/'statistics.json').read_text())
    for condition, title in [('nominal', 'Fresh physical coordinates at 1 V'), ('supply_097', 'Same physical coordinates at 0.97 V')]:
        fig, ax = plt.subplots(figsize=(6.4, 4.4), layout='constrained')
        for n in NAMES:
            values = sorted(float(r['error_v'])*1000 for r in rows if r['candidate'] == n and r['condition'] == condition and r['numerical_usable'] == 'True' and r['terminal_supported'] == 'True')
            ax.step(values, np.arange(1, len(values)+1)/len(values), where='post', color=COLORS[n], label=n,
                lw=1.8 if n in ['AL37W', 'R37'] else 1.1, ls='--' if n in ['A37', 'AI37'] else '-')
        ax.axvspan(-25, 25, color='#999999', alpha=.13, label='±25 mV center band'); ax.axvline(0, color='#444444', lw=.7)
        ax.set(xlabel='Signed output error about fixed 0.5 V (mV)', ylabel='Empirical cumulative fraction', ylim=(0, 1.02), title=title)
        ax.legend(frameon=False, fontsize=9, ncol=2, loc='lower right'); ax.grid(alpha=.15)
        save(fig, 'confirmation-'+('nominal' if condition == 'nominal' else 'supply'))
    pair_names = [('A37', 'AL37W'), ('AR37', 'AL37W'), ('AI37', 'AL37W'), ('AB37', 'AL37W'), ('AL37W', 'R37')]
    fig, ax = plt.subplots(figsize=(6.8, 4.8), layout='constrained')
    for c, color, offset, label in [('nominal', '#0072B2', -.12, '1 V'), ('supply_097', '#D55E00', .12, '0.97 V')]:
        selected = [next(r for r in stats['contrasts'] if r['condition'] == c and (r['left'], r['right']) == pair) for pair in pair_names]
        means = np.array([r['mean_mae_reduction_v']*1000 for r in selected]); intervals = np.array([r['interval_v'] for r in selected])*1000
        ax.errorbar(means, np.arange(5)+offset, xerr=[means-intervals[:, 0], intervals[:, 1]-means], fmt='o', capsize=3, color=color, label=label)
    for i, value in enumerate([8, 5, 5, 5, 20]):
        ax.plot([value, value], [i-.32, i+.32], color='#666666', ls=':', lw=1.8, label='Frozen useful-effect threshold' if i == 0 else None)
    ax.axvline(0, color='#444444', lw=.8); ax.set_yticks(np.arange(5), [a+' → '+b for a, b in pair_names]); ax.invert_yaxis()
    ax.set(xlabel='Mean absolute-error reduction favoring right candidate (mV)', title='Finite matched choices and the resistor baseline')
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
    print('Generated seven SVG/PNG figures from public compact data; no SPICE.')


if __name__ == '__main__':
    main()
