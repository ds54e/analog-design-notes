"""MIT: five ordinary plots from compact sequencing observations; no SPICE."""
import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES = ['A75', 'R75', 'F75L']
COLORS = dict(zip(NAMES, ['#777777', '#0072B2', '#009E73']))


def main():
    p = argparse.ArgumentParser(); p.add_argument('--data', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.hashsalt': 'adn-sequencing-v1'})
    def read(name):
        with (a.data/(name+'.csv')).open() as f: return list(csv.DictReader(f))
    def save(fig, name):
        fig.savefig(a.output/(name+'.svg'), metadata={'Date': None}, bbox_inches='tight')
        fig.savefig(a.output/(name+'.png'), dpi=160, bbox_inches='tight'); plt.close(fig)
    def finish(ax): ax.grid(alpha=.18); ax.set_axisbelow(True)
    static = read('static'); histories = read('histories'); waves = read('waveforms')
    def trace(cid):
        rows = [r for r in waves if r['case_id'] == cid]
        assert rows
        return {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if k != 'case_id'}
    states = ['input_only', 'supply_only', 'both_on']; labels = ['Input only\nVDD=0, Vin=0.45 V', 'Supply only\nVDD=1, Vin=0 V', 'Both on\nVDD=1, Vin=0.45 V']
    fig, axes = plt.subplots(2, 1, figsize=(7, 6.8), layout='constrained')
    x = np.arange(3)
    for i, name in enumerate(NAMES):
        rows = [next(r for r in static if r['candidate'] == name and r['state'] == s) for s in states]
        offset = (i-1)*.24
        axes[0].bar(x+offset, [float(r['out_v']) for r in rows], width=.22, color=COLORS[name], label=name)
        axes[1].bar(x+offset, [float(r['input_delivered_a'])*1e6 for r in rows], width=.22, color=COLORS[name], label=name)
    axes[0].set(ylabel='Stationary output (V)', title='Zero volts means a driven/clamped source', ylim=(0, 1.15))
    axes[1].set(ylabel='Input current delivered (µA)', title='Negative input current means absorption')
    axes[0].legend(frameon=False, ncol=3); axes[1].axhline(0, lw=.8, color='#333333')
    for ax in axes: ax.set_xticks(x, labels, fontsize=9); finish(ax)
    save(fig, 'static-interface')

    fig, axes = plt.subplots(3, 1, figsize=(6.7, 8.5), layout='constrained')
    for ax, seq, title in zip(axes, ['input_first', 'supply_first', 'together'],
                              ['Input already at 0.45 V; supply rises', 'Supply already at 1 V; input rises', 'Both sources rise together']):
        for name in NAMES:
            w = trace(name+'-'+seq+'-fast-8p'); t = (w['time_s']-2e-6)*1e6; mask = (t >= -.1) & (t <= 6)
            ax.plot(t[mask], w['output_v'][mask], color=COLORS[name], label=name, lw=1.7)
        ax.axvspan(0, .01, color='#444444', alpha=.12); ax.axhline(.5, color='#666666', ls=':', lw=.8)
        ax.set(xlabel='Time since 10-ns ramp begins (µs)', ylabel='Output (V)', title=title, xlim=(-.1, 6), ylim=(-.025, .98)); finish(ax)
    axes[0].legend(frameon=False, ncol=3); fig.suptitle('Finite startup histories at 8 pF', fontsize=14)
    save(fig, 'startup-histories')

    fig, axes = plt.subplots(2, 1, figsize=(6.7, 6.8), layout='constrained')
    for ax, cid, title in [(axes[0], 'F75L-supply_first-fast-8p', 'F75L: supply first, then a 10-ns input ramp'),
                            (axes[1], 'F75L-together-fast-8p', 'F75L: both sources rise in 10 ns')]:
        w = trace(cid); t = (w['time_s']-2e-6)*1e6; mask = (t >= -.03) & (t <= .08)
        for key, label, color, style in [('input_delivered_a', 'Measured total', '#222222', '-'),
            ('input_series_term_a', '(Vin − Vo) / 110 kΩ', '#009E73', '--'),
            ('input_gate_term_a', '(100/110) × total gate current', '#D55E00', ':')]:
            ax.plot(t[mask], w[key][mask]*1e6, label=label, color=color, ls=style, lw=1.8)
        ax.axvspan(0, .01, color='#999999', alpha=.15); ax.axhline(0, color='#666666', lw=.6)
        ax.set(xlabel='Time since ramp begins (µs)', ylabel='Input current delivered (µA)', title=title); finish(ax)
    axes[0].legend(frameon=False, fontsize=9, loc='lower right')
    save(fig, 'input-current-paths')

    fig, axes = plt.subplots(2, 1, figsize=(7, 6.8), layout='constrained')
    labels = ['Input first', 'Supply first', 'Together']
    for i, name in enumerate(NAMES):
        group = [next(r for r in histories if r['case_id'] == name+'-'+s+'-fast-8p') for s in ['input_first', 'supply_first', 'together']]
        positions = x+(i-1)*.24
        axes[0].bar(positions, [float(r['Vdd_positive_delivery_peak_a'])*1e6 for r in group], width=.22, color=COLORS[name], label=name)
        axes[1].bar(positions, [float(r['Vin_positive_delivery_peak_a'])*1e6 for r in group], width=.22, color=COLORS[name])
        axes[1].bar(positions, [-float(r['Vin_absorption_peak_a'])*1e6 for r in group], width=.22, color=COLORS[name], hatch='//', alpha=.65)
    axes[0].set(ylabel='Peak supply delivery (µA)', title='8-pF load and 10-ns source ramp')
    axes[0].axhline(75, color='#555555', ls=':', lw=1.2, label='75 µA nominal DC reference')
    axes[0].legend(frameon=False, fontsize=9, ncol=2)
    axes[1].axhline(0, color='#333333', lw=.8)
    axes[1].set(ylabel='Input delivery / absorption (µA)', title='Hatching below zero: peak absorption')
    for ax in axes: ax.set_xticks(x, labels); finish(ax)
    save(fig, 'peak-demands')

    fig, axes = plt.subplots(3, 1, figsize=(6.7, 8.4), layout='constrained')
    for name in NAMES:
        w = trace(name+'-supply-cycle-8p'); t = w['time_s']*1e6
        axes[1].plot(t, w['output_v'], label=name, color=COLORS[name], lw=1.7)
        axes[2].plot(t, w['input_delivered_a']*1e6, label=name, color=COLORS[name], lw=1.7)
    axes[0].plot(t, w['supply_v'], color='#444444'); axes[0].set(ylabel='VDD (V)', title='Input held at 0.45 V; 1-µs supply ramps, 8-pF load')
    axes[1].set(ylabel='Output (V)', title='Returns to the observed on-state after the cycle', ylim=(-.025, .68)); axes[1].legend(frameon=False, ncol=3, loc='upper center')
    axes[2].set(ylabel='Input delivery (µA)', title='Input still delivers current during the zero-supply hold')
    for ax in axes:
        ax.set(xlabel='Time (µs)', xlim=(0, 20)); ax.axvspan(3, 12, color='#999999', alpha=.10); finish(ax)
    save(fig, 'supply-cycle')
    print('Generated five SVG/PNG figures from compact data; no SPICE.')


if __name__ == '__main__': main()
