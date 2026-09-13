"""MIT: prior capacitor choice and the two later initial-D1 observations."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from load_budget import ROOT, read, calculate
from measurement_v1 import columns

matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
    'axes.titlesize': 12, 'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'svg.hashsalt': 'learning-initial-d1-load-budget-v1',
    'savefig.dpi': 160})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    result = calculate(); plan = read(ROOT / 'load-budget-plan.json')
    cfg = read(ROOT / 'data/interface-inputs.json')
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7), layout='constrained')
    cap = np.linspace(5, 25, 101)
    coefficient = max(plan['prospective_estimates']['seconds_per_farad'].values())
    axes[0].plot(cap, coefficient * cap * 1e-6, '--', color='#666666',
        label='Prior time/capacitance estimate')
    training = list(plan['known_evidence']['training_cases'].values())
    axes[0].plot([x['cload_f'] * 1e12 for x in training],
        [max(x['entry_upper_s'].values()) * 1e6 for x in training], 'x',
        color='#666666', markersize=8, label='Known 5/25-pF upper brackets')
    axes[0].axhline(.5, color='black', linestyle=':', linewidth=1.2,
        label='Illustrative 0.5-µs deadline')
    for name, color in [('D1-budget-12p', '#126779'), ('D1-budget-15p', '#bd5606')]:
        measured = result['cases'][name]['measured']; info = cfg['cases'][name]
        c = info['case']['cload_f'] * 1e12
        entries = measured['transient']['entries']
        worst = max(entries.values(), key=lambda x: x['bracket_s'][1])
        low, high = np.array(worst['bracket_s']) * 1e6
        axes[0].errorbar([c], [high], yerr=np.array([[high - low], [0.]]),
            fmt='o', color=color, capsize=5, markersize=6, label=f'New {c:g}-pF observation')
        with np.load(ROOT / 'data' / info['array_file'], allow_pickle=False) as z:
            arrays = {k: z[k] for k in z.files}
        t, v = columns(info, arrays, 4)
        for edge, linestyle in [('rise', '-'), ('fall', '--')]:
            entry = entries[edge]; x = (t - entry['origin_s']) * 1e6
            y = abs(v['v(out)'] - entry['target_v']) / entry['band_v']
            selected = (x >= .25) & (x <= .7)
            axes[1].semilogy(x[selected], y[selected], linestyle=linestyle,
                color=color, linewidth=1.6, label=f'{c:g} pF {edge}')
    axes[0].set(xlabel='Fixed external capacitance (pF)', ylabel='Worst-direction entry time (µs)',
        title='Use known loads to choose two later conditions', xlim=(4, 26), ylim=(0, 1.12))
    axes[0].legend(loc='upper left', fontsize=8.5, frameon=True)
    axes[1].axhline(1, color='#555555', linestyle=':', linewidth=1)
    axes[1].axvline(.5, color='black', linestyle=':', linewidth=1.2)
    axes[1].set(xlabel='Time after corresponding input 50% crossing (µs)',
        ylabel='Remaining error / actual-step 1% band', xlim=(.25, .7), ylim=(.1, 12),
        title='12 pF enters before the deadline; 15 pF enters after it')
    axes[1].legend(loc='upper right', fontsize=9, ncol=2)
    for ax in axes:
        ax.grid(alpha=.18); ax.set_axisbelow(True)
    for extension in ['svg', 'png']:
        kwargs = {'metadata': {'Date': None}} if extension == 'svg' else {}
        fig.savefig(a.output / ('initial-d1-load-budget.' + extension), bbox_inches='tight', **kwargs)
    plt.close(fig)
    print('PASS: capacitor/deadline SVG and PNG from the frozen forecasts and two later observations')


if __name__ == '__main__':
    main()
