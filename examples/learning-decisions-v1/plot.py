"""MIT: two ordinary graphs from the finite, known supply/loss observations."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analyze import calculate, CANDIDATES

matplotlib.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.titlesize': 12,
    'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'svg.hashsalt': 'learning-supply-loss-v1',
    'savefig.dpi': 160,
})
BLUE = '#126779'; ORANGE = '#bd5606'


def save(fig, output, name):
    fig.savefig(output / (name + '.svg'), bbox_inches='tight', metadata={'Date': None})
    fig.savefig(output / (name + '.png'), bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--feedback-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    result = calculate(a.feedback_zip)
    paths = result['paths']; x = np.arange(len(paths))
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7.1), layout='constrained')
    direct = np.array([r['direct_supply_injection_a_per_v'] * 1e6 for r in paths])
    reference = np.array([r['reference_supply_injection_a_per_v'] * 1e6 for r in paths])
    axes[0].bar(x, direct, color=ORANGE, width=.55, label='Direct drain/resistor path')
    axes[0].bar(x, reference, bottom=direct, color=BLUE, width=.55, label='Finite-reference gate path')
    axes[0].set(ylabel='Supply-induced current (µA / V)', ylim=(0, 195),
        title='Output feedback and reference error change different factors')
    for i, r in enumerate(paths):
        axes[0].text(i, direct[i] + reference[i] + 4,
            f"{direct[i] + reference[i]:.1f}\n$Z_O$ ≈ {r['injection_ohm_estimate'] / 1e3:.1f} kΩ",
            ha='center', va='bottom', fontsize=9)
    axes[0].legend(loc='upper right', fontsize=9, frameon=False)
    d = np.array([r['direct_supply_path_v_per_v'] for r in paths])
    g = np.array([r['reference_supply_path_v_per_v'] for r in paths])
    axes[1].bar(x, d, color=ORANGE, width=.55)
    axes[1].bar(x, g, bottom=d, color=BLUE, width=.55)
    axes[1].plot(x, [r['supply_v_per_v_observed'] for r in paths], 'k_',
        markersize=18, markeredgewidth=2, label='Original measured supply slope')
    axes[1].set(ylabel='Output / supply change (V / V)', ylim=(0, 2.45),
        title='Multiply each supply-induced current by its output resistance')
    axes[1].legend(loc='upper right', fontsize=9, frameon=False)
    for i, r in enumerate(paths):
        axes[1].text(i, d[i] + g[i] + .07, f"{r['supply_v_per_v_estimate']:.3f}", ha='center', fontsize=10)
    for ax in axes:
        ax.set_xticks(x, CANDIDATES); ax.grid(axis='y', alpha=.18); ax.set_axisbelow(True)
    save(fig, a.output, 'supply-reference-paths')

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.8), sharex=True, layout='constrained')
    for ax, condition, title in zip(axes, ['nominal', 'supply_097'], ['1.00-V supply', '0.97-V supply']):
        ax.axvspan(-25, 25, color='#c8ded5', alpha=.65, label='Original ±25-mV center band')
        for c, color in [('F75S', ORANGE), ('F75L', BLUE)]:
            errors = np.sort([r['error_v'] * 1e3 for r in result['errors']
                if r['candidate'] == c and r['condition'] == condition])
            ax.step(errors, np.arange(1, len(errors) + 1) / len(errors), where='post',
                color=color, label=c, linewidth=1.6)
            summary = next(r for r in result['loss'] if r['candidate'] == c and r['condition'] == condition)
            ax.axvline(summary['mean_error_v'] * 1e3, color=color, linestyle=':', alpha=.85, linewidth=1)
        ax.set(title=title + ' · same 521 known coordinates', ylim=(0, 1.03),
            ylabel='Fraction at or below error')
        ax.grid(alpha=.2)
    axes[0].legend(loc='upper left', fontsize=9, frameon=True)
    axes[1].set(xlabel='Signed output error from fixed 0.5-V target (mV)', xlim=(-150, 115))
    axes[1].text(.97, .06, 'Dotted lines: sample means\nNo recentering or normal-distribution fit',
        transform=axes[1].transAxes, ha='right', va='bottom', fontsize=9)
    save(fig, a.output, 'fixed-target-loss')
    print('PASS: two supply/loss SVG and PNG graphs from known public inputs')


if __name__ == '__main__':
    main()
