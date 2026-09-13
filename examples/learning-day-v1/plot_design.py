"""MIT: deterministic ordinary plots of the selected saved design observations."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from design import calculate


def main():
    p = argparse.ArgumentParser(); p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    assert not a.output.exists(), 'Use a new figure destination'
    data = calculate(); a.output.mkdir(parents=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'none',
                         'svg.hashsalt': 'adn-learning-design-v1', 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.grid': True, 'grid.alpha': .22})

    def save(fig, name):
        fig.savefig(a.output/(name+'.svg'), metadata={'Date': None}, bbox_inches='tight')
        fig.savefig(a.output/(name+'.png'), dpi=150, bbox_inches='tight')
        plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.6), sharex=True, constrained_layout=True)
    for gate, color, marker in [(.45, '#005b96', 'o'), (.50, '#ad4b00', 's'), (.55, '#006c58', '^')]:
        rows = [r for r in data['feasibility'] if abs(r['gate_v']-gate)<1e-12]
        x = [r['l_um'] for r in rows]
        axes[0].plot(x, [-r['optimistic_gain_32p5u'] for r in rows], color=color, marker=marker, label=f'Gate {gate:.2f} V')
        axes[1].plot(x, [r['estimated_total_width_32p5u_um'] for r in rows], color=color, marker=marker)
    axes[0].axhline(6, color='black', linestyle='--', linewidth=1, label='Gain magnitude required: 6')
    axes[0].set(ylabel='Optimistic gain magnitude (V/V)', ylim=(3.2, 6.7))
    axes[1].set(ylabel='Estimated total width (µm)', xlabel='Unit length (µm)')
    axes[0].legend(fontsize=9, ncol=2, loc='upper left', bbox_to_anchor=(0, -.03), framealpha=1)
    axes[0].set_ylim(2.5, 6.7)
    fig.suptitle('Current efficiency and width trade places\nGrounded source; width scaled to 32.5 µA, before finite Rs', fontsize=12)
    save(fig, 'sizing-estimate')

    m = data['mirror']; x = np.array(m['voltage_v']); err = 100*np.array(m['fractional_ratio_error'])
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    ax.axhspan(-1, 1, color='#dfede4', label='Illustrative ±1% copying tolerance')
    ax.axhline(0, color='gray', linewidth=.8)
    ax.plot(x, err, color='#005b96', marker='o', markevery=5, label='Same two saved physical units')
    ax.set(xlabel='Clamped output voltage (V)', ylabel='Current copying error (%)', xlim=(.29, .81), ylim=(-1.3, 3.0))
    ax.legend(fontsize=9, loc='upper left')
    ax.set_title('Equal gate voltage does not fix output current\n30-µA reference; SOURCE_TRANSFER_HYPOTHESIS realization', fontsize=11)
    save(fig, 'mirror-compliance')
    print('Regenerated sizing-estimate and mirror-compliance SVG/PNG from saved inputs.')


if __name__ == '__main__':
    main()
