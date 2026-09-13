"""MIT: current/headroom and local-response comparison from selected observations."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analyze import ROOT, bandwidth, columns, inputs, read, sha


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    plan, prior = inputs(a.source_zip)
    summary = read(ROOT/'data/summary.json')
    index = read(ROOT/'data/index.json')
    assert summary['plan_sha256'] == index['plan_sha256'] == sha(ROOT/'plan.json')
    new = summary['cases']['C-high-common-mode']
    old = read(ROOT/'data/prior.json')['cases']['C-source-0p70']
    plt.rcParams.update({'font.size': 10, 'svg.hashsalt': 'learning.initial-d1-common-mode.v1',
                         'font.family': 'DejaVu Sans', 'axes.spines.top': False, 'axes.spines.right': False})
    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.6), layout='constrained')
    x = np.arange(4)
    olds = np.array([old['tail_delivered_a'], old['source_branch_a'], old['sink_branch_a'], old['net_after_load_a']])*1e6
    news = np.array([new['tail_delivered_a'], new['source_branch_a'], new['sink_branch_a'], new['net_after_resistor_a']])*1e6
    for offset, values, label, color in [(-.19, olds, 'Known: CM 0.3 V', '#196b80'), (.19, news, 'New: CM 0.7 V', '#bb580b')]:
        bars = left.bar(x+offset, values, width=.36, label=label, color=color)
        left.bar_label(bars, fmt='%.3f', fontsize=8, padding=3)
    left.set_xticks(x, ['Tail P5', 'P2 source', 'N4 sink', 'Net after\n1 MΩ'])
    left.set_ylabel('Current magnitude (µA)')
    left.set_ylim(0, 18)
    left.set_title('C: fixed output 0.7 V and input difference +0.1 V')
    left.legend(loc='upper right', frameon=False, fontsize=9)
    left.grid(axis='y', alpha=.2)
    left.set_axisbelow(True)
    oldcase = prior['D1-5p']['case']
    f, v = columns(oldcase['analyses'][3]['recipe'], prior['D1-5p']['arrays']['analysis3'])
    h = v['v(out)']/v['v(source)']
    old_gain, old_bw = float(abs(h[0])), bandwidth(f, h)['hz']/1e6
    right.semilogx(f/1e6, abs(h), color='#196b80', label='Known: command 0.3 V')
    entry = index['cases']['D1-high-command']
    assert sha(ROOT/'data'/entry['file']) == entry['sha256']
    with np.load(ROOT/'data'/entry['file'], allow_pickle=False) as z:
        f, v = columns(plan['cases']['D1-high-command']['analyses'][1], z['analysis1'])
    h = v['v(out)']/v['v(source)']
    right.semilogx(f/1e6, abs(h), color='#bb580b', label='New: command 0.7 V')
    for bw, gain, color, label in [(old_bw, old_gain, '#196b80', f'{old_bw:.3f} MHz'),
                                  (summary['cases']['D1-high-command']['bandwidth']['hz']/1e6,
                                   abs(complex(*summary['cases']['D1-high-command']['gain_one_khz'])), '#bb580b', '0.863 MHz')]:
        y = gain*10**(-3/20)
        right.scatter([bw], [y], color=color, s=22)
        right.annotate(label, (bw, y), xytext=(9, 9), textcoords='offset points', color=color, fontsize=9)
    right.set_ylim(0, 1.12)
    right.set_xlim(.001, 1e3)
    right.set_xlabel('Frequency (MHz)')
    right.set_ylabel('Closed-loop signal gain |Vout/Vsource| (V/V)')
    right.set_title('D1: same six-unit core and 5-pF load')
    right.legend(loc='upper right', frameon=False, fontsize=9)
    right.grid(alpha=.2)
    a.output.mkdir(parents=True)
    fig.savefig(a.output/'common-mode-interface.svg', metadata={'Date': None})
    fig.savefig(a.output/'common-mode-interface.png', dpi=160)
    plt.close(fig)
    print('Saved selected common-mode/current and local AC comparison.')


if __name__ == '__main__':
    main()
