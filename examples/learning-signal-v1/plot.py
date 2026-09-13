"""MIT: two figures from the finite saved-cycle explanation."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analyze import calculate


def main():
    p = argparse.ArgumentParser(); p.add_argument('--feedback-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    result, harmonics, waves = calculate(a.feedback_zip)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
        'axes.labelsize':10,'legend.fontsize':9,'svg.fonttype':'none',
        'svg.hashsalt':'learning.nonlinearity.v1','axes.spines.top':False,
        'axes.spines.right':False})
    colors = ['#176c80','#c65b0b']
    rows = {r['case']:r for r in result['cases']}
    estimates = {r['candidate']:r for r in result['retrospective_extrapolations']}

    def save(fig, name):
        fig.savefig(a.output/(name+'.svg'),metadata={'Date':None})
        fig.savefig(a.output/(name+'.png'),dpi=160,metadata={'Software':'Analog Design Studies'})
        plt.close(fig)

    fig, axes = plt.subplots(1,2,figsize=(8.6,3.9))
    estimate = estimates['R75']; slope = estimate['linear_slope_from_old_ac']
    for key, color in zip(['R75-10m','R75-20m'],colors):
        w = waves[key]; r = rows[key]; x = w['input_v']-.45
        error = w['output_v']-r['independent_op_output_v']-slope*x
        axes[0].plot(x*1e3,error*1e3,color=color,lw=1.5,label=f"Known {r['amplitude_peak_v']*1000:.0f}-mV cycle")
    x = np.linspace(-.02,.02,201)
    axes[0].plot(x*1e3,.5*estimate['static_second_derivative_per_v']*x*x*1e3,
        '--',color='#222222',lw=1,label='Saved local DC curvature estimate')
    axes[0].set(xlabel='External input excursion (mV)',ylabel='Output minus OP and local linear term (mV)',
        title='R75: finite transfer bends below its tangent')
    axes[0].legend(loc='lower center',fontsize=8)
    orders = np.arange(2,6)
    for j,(key,color) in enumerate(zip(['R75-10m','R75-20m'],colors)):
        table = {r['harmonic']:r for r in harmonics if r['case']==key}
        axes[1].bar(orders+(j-.5)*.32,[100*table[n]['fraction_of_fundamental'] for n in orders],
            width=.3,color=color,label=f'{rows[key]["amplitude_peak_v"]*1000:.0f} mV input')
    axes[1].set(xticks=orders,xlabel='Output harmonic order',ylabel='Peak / fundamental (%)',
        title='Second harmonic dominates the distortion',yscale='log',ylim=(1e-5,10))
    axes[1].legend(loc='upper right')
    for ax in axes:ax.grid(alpha=.18);ax.set_axisbelow(True)
    fig.subplots_adjust(left=.095,right=.985,top=.88,bottom=.16,wspace=.36)
    save(fig,'local-gain-and-curvature')

    fig,axes = plt.subplots(2,1,figsize=(7.1,6.0))
    names = ['R75','F75L']; x = np.arange(4)
    for j,(name,color) in enumerate(zip(names,colors)):
        r = rows[name+'-20m']; heights = [r[k]*1e3 for k in ['observed_input_fundamental_peak_v',
            'gate_fundamental_peak_v','source_fundamental_peak_v','vgs_fundamental_peak_v']]
        axes[0].bar(x+(j-.5)*.34,heights,width=.32,color=color,label=name)
        for xx,y in zip(x+(j-.5)*.34,heights):axes[0].text(xx,y+.35,f'{y:.2f}',ha='center',va='bottom',fontsize=8)
        w = waves[name+'-20m']
        axes[1].plot((w['time_s']-.0009)*1e6,w['input_delivered_current_a']*1e6,color=color,label=name,lw=1.6)
    axes[0].set(xticks=x,xticklabels=['External input','MOS gate','MOS source','Gate − source'],
        ylabel='Fundamental peak (mV)',title='Same 20-mV source; different internal voltage motion',ylim=(0,25))
    axes[0].legend(loc='upper right',ncol=2)
    axes[1].axhline(0,color='#444444',lw=.8)
    axes[1].set(xlabel='Time within the saved final cycle (µs)',ylabel='Input-source delivered current (µA)',
        title='Positive supplies current; negative absorbs it',xlim=(0,100),ylim=(-1.9,1.1))
    axes[1].legend(loc='lower left')
    for ax in axes:ax.grid(alpha=.18);ax.set_axisbelow(True)
    fig.subplots_adjust(left=.13,right=.97,top=.92,bottom=.1,hspace=.46)
    save(fig,'feedback-gate-motion')
    print('PASS: two figures from six known nominal cycles')


if __name__ == '__main__':
    main()
