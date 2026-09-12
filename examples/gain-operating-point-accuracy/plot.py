"""Original MIT code: regenerate figures from compact CSV, without SPICE."""
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES=['R37','A37','R75','A75']
COLORS=dict(zip(NAMES,['#0072B2','#D55E00','#009E73','#CC79A7']))

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'ads-gain-v1'})
    def read(n):
        with (a.data/(n+'.csv')).open() as f:return list(csv.DictReader(f))
    def save(fig,n):
        panels={'nominal-response':['gain','noise'],'error-paths':['current','supply','budget'],'confirmation':['nominal','supply']}
        if n in panels:
            fig.canvas.draw()
            for ax,label in zip(fig.axes,panels[n]):
                bbox=ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).expanded(1.04,1.06)
                fig.savefig(a.output/(n+'-'+label+'.svg'),metadata={'Date':None},bbox_inches=bbox)
        fig.savefig(a.output/(n+'.svg'),metadata={'Date':None},bbox_inches='tight')
        fig.savefig(a.output/(n+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
    ac=read('ac');noise=read('noise');fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    for n in NAMES:
        r=[x for x in ac if x['candidate']==n];f=np.array([float(x['frequency_hz']) for x in r]);h=np.array([complex(float(x['gain_real']),float(x['gain_imag'])) for x in r])
        axes[0].semilogx(f,abs(h),label=n,color=COLORS[n],lw=2)
        r=[x for x in noise if x['candidate']==n];f=np.array([float(x['frequency_hz']) for x in r]);v=np.sqrt([float(x['input_psd_v2_per_hz']) for x in r])*1e9
        axes[1].loglog(f,v,label=n,color=COLORS[n],lw=2)
    axes[0].set(xlabel='Frequency (Hz)',ylabel='External voltage gain magnitude (V/V)',title='Nominal signal transfer',ylim=(0,6.8))
    axes[1].set(xlabel='Frequency (Hz)',ylabel='Input noise density (nV / sqrt(Hz))',title='Nominal noise; external input reference')
    for ax in axes:ax.grid(alpha=.2);ax.legend(ncol=2,frameon=False)
    save(fig,'nominal-response')
    paths=read('sensitivities');budget=read('source-budget');fig,axes=plt.subplots(1,3,figsize=(12,3.9),layout='constrained');x=np.arange(4)
    z=[float(next(r for r in paths if r['candidate']==n)['injection_ohm'])/1000 for n in NAMES]
    supply=[float(next(r for r in paths if r['candidate']==n)['supply_v_per_v']) for n in NAMES]
    axes[0].bar(x,z,color=[COLORS[n] for n in NAMES]);axes[1].bar(x,supply,color=[COLORS[n] for n in NAMES])
    bottom=np.zeros(4)
    for role,color in [('input','#56B4E9'),('reference','#E69F00'),('load','#999999')]:
        v=np.array([sum(float(r['variance_mV2']) for r in budget if r['candidate']==n and r['role']==role) for n in NAMES])
        axes[2].bar(x,v,bottom=bottom,label=role,color=color);bottom+=v
    axes[0].set(ylabel='dVout / dIinjected (kohm)',title='Current-error transmission')
    axes[1].set(ylabel='dVout / dVDD (V/V)',title='Supply transmission')
    axes[2].set(ylabel='Output variance (mV squared)',title='First-order Local source budget');axes[2].legend(frameon=False)
    for ax in axes:ax.set_xticks(x,NAMES);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    save(fig,'error-paths')
    rows=read('confirmation');extent=max(300,25*np.ceil(max(abs(float(r['error_v']))*1000 for r in rows if r['error_v'])/25));fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,c,title in zip(axes,['nominal','supply_097'],['Nominal, 1.00 V supply','Transfer, 0.97 V supply']):
        for n in NAMES:
            vals=sorted(abs(float(r['error_v']))*1000 for r in rows if r['candidate']==n and r['condition']==c and r['numerical_usable']=='True' and r['terminal_supported']=='True')
            ax.step(vals,np.arange(1,len(vals)+1)/len(vals),where='post',color=COLORS[n],label=f'{n} (n={len(vals)})',lw=1.8)
        ax.axvline(25,color='#333333',ls='--',lw=1,label='25 mV center limit');ax.set(xlabel='Absolute error from fixed 0.5 V (mV)',ylabel='Cumulative fraction of supported observations',title=title,ylim=(0,1.02),xlim=(0,extent));ax.grid(alpha=.2);ax.legend(fontsize=9,frameon=False,loc='lower right')
    save(fig,'confirmation')
    contrasts=read('contrasts');r=[x for x in contrasts if x['condition']=='nominal'];fig,ax=plt.subplots(figsize=(5.4,3.8),layout='constrained')
    m=np.array([float(x['mean_mae_reduction_v']) for x in r])*1000;lo=np.array([float(x['ci_low_v']) for x in r])*1000;hi=np.array([float(x['ci_high_v']) for x in r])*1000
    ax.errorbar(m,np.arange(len(r)),xerr=[m-lo,hi-m],fmt='o',color='#0072B2',capsize=4)
    ax.axvline(0,color='#777777',lw=1);ax.set_yticks(np.arange(len(r)),[x['left']+' to '+x['right'] for x in r]);ax.invert_yaxis()
    ax.set(xlabel='Reduction in mean absolute center error (mV)',title='Nominal paired comparisons\nFamily 95% t intervals');ax.grid(axis='x',alpha=.2)
    save(fig,'paired-comparisons')
    print('Generated combined figures and eight responsive SVG panels from compact CSV; no simulator invoked.')

if __name__=='__main__':main()
