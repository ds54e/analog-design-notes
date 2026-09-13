"""MIT: ordinary reproducible figures for the unchanged initial-D1 interface."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from interface import inputs,columns,calculate


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists(), 'Keep prior plots; use a new destination'
    cfg,plan,arrays=inputs();summary=calculate();a.output.mkdir(parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
        'svg.hashsalt':'adn-initial-d1-interface-v1','axes.spines.top':False,'axes.spines.right':False,
        'axes.grid':True,'grid.alpha':.22})
    colors={'D1-5p':'#005b96','D1-25p':'#ad4b00'}
    labels={'D1-5p':'5 pF','D1-25p':'25 pF'}
    def save(fig,name):
        fig.savefig(a.output/(name+'.svg'),metadata={'Date':None},bbox_inches='tight')
        fig.savefig(a.output/(name+'.png'),dpi=150,bbox_inches='tight');plt.close(fig)

    fig,axs=plt.subplots(2,1,figsize=(6.4,6.6),constrained_layout=True)
    for name,color in colors.items():
        f,v=columns(cfg['cases'][name],arrays[name],3);h=v['v(out)']/v['v(source)']
        axs[0].semilogx(f,abs(h),color=color,label=labels[name])
        t,v=columns(cfg['cases'][name],arrays[name],4);mask=(t>=.48e-6)&(t<=1.75e-6)
        axs[1].plot((t[mask]-.50625e-6)*1e6,v['v(out)'][mask],color=color,label=labels[name])
    axs[0].set(xlim=(1e3,1e8),ylim=(0,1.05),xlabel='Frequency (Hz)',ylabel='|Vout / Vsource| (V/V)')
    axs[0].legend();axs[0].set_title('Same initial six-unit core; changed capacitor')
    target=summary['cases']['D1-25p']['dc'][2]['output_v'];band=summary['cases']['D1-25p']['transient']['step_band_v']
    axs[1].axhspan(target-band,target+band,color='#dce9dd',label='1% band around independent DC endpoint')
    axs[1].axhline(.35,color='gray',linestyle=':',label='Final input source: 0.35 V')
    axs[1].set(xlim=(-.025,1.2),xlabel='Time from rising input 50% crossing (µs)',ylabel='Output voltage (V)')
    axs[1].legend(fontsize=8,loc='lower right')
    save(fig,'initial-d1-load-response')

    fig,axs=plt.subplots(2,1,figsize=(6.4,6.6),constrained_layout=True)
    for name,color in colors.items():
        t,v=columns(cfg['cases'][name],arrays[name],4);m=summary['cases'][name];step=m['dc'][2]['output_v']-m['dc'][0]['output_v']
        net=-v['i(VtM2d)']-v['i(VtM4d)']-v['v(out)']/1e6+v['i(Vfeedback)']
        for edge,origin,target,style in [('rise',.50625e-6,m['dc'][2]['output_v'],'-'),('fall',4.50625e-6,m['dc'][0]['output_v'],'--')]:
            mask=(t>=origin)&(t<=origin+1.2e-6)
            x=(t[mask]-origin)*1e6
            axs[0].semilogy(x,100*abs(v['v(out)'][mask]-target)/abs(step),color=color,linestyle=style,label=labels[name]+' '+edge)
            axs[1].plot(x,net[mask]*1e6,color=color,linestyle=style,label=labels[name]+' '+edge)
    axs[0].axhline(1,color='black',linewidth=1,linestyle=':',label='1% of actual DC output step')
    axs[0].set(ylim=(.01,110),ylabel='Remaining endpoint error (% of step)')
    axs[0].legend(fontsize=8,ncol=2)
    axs[1].axhline(0,color='gray',linewidth=.8)
    axs[1].set(ylabel='Net capacitor current (µA)')
    for ax in axs:ax.set(xlim=(0,1.2),xlabel='Time from corresponding input 50% crossing (µs)')
    axs[0].set_title('Charging current shrinks as feedback reduces error')
    save(fig,'initial-d1-step-current')

    fig,ax=plt.subplots(figsize=(6.4,4.5),constrained_layout=True)
    ax.axvspan(.25,.35,color='#eeeeee',label='Functional example output neighborhood')
    for direction,color,marker in [('source','#005b96','o'),('sink','#ad4b00','s')]:
        names=[n for n in cfg['cases'] if n.startswith('C-'+direction)]
        points=sorted((summary['cases'][n]['op']['v(out)'],summary['cases'][n]['available_current_a']*1e6) for n in names)
        xy=np.array(points);ax.plot(xy[:,0],xy[:,1],marker=marker,color=color,linestyle='--',
                                   label='Full input difference '+('+' if direction=='source' else '−')+'0.10 V')
    ax.axhline(0,color='gray',linewidth=.8)
    ax.set(xlabel='Clamped output voltage (V)',ylabel='OTA available current before Rload (µA)',
           xlim=(0,.75),ylim=(-12,12),title='Headroom changes the current supplied to a port\nC fixture: fixed 0.30-V common mode; six discrete OP clamps')
    ax.legend(fontsize=8,loc='center right')
    save(fig,'initial-d1-clamped-current')
    print('Regenerated three initial-D1 interface SVG/PNG figures from selected observations.')


if __name__=='__main__':main()
