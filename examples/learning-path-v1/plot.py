"""MIT: deterministic ordinary plots from the selected learning inputs only."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from analyze import inputs, column, calculate

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,
    'axes.titlesize':12,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,
    'svg.hashsalt':'adn-learning-v1','axes.spines.top':False,'axes.spines.right':False})
COLORS={'R37':'#176b83','R75':'#ad4b16','left':'#176b83','right':'#ad4b16'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists(), 'Use a new figure directory'
    a.output.mkdir(parents=True)
    here=Path(__file__).resolve().parent;cfg,arrays=inputs(here);summary=calculate(here)
    c=lambda k,i,n:column(cfg,arrays,k,i,n)
    axis=lambda k,i:arrays[k]['analysis'+str(i)][:,0]
    def save(fig,name):
        fig.tight_layout(pad=1.1)
        fig.savefig(a.output/(name+'.svg'),metadata={'Date':None})
        fig.savefig(a.output/(name+'.png'),dpi=160)
        plt.close(fig)
    def read(name):
        with (here/'data'/name).open() as f:return list(csv.DictReader(f))
    fig,ax=plt.subplots(figsize=(5.4,3.5))
    rows=read('gain-ac.csv')
    for key in ['R37','R75']:
        r=[x for x in rows if x['candidate']==key];f=np.array([float(x['frequency_hz']) for x in r]);h=np.array([complex(float(x['gain_real']),float(x['gain_imag'])) for x in r])
        ax.semilogx(f/1e6,abs(h),label=key,color=COLORS[key],lw=2)
    ax.axhline(6*10**(-3/20),color='#777777',ls=':',lw=1,label='−3 dB from 6')
    ax.set(xlabel='Frequency (MHz)',ylabel='Voltage gain magnitude (V/V)',xlim=(.001,100),ylim=(0,6.5))
    ax.legend(loc='lower left');ax.grid(alpha=.2);save(fig,'resistor-gain')
    fig,ax=plt.subplots(figsize=(5.4,3.5));rows=read('gain-noise.csv')
    for key in ['R37','R75']:
        r=[x for x in rows if x['candidate']==key]
        ax.loglog([float(x['frequency_hz'])/1e3 for x in r],[np.sqrt(float(x['input_psd_v2_per_hz']))*1e9 for x in r],label=key,color=COLORS[key],lw=2)
    ax.set(xlabel='Frequency (kHz)',ylabel='Input noise ASD (nV/√Hz)',xlim=(1,1000));ax.legend();ax.grid(alpha=.2);save(fig,'resistor-noise')
    fig,ax=plt.subplots(figsize=(5.4,3.6));x=axis('A',1)*1e3
    ax.plot(x,c('A',1,'i(VtM1d)')*1e6,label='Left drain',color=COLORS['left'],lw=2)
    ax.plot(x,c('A',1,'i(VtM2d)')*1e6,label='Right drain',color=COLORS['right'],lw=2)
    ax.set(xlabel='Differential input (mV)',ylabel='Drain current (µA)',ylim=(0,10));ax.legend();ax.grid(alpha=.2);save(fig,'differential-steering')
    fig,ax=plt.subplots(figsize=(5.4,3.5));x=axis('B',1)
    ax.plot(x,c('B',1,'i(VtMtd)')*1e6,color=COLORS['left'],lw=2,label='Finite tail B')
    ax.axhline(10,color='#777777',ls=':',label='Nominal target 10 µA')
    ax.axvline(.6,color='#777777',ls='--',lw=1)
    ax.set(xlabel='Input common mode (V)',ylabel='Tail drain current (µA)',ylim=(4.5,11));ax.legend(loc='lower right');ax.grid(alpha=.2);save(fig,'finite-tail')
    fig,ax=plt.subplots(figsize=(5.4,3.7));x=axis('C',1)*1e3
    available=-c('C',1,'i(VtM2d)')-c('C',1,'i(VtM4d)')
    ax.plot(x,available*1e6,color=COLORS['left'],lw=2,label='OTA available current')
    ax.axhline(0,color='#777777',ls=':',label='Unloaded requirement')
    ax.axhline(.3,color=COLORS['right'],ls='--',label='1 MΩ at 0.3 V')
    for n,yy,color in [('unloaded_nulls_V',0,COLORS['left']),('loaded_nulls_V',.3,COLORS['right'])]:
        xx=summary['cases']['C'][n][0]*1e3;ax.plot(xx,yy,'o',color=color);ax.annotate(f'{xx:+.3f} mV',(xx,yy),xytext=(0,12 if yy else -22),textcoords='offset points',ha='center',fontsize=10)
    ax.set(xlabel='Differential input (mV)',ylabel='Current at clamped output (µA)',xlim=(-4,4),ylim=(-.5,.8));ax.legend(loc='upper left');ax.grid(alpha=.2);save(fig,'ota-current-balance')
    fig,ax=plt.subplots(figsize=(5.4,3.5));x=axis('D1',1)[1:-1];out=c('D1',1,'v(out)')[1:-1]
    ax.plot(x,(out-x)*1e3,lw=2,color=COLORS['left']);ax.axhline(0,color='#777777',ls=':')
    ax.set(xlabel='External input source (V)',ylabel='Output − input (mV)',title='Initial D1 · 1 MΩ ∥ 5 pF');ax.grid(alpha=.2);save(fig,'buffer-accuracy')
    fig,ax=plt.subplots(figsize=(5.4,3.5));f=axis('D1',3);h=c('D1',3,'v(out)')/c('D1',3,'v(source)')
    ax.semilogx(f/1e6,20*np.log10(abs(h)),lw=2,color=COLORS['left'])
    ax.axhline(20*np.log10(abs(h[0]))-3,color='#777777',ls=':',label='−3 dB from 1-kHz gain')
    ax.set(xlabel='Frequency (MHz)',ylabel='Voltage gain (dB)',xlim=(.001,100),ylim=(-30,1));ax.legend(loc='lower left');ax.grid(alpha=.2);save(fig,'buffer-response')
    fig,axes=plt.subplots(2,1,figsize=(5.4,6.2),sharex=True)
    t=axis('D1',12);y=c('D1',12,'v(out)')
    steps=summary['cases']['D1']['dc_endpoint_steps']
    for ax,(name,origin) in zip(axes,[('rise',.50625e-6),('fall',4.50625e-6)]):
        x=(t-origin)*1e9;use=(x>=-15)&(x<=300);target=steps['directions'][name]['target_V']
        ax.plot(x[use],(y[use]-target)*1e3,color=COLORS['left'],lw=2,label='Output − DC endpoint')
        ax.axhspan(-steps['band_V']*1e3,steps['band_V']*1e3,color='#d6e8c4',label='1% of actual DC step')
        ax.set(title='Input '+name,ylabel='Endpoint error (mV)',ylim=(-105,5) if name=='rise' else (-5,105));ax.grid(alpha=.2)
    axes[0].legend(loc='lower right');axes[1].set_xlabel('Time from input 50% crossing (ns)');save(fig,'buffer-step')
    print('Regenerated 8 SVG figures and PNG inspection copies from selected public inputs.')


if __name__=='__main__':main()
