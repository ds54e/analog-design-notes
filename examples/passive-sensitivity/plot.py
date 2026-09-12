"""MIT: ordinary plots from the compact passive-study data; no SPICE."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAMES=['R37','R75','A37','A75']
COLORS=dict(zip(NAMES,['#0072B2','#009E73','#D55E00','#CC79A7']))


def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'ads-passive-v1'})
    def read(name):
        with (a.data/(name+'.csv')).open() as f:return list(csv.DictReader(f))
    def save(fig,name):
        fig.savefig(a.output/(name+'.svg'),metadata={'Date':None},bbox_inches='tight')
        fig.savefig(a.output/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
    slopes=read('sensitivities');s={n:{r['parameter']:float(r['observed_v_per_logR']) for r in slopes if r['candidate']==n} for n in NAMES}
    fig,ax=plt.subplots(figsize=(6.2,4.5),layout='constrained');x=np.arange(4)
    ax.barh(x-.18,[10*s[n]['rs'] for n in NAMES],height=.28,label='Rs: weakens NMOS sink',color='#0072B2')
    ax.barh(x+.15,[10*s[n].get('rd',s[n].get('rbias')) for n in NAMES],height=.28,label='Rd / Rbias: weakens supplied current',color='#D55E00')
    ax.axvline(0,color='#444444',lw=.8);ax.set_yticks(x,NAMES);ax.invert_yaxis()
    ax.set(xlabel='Output change for 0.01 increase in ln(R/R0) (mV)',title='Opposite internal-resistor error paths')
    ax.legend(frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(0,-.22));ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    save(fig,'internal-sensitivity')
    boundary=read('nominal-boundaries')
    for mode,label in [('common_internal','Both internal resistors change together'),('ratio_internal','Rs / RX changes; geometric mean fixed')]:
        fig,ax=plt.subplots(figsize=(6.0,3.9),layout='constrained')
        for n in NAMES:
            r=sorted([r for r in boundary if r['candidate']==n and r['mode']==mode],key=lambda r:float(r['coordinate_fraction']))
            xv=[100*float(t['coordinate_fraction']) for t in r];y=[1e3*(float(t['out_v'])-.5) for t in r]
            ax.plot(xv,y,'o',color=COLORS[n],label=n,lw=1.8,ms=4,
                    linestyle='--' if n=='A75' else '-',zorder=3 if n=='A75' else 2)
        ax.axhspan(-25,25,color='#999999',alpha=.12,label='±25 mV center band');ax.axhline(0,color='#555555',lw=.7)
        ax.set(xlabel='Multiplicative resistance / ratio change (%)',ylabel='Output error about fixed 0.5 V (mV)',title=label)
        ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9,ncol=3)
        save(fig,'nominal-'+mode.split('_')[0])
    stats=json.loads((a.data/'statistics.json').read_text())
    for kind in ['common','ratio']:
        fig,ax=plt.subplots(figsize=(6.0,4.1),layout='constrained');z=np.linspace(0,.12,121)
        for n in ['R37','R75']:
            b=stats['baseline'][n]['mse_v2'];coefficient=stats['slope_moments'][n]['nominal_'+kind+'_v']
            ax.plot(z*100,np.sqrt(b+(coefficient*z)**2)*1000,color=COLORS[n],label=n+' nominal-gradient estimate')
            group=[g for g in stats['groups'] if g['candidate']==n and g['kind']==kind]
            ax.scatter([100*g['log_amplitude'] for g in group],[1000*g['rms_v'] for g in group],color=COLORS[n],marker='s',s=46,label=n+' native two-sign result',zorder=3)
        if kind=='ratio':
            cross=stats['linear_ratio_break_even_log_amplitude']*100
            ax.axvline(cross,color='#555555',ls=':',label=f'Linear crossing: {cross:.2f}')
        ax.set(xlabel='Two-point log-error amplitude × 100 (not percent R)',ylabel='RMS output error about fixed 0.5 V (mV)',
               title='Known MOS cohort + synthetic '+kind+' resistance error')
        ax.grid(alpha=.2);ax.legend(fontsize=8.5,frameon=False)
        save(fig,'finite-'+kind)
    fig,ax=plt.subplots(figsize=(6.0,4.0),layout='constrained')
    for n in ['R37','R75']:
        g=[g for g in stats['groups'] if g['candidate']==n]
        xpos=np.arange(3)+(-.16 if n=='R37' else .16)
        ax.bar(xpos,[100*g0['dc_ac_current_pass']/g0['observations'] for g0 in g],width=.3,color=COLORS[n],label=n)
    ax.set_xticks(np.arange(3),['Common ±0.10','Ratio ±0.06','Ratio ±0.10'])
    ax.set(ylabel='Assessed DC/AC/current pass fraction (%)',ylim=(0,100),title='Accuracy ranking does not replace function checks')
    ax.set_xlabel('Synthetic equal-weight log-resistance conditions')
    ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True);ax.legend(frameon=False)
    save(fig,'finite-function')
    print('Generated six SVG/PNG figures from compact CSV/statistics; no SPICE.')


if __name__=='__main__':main()
