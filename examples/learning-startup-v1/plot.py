"""MIT: one current/voltage figure from the bounded saved R75 calculation."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--analysis',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists()
    summary=json.loads((a.analysis/'summary.json').read_text())
    assert summary['schema']=='learning.startup-current.analysis.v1'
    with (a.analysis/'waveforms.csv').open() as f:rows=list(csv.DictReader(f))
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans',
        'svg.hashsalt':'learning.startup-current.v1','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(10.5,7.2),layout='constrained',sharex='col')
    for col,(speed,title,limit) in enumerate([('fast','10 ns supply ramp',.55),('slow','1 µs supply ramp',1.5)]):
        name='R75-input_first-'+speed+'-8p';r=[v for v in rows if v['case']==name]
        values={k:np.array([float(v[k]) for v in r]) for k in r[0] if k!='case'}
        x=values['relative_time_s']*1e6;top,bottom=axes[:,col]
        top.plot(x,values['output_v'],color='#196b80',lw=2,label='Observed R75 output')
        top.plot(x,values['passive_output_v'],color='#b76314',ls='--',label='RC screen, NMOS omitted')
        top.axhline(summary['reference']['final_output_v'],color='#555555',lw=.7,ls=':',label='Independent final OP')
        top.set_title(title);top.set_ylim(-.02,1.02);top.set_ylabel('Output voltage (V)')
        if col==0:top.legend(loc='lower right',frameon=False,fontsize=8)
        for key,label,color,style,width in [('supply_a','Supply delivery','#222222','-',2.),
            ('capacitor_a','External capacitor','#196b80','-',1.5),
            ('nmos_drain_a','Sum of NMOS drain terminals','#bb580b','-',1.5),
            ('load_a','100-kΩ load','#527b3c','--',1.2)]:
            bottom.plot(x,values[key]*1e6,label=label,color=color,ls=style,lw=width)
        peak=summary['cases'][name]['current_peak']
        bottom.scatter(peak['relative_to_ramp_start_s']*1e6,peak['supply_delivered_a']*1e6,color='#222222',s=24,zorder=5)
        bottom.annotate(f"{peak['supply_delivered_a']*1e6:.3f} µA",(peak['relative_to_ramp_start_s']*1e6,peak['supply_delivered_a']*1e6),
            xytext=(10,6) if col==0 else (-72,9),textcoords='offset points',fontsize=9)
        bottom.set_xlabel('Time from supply-ramp start (µs)');bottom.set_ylabel('Current (µA)')
        bottom.set_ylim(-4,159)
        if col==1:bottom.legend(loc='upper left',frameon=False,fontsize=8)
        for ax in [top,bottom]:
            ax.set_xlim(0,limit);ax.grid(alpha=.2);ax.axvline(summary['cases'][name]['ramp_time_s']*1e6,color='#888888',lw=.8,ls=':')
    fig.suptitle('Same R75, 0.45-V input already present, 8-pF output load',fontsize=12)
    a.output.mkdir(parents=True)
    fig.savefig(a.output/'startup-current-and-charge.svg',metadata={'Date':None})
    fig.savefig(a.output/'startup-current-and-charge.png',dpi=160);plt.close(fig)
    # The same figure needs readable labels on a narrow screen; stack its panels.
    plt.rcParams.update({'font.size':11})
    fig,axes=plt.subplots(4,1,figsize=(4.2,12.2),layout='constrained')
    for col,(speed,title,limit) in enumerate([('fast','10 ns supply ramp',.55),('slow','1 µs supply ramp',1.5)]):
        name='R75-input_first-'+speed+'-8p';r=[v for v in rows if v['case']==name]
        values={k:np.array([float(v[k]) for v in r]) for k in r[0] if k!='case'}
        x=values['relative_time_s']*1e6;top,bottom=axes[2*col:2*col+2]
        top.plot(x,values['output_v'],color='#196b80',lw=2,label='Observed R75 output')
        top.plot(x,values['passive_output_v'],color='#b76314',ls='--',label='RC: NMOS omitted')
        top.axhline(summary['reference']['final_output_v'],color='#555555',lw=.7,ls=':',label='Independent final OP')
        top.set_title(title);top.set_ylim(-.02,1.02);top.set_ylabel('Output voltage (V)')
        top.tick_params(labelbottom=False)
        if col==0:top.legend(loc='lower right',frameon=False,fontsize=10)
        for key,label,color,style,width in [('supply_a','Supply delivery','#222222','-',2.),
            ('capacitor_a','External capacitor','#196b80','-',1.5),
            ('nmos_drain_a','NMOS drain terminals','#bb580b','-',1.5),
            ('load_a','100-kΩ load','#527b3c','--',1.2)]:
            bottom.plot(x,values[key]*1e6,label=label,color=color,ls=style,lw=width)
        peak=summary['cases'][name]['current_peak']
        bottom.scatter(peak['relative_to_ramp_start_s']*1e6,peak['supply_delivered_a']*1e6,color='#222222',s=24,zorder=5)
        bottom.annotate(f"{peak['supply_delivered_a']*1e6:.3f} µA",(peak['relative_to_ramp_start_s']*1e6,peak['supply_delivered_a']*1e6),
            xytext=(10,5) if col==0 else (10,8),textcoords='offset points',fontsize=10)
        bottom.set_xlabel('Time from ramp start (µs)');bottom.set_ylabel('Current (µA)');bottom.set_ylim(-4,184)
        bottom.legend(loc='upper right' if col==0 else 'upper left',frameon=False,fontsize=10)
        for ax in [top,bottom]:
            ax.set_xlim(0,limit);ax.grid(alpha=.2);ax.axvline(summary['cases'][name]['ramp_time_s']*1e6,color='#888888',lw=.8,ls=':')
    fig.suptitle('R75: input already at 0.45 V\nSame six NMOS units and 8-pF load',fontsize=11)
    fig.savefig(a.output/'startup-current-and-charge-narrow.svg',metadata={'Date':None})
    fig.savefig(a.output/'startup-current-and-charge-narrow.png',dpi=160);plt.close(fig)
    print('Saved wide and narrow layouts of one voltage/current comparison; no measurement uses the drawing grid.')


if __name__=='__main__':main()
