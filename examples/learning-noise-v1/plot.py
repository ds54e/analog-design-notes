"""MIT: one figure for coherent resistor paths and an explicit partial noise budget."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from analyze import calculate


def main():
    p=argparse.ArgumentParser();p.add_argument('--feedback-zip',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);result,tables=calculate(a.feedback_zip)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'legend.fontsize':9,'svg.fonttype':'none','svg.hashsalt':'learning.resistor-noise-path.v1','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(9.6,4.5));r=result['cases'][1]['one_khz']
    values=[r['gate_path_real_ohm'],r['opposite_output_path_real_ohm'],r['combined_branch_real_ohm']]
    labels=['Through gate\nRt × H','Opposite output\n−Zout','Sum for one source\nRt × H − Zout'];colors=['#176c80','#c65b0b','#40566d']
    axes[0].bar(np.arange(3),np.array(values)/1000,color=colors,width=.58)
    for x,y in enumerate(values):axes[0].text(x,y/1000-2,f'{y/1000:.2f}',ha='center',va='top',fontsize=9)
    axes[0].set(xticks=np.arange(3),xticklabels=labels,ylabel='Real output transfer per branch current (kΩ)',title='F75L at 1 kHz: the two effects add',ylim=(-94,5))
    axes[0].axhline(0,color='#555',lw=.7)
    axes[0].text(.03,.06,'Small imaginary part retained in the calculation.',transform=axes[0].transAxes,fontsize=8)
    roles=['series_input','drain_resistor','feedback_resistor','output_load','unassigned_remainder'];names=['Series input resistors','R75 drain resistor','F75L feedback resistor','Output load','Unassigned remainder']
    colors=['#176c80','#4b8d79','#c65b0b','#c59844','#bec5cc'];bottom=np.zeros(2);data={(r['candidate'],r['component']):r for r in tables['noise-budget.csv']}
    for role,label,color in zip(roles,names,colors):
        heights=np.array([data[(n,role)]['variance_v2']*1e12 for n in ['R75','F75L']]);axes[1].bar([0,1],heights,bottom=bottom,color=color,width=.52,label=label);bottom+=heights
    for i,c in enumerate(result['cases']):axes[1].text(i,bottom[i]+8,f"Total {c['original_input_noise_rms_v']*1e6:.3f} µV RMS",ha='center',va='bottom',fontsize=8)
    axes[1].set(xticks=[0,1],xticklabels=['R75','F75L'],ylabel='Input-noise variance over 1 kHz–1 MHz (µV²)',title='Add powers, then take the square root',ylim=(0,490),xlim=(-.7,1.7))
    axes[1].legend(loc='upper left',fontsize=8)
    for ax in axes:ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
    fig.subplots_adjust(left=.09,right=.99,bottom=.23,top=.87,wspace=.35)
    fig.savefig(a.output/'resistor-noise-paths.svg',metadata={'Date':None});fig.savefig(a.output/'resistor-noise-paths.png',dpi=160,metadata={'Software':'Analog Design Studies'});plt.close(fig)
    print('PASS: selected resistor paths and unassigned remainder; no complete MOS budget')


if __name__=='__main__':main()
