"""MIT: recompute nominal/transfer MAE and paired intervals from public CSV."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.stats import t

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    with (a.data/'confirmation.csv').open() as f:rows=list(csv.DictReader(f))
    published=json.loads((a.data/'statistics.json').read_text());answer=[]
    for c in ['nominal','supply_097']:
        groups={}
        for n in ['R37','A37','R75','A75']:
            r=[x for x in rows if x['candidate']==n and x['condition']==c and x['numerical_usable']=='True' and x['terminal_supported']=='True']
            groups[n]={int(x['index']):float(x['error_v']) for x in r}
            v=np.array(list(groups[n].values()));old=published['summaries'][n+'/'+c]
            assert len(v)==old['usable_supported'] and abs(np.mean(abs(v))-old['mae_v'])<1e-15 and abs(v.std(ddof=1)-old['sd_error_v'])<1e-15
        for left,right in [('R37','R75'),('A37','A75'),('A37','R37'),('A75','R75')]:
            ids=sorted(groups[left].keys()&groups[right].keys());v=np.array([abs(groups[left][i])-abs(groups[right][i]) for i in ids])
            half=t.ppf(1-.05/8,len(ids)-1)*v.std(ddof=1)/np.sqrt(len(ids));ci=[float(v.mean()-half),float(v.mean()+half)]
            old=next(r for r in published['contrasts'] if r['condition']==c and r['left']==left and r['right']==right)
            assert np.max(abs(np.array(ci)-old['family95_t_interval_v']))<1e-15
            answer.append(dict(condition=c,left=left,right=right,n=len(ids),reduction_v=float(v.mean()),interval_v=ci))
    with a.output.open('x') as f:json.dump(dict(status='PASS',paired_comparisons=answer),f,indent=2);f.write('\n')
    print('All eight group summaries and eight paired comparisons match the published statistics to 1e-15 V.')

if __name__=='__main__':main()
