#!/usr/bin/env python3
"""Post-hoc *sensitivity* of 6-header EERES 2021 source ratios to eight tested outcomes.

Uses analytic Student-t two-sided p for df=5; publication decisions are not
preregistered. The two eight-test families are overlapping in source and are
not independent repetitions or confirmation tests.
"""
import argparse,csv,json,math
from pathlib import Path

def t5_two_sided_p(est,se):
    """Analytic t5 tail via beta(5/2,1/2) trigonometric integral."""
    if se<=0 or not math.isfinite(se):raise ValueError('SOURCE_STOP_NONPOSITIVE_SE')
    t=abs(est)/se
    theta=math.asin(math.sqrt(5/(5+t*t)))
    p=(2*theta-4*math.sin(2*theta)/3+math.sin(4*theta)/6)/math.pi
    return max(0.0,min(1.0,p))

def corrections(pvals):
    m=len(pvals);order=sorted(range(m),key=lambda i:pvals[i]);holm=[None]*m;bh=[None]*m
    ma=0.0
    for r,i in enumerate(order):
        ma=max(ma, min(1.0,pvals[i]*(m-r)))
        holm[i]=ma
    mi=1.0
    for r in range(m-1,-1,-1):
        i=order[r]
        mi=min(mi,pvals[i]*m/(r+1),1.0)
        bh[i]=mi
    return holm,bh

def derive(data):
    original=[]
    for endpoint,metrics in data['effect_summaries'].items():
        for metric,rows in metrics.items():
            for v in rows:
                original.append({'family':'original_8_effects','term':f'{endpoint}|{metric}|{v["treatment"]}/{v["control"]}',
                    'ln_point':v['lnRR'],'SE_block':v['SE_paired_HT_delta']})
    comparison=[]
    for v in data['joint_shared_HT']['stage_contrasts']+data['joint_shared_HT']['light_by_flow_modifiers']:
        comparison.append({'family':'joint_8_posthoc_contrasts','term':v['kind']+'|'+v['numerator_lnRR']+'−'+v['denominator_lnRR'],
             'ln_point':v['difference_lnRR'],'SE_block':v['SE_shared_HT']})
    if len(original)!=8 or len(comparison)!=8 or data['joint_shared_HT']['independent_header_blocks']!=6:
        raise ValueError('SOURCE_STOP_NOT_8_AND_8_OR_NOT_SIX_HEADERS')
    result=[]
    for rows in (original,comparison):
        pvals=[t5_two_sided_p(v['ln_point'],v['SE_block']) for v in rows]
        holm,bh=corrections(pvals)
        result.extend([{**r,'df':5,'p_unadjusted':pvals[i], 'p_holm_within_8':holm[i], 'q_BH_within_8':bh[i],
                'status':'EXPLORATORY_NO_CONFIRMATORY_MULTIPLICITY_POLICY'} for i,r in enumerate(rows)])
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--csv',required=True);args=p.parse_args()
    data=json.loads(Path(args.input).read_text(encoding='utf-8'))
    rows=derive(data);d=Path(args.csv);d.parent.mkdir(parents=True,exist_ok=True)
    with d.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('MULTIPLICITY_AUDIT',len(rows),'original min Holm',min(r['p_holm_within_8'] for r in rows[:8]),'joint min BH',min(r['q_BH_within_8'] for r in rows[8:]))

if __name__=='__main__':main()
