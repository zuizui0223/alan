#!/usr/bin/env python3
"""Source-verified, within-one-pond ALAN colour x habitat sensitivity.

No added independent study, no causal wet-vs-dry plant moderator. Read existing
30 authentic trap units and matching 193 family rows; never use pooled effects
as input. No pseudo-count or taxon absence imputation.
"""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, math
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'database/alan_experiments/parkinson2023_sources'
COLORS=('dark','green','purple','red','white')
HABITATS=('littoral','riparian')

def read(path):
    with path.open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))

def source():
    u=read(BASE/'trap_units.csv'); tx=read(BASE/'taxon_counts.csv')
    if len(u)!=30 or len(tx)!=193: raise ValueError('source rows unexpectedly changed')
    per=defaultdict(dict); group_totals=defaultdict(float)
    for t in tx:
        key=t['unit_id']; fam=t['family_source_spelling']; count=float(t['count'])
        if not math.isfinite(count) or count<0 or fam in per[key]:raise ValueError('invalid duplicate family/negative count')
        per[key][fam]=count;group_totals[fam]+=count
    keys=set()
    for row in u:
        k=row['unit_id'];h=row['habitat'];c=row['colour'];n=float(row['all_invertebrate_capture'])
        if k in keys or k not in per or not math.isfinite(n) or abs(sum(per[k].values())-n)>1e-9:
            raise ValueError('taxon count mismatch or duplicate trap')
        keys.add(k)
        if (h,c) not in itertools.product(HABITATS,COLORS):raise ValueError('unexpected experimental arm')
    if set(keys)!=set(per):raise ValueError('unmatched taxon rows')
    if any(sum(x['habitat']==h and x['colour']==c for x in u)!=3 for h,c in itertools.product(HABITATS,COLORS)):
        raise ValueError('expected randomised balanced n=3 within 10 cells')
    return u,per,dict(group_totals)

def mean(x):return sum(x)/len(x)

def interaction(u,per,excluded=(),transform='raw'):
    rows=[]
    for x in u:
        y=float(x['all_invertebrate_capture'])-sum(per[x['unit_id']].get(t,0) for t in excluded)
        if y<0:raise ValueError('negative residual count')
        if transform=='log1p':y=math.log1p(y)
        elif transform!='raw':raise ValueError('unsupported transformation')
        rows.append((x['habitat'],x['colour'],y))
    grand=mean([z[2] for z in rows]); hm={h:mean([z[2] for z in rows if z[0]==h]) for h in HABITATS}
    cm={c:mean([z[2] for z in rows if z[1]==c]) for c in COLORS}
    si=se=0.
    for h,c in itertools.product(HABITATS,COLORS):
        y=[z[2] for z in rows if z[0]==h and z[1]==c]
        m=mean(y);si+=3*(m-hm[h]-cm[c]+grand)**2;se+=sum((v-m)**2 for v in y)
    if se<=0:raise ValueError('nonpositive residual SS')
    f=(si/4)/(se/20)
    x=5/(f+5)
    # Exact F(4,20) survival for balanced two-way ANOVA, beta(10,2) tail.
    p=x**10*(11-10*x)
    return {'F':f,'p_nominal':p,'df_num':4,'df_den':20}

def violet_moderator(u,per,ex=()):
    means={}
    for h,c in itertools.product(HABITATS,('purple','dark')):
        yy=[float(x['all_invertebrate_capture'])-sum(per[x['unit_id']].get(t,0) for t in ex)
            for x in u if x['habitat']==h and x['colour']==c]
        means[f'{h}_{c}']=mean(yy)
    near=means['littoral_purple']/means['littoral_dark']
    far=means['riparian_purple']/means['riparian_dark']
    return {'violet_littoral_RR':near,'violet_riparian_RR':far,'ratio_of_RRs':near/far,'cell_means':means}

def run():
    u,per,total=source(); sorted_f=sorted(total,key=lambda x:(-total[x],x));top2=sorted_f[:2]
    if top2!=['caenidae','leptoceridae']:raise ValueError('unexpected dominant family ranks')
    profiles={key:interaction(u,per,ex,tr) for key,ex,tr in (
        ('full_raw',(),'raw'),('full_log1p',(),'log1p'),
        ('top2_excluded_raw',top2,'raw'),('top2_excluded_log1p',top2,'log1p'))}
    single=[dict(family=f,**interaction(u,per,[f],'log1p')) for f in sorted(total)]
    pairs=[]
    for f1,f2 in itertools.combinations(sorted(total),2):
        pairs.append(dict(family_1=f1,family_2=f2,**interaction(u,per,[f1,f2],'log1p')))
    pairs.sort(key=lambda x:(x['F'],x['family_1'],x['family_2']))
    non_sig=[x for x in pairs if x['p_nominal']>=0.05]
    inputhash={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (BASE/'trap_units.csv',BASE/'taxon_counts.csv')}
    return {
        'study_doi':'10.1098/rstb.2022.0361','programme_id':'SevenPonds2022','year_of_original_field':2022,
        'true_independent_lighting_units':30,'one_pond_one_night':True,
        'original_taxa_rows':sum(len(v) for v in per.values()),'families':len(total),
        'original_invertebrates':sum(float(x['all_invertebrate_capture']) for x in u),
        'top_two_families':[{ 'name':x,'count':total[x]} for x in top2],
        'top_two_fraction':sum(total[x] for x in top2)/sum(total.values()),
        'F4_20_profiles':profiles,'violet_full':violet_moderator(u,per),
        'violet_excluding_top_two':violet_moderator(u,per,top2),
        'leave_one_family_out':sorted(single,key=lambda x:x['F']),
        'exhaustive_two_family_exclusion_scan':pairs,
        'pairs_p_at_least_0_05':len(non_sig),'pairs_total':len(pairs),
        'pairs_nonsignificant':[[x['family_1'],x['family_2']] for x in non_sig],
        'source_hashes':inputhash,
        'null_caveat':'All F p-values are original ANOVA-model/normality tests: not randomization-exact, not family-wise adjusted; pair removal post hoc. Species-trait moisture and independent habitat replication unavailable.',
        'carannante_design_idealised_small_n':{'study':'10.1002/2688-8319.12053','independent_lit_dark_pairs':4,'distance_classes':[3,40],'spatial_blocks':2,'two_sided_min_p_if_two_blocks_exchangeable':2/2**2,'two_sided_min_p_if_four_pairs_fully_exchangeable':2/math.comb(4,2),'note':'Hypothetical randomization-support bound, not an observed exact p or proof that distance was randomized. Cannot confirm a habitat moderator with only 2 blocks.'}
    }

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=ROOT/'database/alan_experiments/arthropod_expansion_20261122/sevenponds_family_sensitivity');args=ap.parse_args()
    d=run();args.outdir.mkdir(parents=True,exist_ok=True)
    summary={k:v for k,v in d.items() if k not in ('exhaustive_two_family_exclusion_scan','leave_one_family_out')}
    (args.outdir/'summary_20261123.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for name,records in [('every_two_family_removal_20261123.csv',d['exhaustive_two_family_exclusion_scan']),('single_family_removal_20261123.csv',d['leave_one_family_out'])]:
        with (args.outdir/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader()
            for row in records:w.writerow(row)
    print('SOURCE_REPLAY',json.dumps({'n_traps':30,'n_families':d['families'],'n_pairs':d['pairs_total'], 'n_nonsignificant':d['pairs_p_at_least_0_05'], 'full_F':d['F4_20_profiles']['full_raw']['F'],'filtered_F':d['F4_20_profiles']['top2_excluded_log1p']['F']}))

if __name__=='__main__':main()
