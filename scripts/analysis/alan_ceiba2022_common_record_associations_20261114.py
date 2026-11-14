#!/usr/bin/env python3
"""Descriptive common-record associations, not a causal SEM or trial reanalysis.

No p-value or independent-tree CI is inferred without verified design/cluster
metadata. Zero visits remain observed zeros; missing visit durations stay NA.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, math, statistics
from pathlib import Path
from scripts.analysis.alan_crossstudy_source_intake_20261114 import get
from scripts.analysis.alan_ceiba2022_original_source_20261114 import URL

SHA256='f5b4dca48da2fc9736e8e2441491f365d9916f595ad129eebf3dd5809a4190d5'
ENDPOINTS=('VisitsN','Pload','TubesP','GermPolen','Fruitset','SeedSet')
EXPOSURES=('Lux','PFD')


def ranks(values):
    if not values or not all(math.isfinite(x) for x in values):raise ValueError('Invalid numeric data')
    order=sorted(range(len(values)),key=lambda i:values[i]);out=[0.0]*len(values);start=0
    while start<len(values):
        end=start+1
        while end<len(values) and values[order[end]]==values[order[start]]:end+=1
        rank=(start+1+end)/2
        for i in order[start:end]:out[i]=rank
        start=end
    return out


def spearman(x,y):
    if len(x)!=len(y) or len(x)<3:raise ValueError('Need matched samples')
    if len(set(x))<2 or len(set(y))<2:raise ValueError('Constant data')
    return statistics.correlation(ranks(x),ranks(y))


def numerical(rows,key):
    vals=[float(r[key]) for r in rows]
    if not all(math.isfinite(v) for v in vals):raise ValueError('Nonfinite source column '+key)
    return vals


def source_ratio_check(rows,numerator,denominator,reported):
    a=numerical(rows,numerator);b=numerical(rows,denominator);c=numerical(rows,reported)
    if any(v<=0 for v in b):raise ValueError('Nonpositive denominator')
    errors=[abs(x/y-z) for x,y,z in zip(a,b,c)]
    return {'numerator':numerator,'denominator':denominator,'reported':reported,
            'maximum_absolute_error':max(errors),'within_0p001':sum(e<=.001 for e in errors),
            'records':len(rows),'identity_does_not_prove_shared_flower_sampling':True}


def association(x,y):
    r=spearman(x,y);loo=[]
    for i in range(len(x)):
        xx=x[:i]+x[i+1:];yy=y[:i]+y[i+1:]
        loo.append(spearman(xx,yy))
    return {'n':len(x),'spearman_rho':r,'leave_one_record_out_range':[min(loo),max(loo)],
            'leave_one_out_not_independent_replications':True,
            'p_and_cluster_CI_not_estimated':True}


def main(out):
    raw=get(URL,1_000_000)
    if hashlib.sha256(raw).hexdigest()!=SHA256:raise ValueError('Original public CSV changed')
    reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
    rows=[{key.strip():v.strip() for key,v in row.items()} for row in reader]
    if len(rows)!=43 or len({r['ID'] for r in rows})!=43:raise ValueError('Original row identity drift')
    result={'schema':'ALAN_Ceiba2022_same43record_light_visits_reproduction_descriptive_v1',
        'original_source_doi':'10.1007/s00114-021-01783-5','source_sha256':SHA256,
        'source_records':43,'unique_source_identifiers':43,
        'raw_records_exported':False,'design':'observational light/noise gradient',
        'new_associations_computed':True,'original_author_SEM_reproduced':False,
        'independent_experimental_unit_not_established_from_ID_alone':True,
        'source_ratio_checks':[
            source_ratio_check(rows,'SeedsN','Ovules','SeedSet'),
            source_ratio_check(rows,'FruitN','FlowersN','Fruitset'),
            source_ratio_check(rows,'TubesP','Pload','GermPolen')],
        'light_outcome_associations':{},'visits_downstream_associations':{},
        'source_outcome_structure':{
            'zero_recorded_visit_rows':sum(float(r['VisitsN'])==0 for r in rows),
            'missing_visit_duration_rows':sum(r['VisitsD'] in {'','NA','NaN'} for r in rows),
            'missing_duration_only_where_no_recorded_visits':all(
                float(r['VisitsN'])==0 for r in rows if r['VisitsD'] in {'','NA','NaN'}),
            'zero_fruitN_with_positive_SeedsN_rows':sum(float(r['FruitN'])==0 and float(r['SeedsN'])>0 for r in rows),
            'do_not_treat_seed_per_fruit_as_total_seeds_or_impute_seedN_zero':True}}
    for exposure in EXPOSURES:
        x=numerical(rows,exposure)
        result['light_outcome_associations'][exposure]={end:association(x,numerical(rows,end)) for end in ENDPOINTS}
    visits=numerical(rows,'VisitsN')
    for endpoint in ENDPOINTS[1:]:
        result['visits_downstream_associations'][endpoint]=association(visits,numerical(rows,endpoint))
    result['exposure_covariation']={col:association(numerical(rows,'Lux'),numerical(rows,col))
        for col in ('PFD','SoundFreq','SoundD','VegetationHa','FlowersTreeNigth')}
    result['limits']=[
        'Correlations use the same 43 archived records; no significance claim, random assignment, cluster replication or global generalization.',
        'Night visit frequency, stigma pollen, pollen tubes, fruit set and seed/ovule proportion are different endpoints.',
        'Noise, tree size, flowering intensity, vegetation and other urban drivers may confound these bivariate associations.',
        'Seed numbers remain positive in records with zero marked-fruit counts: those columns cannot be forced into one same-flower mass balance.',
        'A negative light-visit association with weak reproductive associations is not evidence that bats are dispensable or that light directly improves fitness.',
        'No island-mainland comparison, predator dual role or flowering-visitor mismatch is observed in this CSV.'
    ]
    target=Path(out);target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists():raise ValueError('Refuse overwrite')
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('CEIBA_COMMON_RECORD_ANALYSIS',json.dumps(result,ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    main(p.parse_args().out)
