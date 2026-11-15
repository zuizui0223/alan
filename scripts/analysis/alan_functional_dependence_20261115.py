"""Source-locked diel exclusion contrasts, not an estimated effect of ALAN.

Primary contrast: 1 - day-only/open. Night-only treatments and floral traits
are optional. Paper bootstrap intervals do not propagate within-study error.
"""
from __future__ import annotations
import argparse, hashlib, json, platform
from pathlib import Path
import numpy as np
import pandas as pd

KEYS=['study_ID','plant_species','plant_cultivar','country','locality','Site_period',
      'coordinates_lat','coordinates_lon','year_start','year_end','month_start',
      'month_end','additional_treatment','treatment_effectiveness_metric',
      'effectiveness_units','treatmenteffectiveness_metric_comments']
ARMS=['day pollination','night pollination','open pollination','complete exclusion']
SOURCE_SHA='1d9c4c778206584130308b99ba21d902b0ccf1047118dfcf52abaec9d652e2aa'
TRAIT_SHA='9c4a00d74c954eb448330066071ee0c584a9df66d96911c011cef5cd1afc7173'

def match_panels(df:pd.DataFrame):
    required=set(KEYS+['treatment_condition','effectiveness_value'])
    if not required<=set(df.columns): raise ValueError('Required matching columns missing')
    x=df.copy()
    if 'exclude' in x: x=x[x['exclude'].eq('n')].copy()
    x=x[x.treatment_condition.isin(ARMS)].copy()
    for k in KEYS:x[k]=x[k].fillna('<NA>').astype(str).str.strip()
    x['effectiveness_value']=pd.to_numeric(x.effectiveness_value,errors='coerce')
    audit={'source_rows':int(len(df)),'relevant_arm_rows':int(len(x)),
           'ambiguous_panels':0,'zero_open_panels':0,'missing_primary_pair':0,
           'negative_or_nonfinite_means':0,'source_groups':0}
    out=[]
    for key,g in x.groupby(KEYS,dropna=False,sort=True):
        audit['source_groups']+=1
        if g.treatment_condition.duplicated().any():
            audit['ambiguous_panels']+=1;continue
        vals=g.set_index('treatment_condition').effectiveness_value.to_dict()
        day=vals.get(ARMS[0],np.nan); night=vals.get(ARMS[1],np.nan)
        op=vals.get(ARMS[2],np.nan); cl=vals.get(ARMS[3],np.nan)
        if not np.isfinite(day) or not np.isfinite(op):
            audit['missing_primary_pair']+=1;continue
        if day<0 or op<0 or (np.isfinite(night) and night<0) or (np.isfinite(cl) and cl<0):
            audit['negative_or_nonfinite_means']+=1;continue
        if op==0:audit['zero_open_panels']+=1;continue
        row=dict(zip(KEYS,key))
        row['panel_id']=hashlib.sha256(json.dumps(list(key),ensure_ascii=False).encode()).hexdigest()[:20]
        row.update(day_mean=float(day),open_mean=float(op),night_mean=float(night),closed_mean=float(cl))
        row['night_loss']=1-day/op
        row['night_alone']=night/op if np.isfinite(night) else np.nan
        row['day_loss']=1-night/op if np.isfinite(night) else np.nan
        row['night_share']=night/(day+night) if np.isfinite(night) and day+night>0 else np.nan
        row['closed_fraction']=cl/op if np.isfinite(cl) else np.nan
        row['night_loss_above_closed']=(op-day)/(op-cl) if np.isfinite(cl) and op>cl else np.nan
        row['all_three']=bool(np.isfinite(night))
        out.append(row)
    return pd.DataFrame(out),audit

def join_traits(panels:pd.DataFrame, traits:pd.DataFrame):
    fields=[c for c in ['bloom_period_simple','flower_shape','color','breeding_system'] if c in traits]
    aliases={}
    for _,r in traits.iterrows():
        payload={'trait_name':str(r['plant_species']).strip()}
        payload.update({c:str(r[c]).strip() if pd.notna(r[c]) else 'unknown' for c in fields})
        for name in [r.get('old_name'),r['plant_species']]:
            if pd.isna(name):continue
            name=str(name).strip()
            if name in aliases and aliases[name]!=payload:raise ValueError('Conflicting trait alias: '+name)
            aliases[name]=payload
    result=panels.copy()
    for c in ['trait_name']+fields:
        result[c]=result.plant_species.map(lambda s:aliases.get(s,{}).get(c,'unknown'))
    if 'bloom_period_simple' not in result:result['bloom_period_simple']='unknown'
    result['bloom_period_simple']=result.bloom_period_simple.where(result.bloom_period_simple.isin(['night','day','both']),'unknown')
    return result

def cluster_summary(df:pd.DataFrame,iterations:int=4000,seed:int=15, value:str='night_loss'):
    v=df.groupby('study_ID')[value].median().dropna().to_numpy(float)
    if not len(v):return {'panels':0,'studies':0}
    rng=np.random.default_rng(seed)
    draws=rng.choice(v,(iterations,len(v)),replace=True)
    return {'panels':int(df[value].notna().sum()),'studies':int(len(v)),
            'median_study_median':float(np.median(v)),
            'bootstrap95_median':np.quantile(np.median(draws,axis=1),[.025,.975]).tolist() if len(v)>1 else None,
            'mean_study_median':float(np.mean(v)),
            'bootstrap95_mean':np.quantile(np.mean(draws,axis=1),[.025,.975]).tolist() if len(v)>1 else None,
            'method':'publication-equal, median within each paper; percentile paper bootstrap',
            'not_inverse_variance_meta_or_measurement_error_CI':True}

def group_difference(df,group_a='night',group_b='day',iterations=4000,seed=15):
    g=df[df.bloom_period_simple.isin([group_a,group_b])].groupby(['study_ID','bloom_period_simple']).night_loss.median().unstack()
    if not {group_a,group_b}<=set(g.columns):return {'available':False}
    a=g[group_a].to_numpy();b=g[group_b].to_numpy()
    rng=np.random.default_rng(seed); inds=rng.integers(0,len(g),(iterations,len(g)))
    dif=np.nanmedian(a[inds],axis=1)-np.nanmedian(b[inds],axis=1)
    return {'available':True,'contrast':group_a+' minus '+group_b,
            'difference_of_study_equal_medians':float(np.nanmedian(a)-np.nanmedian(b)),
            'bootstrap95':np.quantile(dif[np.isfinite(dif)],[.025,.975]).tolist(),
            'studies_union':int(len(g)), 'studies_a':int(np.isfinite(a).sum()),'studies_b':int(np.isfinite(b).sum()),
            'both_groups_in_same_paper':int((np.isfinite(a)&np.isfinite(b)).sum()),
            'between_species_and_paper_comparison_not_causal_flowering_effect':True}

def diagnostic(df,threshold):
    b=df[df.night_share.notna() & (df.night_share<=.5)].copy()
    b['flag']=(b.night_loss>threshold).astype(float)
    c=df[df.night_alone.notna() & (df.night_alone>=.8)].copy()
    c['flag']=(c.night_loss<=threshold).astype(float)
    return {'threshold':threshold,
      'day_ge_night':{'panels':len(b),'night_loss_exceeds_threshold':int(b.flag.sum()),'studies':b.study_ID.nunique(),
          'paper_equal_conditional_fraction':float(b.groupby('study_ID').flag.mean().mean()) if len(b) else None},
      'night_alone_ge80pct_open':{'panels':len(c),'night_loss_at_or_below_threshold':int(c.flag.sum()),'studies':c.study_ID.nunique(),
          'paper_equal_conditional_fraction':float(c.groupby('study_ID').flag.mean().mean()) if len(c) else None}}

def endpoint_subtype(metric,unit):
    unit=str(unit).lower().strip()
    if metric=='fruit set':
        if 'proportion' in unit or 'fruit/flower' in unit or 'capsule/flowers' in unit:return 'fruit_proportion'
        if unit not in ['<na>','nan','mean','']:return 'fruit_count'
        return 'fruit_unspecified'
    if metric=='seed set':
        if 'proportion of flowers' in unit:return 'seed_formation_probability'
        if 'proportion' in unit and 'viable' in unit:return 'seed_viability_fraction'
        if 'ovule' in unit or 'fertilized' in unit:return 'seed_per_ovule'
        if 'plant' in unit:return 'seed_per_plant'
        if 'flower' in unit:return 'seed_per_flower'
        if any(x in unit for x in ['fruit','capsule','pod']):return 'seed_per_fruit'
        return 'seed_unspecified'
    return metric.replace(' ','_')

def analyze(source,traits,out,iterations=4000):
    source=Path(source);traits=Path(traits);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if hashlib.sha256(source.read_bytes()).hexdigest()!=SOURCE_SHA:raise ValueError('Diel source SHA mismatch')
    if hashlib.sha256(traits.read_bytes()).hexdigest()!=TRAIT_SHA:raise ValueError('Trait source SHA mismatch')
    raw=pd.read_csv(source);tr=pd.read_csv(traits);p,audit=match_panels(raw);p=join_traits(p,tr)
    p['endpoint_subtype']=[endpoint_subtype(m,u) for m,u in zip(p.treatment_effectiveness_metric,p.effectiveness_units)]
    metrics=['fruit set','seed set']
    result={'schema':'ALAN_diel_marginal_night_dependence_v1','source_sha256':SOURCE_SHA,'traits_sha256':TRAIT_SHA,
      'source_article_doi':'10.1111/ele.70036','source_commit':'1945b7c114ad87000608ce6018d0280001c04ef1',
      'source_rows':len(raw),'source_studies':raw.study_ID.nunique(),'source_species':raw.plant_species.nunique(),
      'matching_audit':audit,'primary_panels_all_endpoints':len(p),'primary_papers_all_endpoints':p.study_ID.nunique(),
      'primary_species_all_endpoints':p.plant_species.nunique(),'three_arm_panels':int(p.all_three.sum()),
      'coverage_by_endpoint':p.groupby('treatment_effectiveness_metric').size().to_dict(),
      'no_ALAN_treatment_effect_fitted':True,'no_future_ALAN_loss_predicted':True,
      'closed_exclusion_is_not_necessarily_autonomous_selfing':True,
      'exploratory_thresholds_not_preregistered':True,'seed':15,'bootstrap_iterations':iterations,
      'endpoint_results':{},'denominator_specific_results':{},
      'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}}
    for metric in metrics:
        q=p[p.treatment_effectiveness_metric.eq(metric)].copy()
        result['endpoint_results'][metric]={
          'species':q.plant_species.nunique(),'overall':cluster_summary(q,iterations),
          'bloom_groups':{b:cluster_summary(q[q.bloom_period_simple.eq(b)],iterations) for b in ['night','day','both','unknown']},
          'bloom_contrast':group_difference(q,iterations=iterations),
          'threshold_diagnostics':[diagnostic(q,t) for t in [.1,.2,.3]],
          'closed_adjusted':cluster_summary(q,iterations,value='night_loss_above_closed'),
          'closed_adjusted_with_closed_below80pct_open':cluster_summary(q[q.closed_fraction<.8],iterations,value='night_loss_above_closed'),
          'exclude_known_ALAN_manipulation_paper':cluster_summary(q[~q.study_ID.eq('Macgregor_2019')],iterations),
          'negative_deficit_panels':int((q.night_loss<0).sum()),
          'open_smaller_than_half_day_panels':int((q.night_loss < -1).sum()),
          'unknown_bloom_panels':int(q.bloom_period_simple.eq('unknown').sum()),
          'source_units':q.effectiveness_units.value_counts().to_dict()}
    for subtype,q in p[p.treatment_effectiveness_metric.isin(metrics)].groupby('endpoint_subtype'):
        result['denominator_specific_results'][subtype]={'summary':cluster_summary(q,iterations),'diagnostics20':diagnostic(q,.2),'species':q.plant_species.nunique()}
    repro=p[p.treatment_effectiveness_metric.isin(metrics)]
    result['reproductive_coverage']={'panels':len(repro),'papers':repro.study_ID.nunique(),'species':repro.plant_species.nunique()}
    keep=['endpoint_subtype','panel_id','study_ID','plant_species','trait_name','treatment_effectiveness_metric','effectiveness_units','bloom_period_simple',
          'flower_shape','color','breeding_system','night_loss','night_alone','day_loss','night_share','closed_fraction','night_loss_above_closed','all_three']
    p[[c for c in keep if c in p]].to_csv(out/'functional_dependence_panels.csv',index=False)
    (out/'functional_dependence_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,ensure_ascii=False,allow_nan=False))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',required=True);parser.add_argument('--traits',required=True);parser.add_argument('--out',required=True);parser.add_argument('--iterations',type=int,default=4000)
    args=parser.parse_args();analyze(args.source,args.traits,args.out,args.iterations)
