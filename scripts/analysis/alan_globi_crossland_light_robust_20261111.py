#!/usr/bin/env python3
"""Exploratory GloBI-positive photographed-plant fraction versus ALAN, island/mainland.

Input is PRIVATE source-derived site-year aggregate with authentic photo/WCVP/
same-year satellite identity. Export AGGREGATES ONLY; do not publish original
site/taxon/photographer locations. The response is NOT observed local moth
visitation, flower population abundance, nocturnal pollination or plant decline.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

SCOPES=('all','non_inat','non_inat_trophich','non_inat_visits')

def cleaned(input_path):
    s=pd.read_csv(input_path)
    needs=['target_id','land_class','n_species','viirs_median','pop_density',
           'photo_attempts','year','region5','bio5','rain','elevation']
    for scope in SCOPES:
        needs.extend([f'{scope}_moth',f'{scope}_butterfly'])
    bad=set(needs)-set(s)
    if bad: raise ValueError(f'Missing source columns: {bad}')
    if not s.target_id.is_unique: raise ValueError('Repeated original photo cell-year')
    if set(s.land_class)!={'island','mainland'}: raise ValueError('Unknown land class')
    for scope in SCOPES:
        for guild in ('moth','butterfly'):
            x=s[f'{scope}_{guild}']
            if not (x.ge(0)&x.le(s.n_species)).all():
                raise ValueError('GloBI-positive photo-taxon count exceeds actual photographed species')
    s['island']=s.land_class.eq('island').astype(float)
    s['region_year']=s.year.astype(str)+'|'+s.region5.astype(str)
    s['loglight']=np.log1p(s.viirs_median.clip(lower=0))
    s['logpop']=np.log1p(s.pop_density.clip(lower=0))
    s['logphotos']=np.log1p(s.photo_attempts)
    s['logspecies']=np.log1p(s.n_species)
    s['lograin']=np.log1p(s.rain.clip(lower=0))
    s['elevation_km']=s.elevation/1000
    s['light_island']=s.loglight*s.island
    return s

def prepare(s, attempt, species, *, same_year_block=True, light_overlap=False):
    g=s.loc[s.photo_attempts.ge(attempt)&s.n_species.ge(species)].copy()
    if same_year_block:
        counts=g.groupby('region_year').land_class.nunique()
        g=g.loc[g.region_year.isin(counts.index[counts.eq(2)])].copy()
    if light_overlap:
        low=g.groupby(['region_year','land_class']).loglight.min().unstack()
        high=g.groupby(['region_year','land_class']).loglight.max().unstack()
        lo,hi=low.max(axis=1),high.min(axis=1)
        valid=lo.index[hi.gt(lo)]
        g=g.loc[g.region_year.isin(valid)].copy()
        g=g.loc[g.loglight.ge(g.region_year.map(lo)) &
                g.loglight.le(g.region_year.map(hi))].copy()
        counts=g.groupby('region_year').land_class.nunique()
        g=g.loc[g.region_year.isin(counts.index[counts.eq(2)])].copy()
    return g

def fit(g,scope,guild,*,balance=False,logspecies=True):
    if guild=='moth_minus_butterfly':
        y=100*(g[f'{scope}_moth']-g[f'{scope}_butterfly'])/g.n_species
    elif guild in ('moth','butterfly'):
        y=100*g[f'{scope}_{guild}']/g.n_species
    else: raise ValueError('Unknown positive evidence guild')
    cols=['island','loglight','light_island','logpop','logphotos']
    if logspecies:cols.append('logspecies')
    cols+=['bio5','lograin','elevation_km']
    a=g[cols+['region_year','region5','land_class']].copy()
    a['response']=y.to_numpy()
    a=a.replace([np.inf,-np.inf],np.nan).dropna().reset_index(drop=True)
    if len(a)<50 or a.region5.nunique()<10 or a.land_class.nunique()<2:
        return {'estimable':False,'reason':'inadequate_complete_cases','n':len(a)}
    counts=a.groupby(['region_year','land_class'],observed=True).size().unstack(fill_value=0)
    valid=counts.index[counts.min(axis=1)>0]
    a=a[a.region_year.isin(valid)].copy()
    if len(a)<50:return {'estimable':False,'reason':'no_complete_case_land_overlap'}
    if balance:
        sizes=a.groupby(['region_year','land_class'],observed=True).size()
        index=pd.MultiIndex.from_frame(a[['region_year','land_class']])
        weights=sizes.groupby(level=0).min().reindex(a.region_year).to_numpy(float)/sizes.reindex(index).to_numpy(float)
    else:
        weights=np.ones(len(a))
    z=a[cols+['response']].to_numpy(float)
    group,unique=pd.factorize(a.region_year,sort=True)
    wsum=np.bincount(group,weights=weights,minlength=len(unique))
    for j in range(z.shape[1]):
        means=np.bincount(group,weights=weights*z[:,j],minlength=len(unique))/wsum
        z[:,j]=(z[:,j]-means[group])*np.sqrt(weights)
    x,y=z[:,:-1],z[:,-1]
    beta,_,rank,_=np.linalg.lstsq(x,y,rcond=1e-12)
    if rank!=len(cols):return {'estimable':False,'reason':'fixed_effect_collinearity'}
    err=y-x@beta; scores=x*err[:,None]
    clusters,geokeys=pd.factorize(a.region5,sort=True);G=len(geokeys)
    sums=np.column_stack([np.bincount(clusters,weights=scores[:,j],minlength=G) for j in range(x.shape[1])])
    bread=np.linalg.inv(x.T@x)
    cov=bread@sums.T@sums@bread*(G/(G-1))*((len(a)-1)/(len(a)-len(cols)))
    ix=cols.index('light_island')
    se=math.sqrt(float(cov[ix,ix])) if cov[ix,ix]>0 else None
    tcrit=float(student_t.ppf(.975,G-1))
    return {'estimable':True,'scope':scope,'guild':guild,'balanced_year_region_landclass':balance,
            'n_site_years':len(a),'n_island_site_years':int(a.land_class.eq('island').sum()),
            'n_mainland_site_years':int(a.land_class.eq('mainland').sum()),
            'source_5degree_clusters':G,'year_region_groups':a.region_year.nunique(),
            'beta_island_minus_mainland_slope_pp_per_log1p':float(beta[ix]),'SE_cluster_5deg':se,
            '95CI_cluster_5deg':[float(beta[ix]-tcrit*se),float(beta[ix]+tcrit*se)] if se else None,
            'p_nominal':float(2*student_t.sf(abs(beta[ix]/se),G-1)) if se else None,
            'beta_mainland_pp_per_log1p':float(beta[cols.index('loglight')]),
            'beta_island_pp_per_log1p':float(beta[cols.index('loglight')]+beta[ix]),
            'GloBI_both_recorded_guilds_cancel_in_contrast':guild=='moth_minus_butterfly',
            'no_true_pollinator_absence_inferred':True,'not_local_effective_nocturnal_pollination':True}

def subset_support(g):
    d=g.replace([np.inf,-np.inf],np.nan).dropna(subset=['loglight','logpop','bio5','lograin','elevation_km','logphotos','logspecies'])
    group=d.groupby('land_class').agg(n=('target_id','size'),median_light=('viirs_median','median'),
        median_accepted_species=('n_species','median'),median_photo_attempts=('photo_attempts','median'))
    support=d.groupby(['region_year','land_class'],observed=True).size().unstack(fill_value=0)
    return {'raw_matched_rows':len(g),'complete_case_rows':len(d),
            'class_counts':group.to_dict(orient='index'),
            'year_region_groups_both_complete':int(support.min(axis=1).gt(0).sum()),
            'groups_with_3_sites_each':int((support.min(axis=1)>=3).sum()),
            'distinct_5deg_regions':int(d.region5.nunique())}

def run(input_path):
    s=cleaned(input_path)
    out={'status':'EXPLORATORY_NO_PROVEN_ALAN_POLLINATOR_MECHANISM',
         'source':'61,967-claim bounded GloBI, original photos, 2019/2021/2023/2024/2025 same-year 4-month VIIRS',
         'source_support':{},'models':[],'radiance_overlap_support':{},
         'identity_level_input_private':True,'nominal_p_values_unadjusted':True}
    for attempt,species in [(3,2),(5,3)]:
        cohort=f'{attempt}photos_{species}species'
        matched=prepare(s,attempt,species)
        out['source_support'][cohort]=subset_support(matched)
        for scope in SCOPES:
            for guild in ('moth','butterfly','moth_minus_butterfly'):
                for balance in (False,True):
                    out['models'].append({'cohort':cohort,'radiance_overlap':False,
                                         **fit(matched,scope,guild,balance=balance)})
        overlap=prepare(s,attempt,species,light_overlap=True)
        out['radiance_overlap_support'][cohort]=subset_support(overlap)
        for scope in ('all','non_inat','non_inat_trophich'):
            for balance in (False,True):
                out['models'].append({'cohort':cohort,'radiance_overlap':True,
                    **fit(overlap,scope,'moth_minus_butterfly',balance=balance)})
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--private-source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():raise ValueError('Refuse to overwrite original results')
    result=run(a.private_source)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'source_support':result['source_support'],'n_fits':len(result['models'])},ensure_ascii=False))

if __name__=='__main__':main()
