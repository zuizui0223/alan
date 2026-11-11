#!/usr/bin/env python3
"""Replay the source-gated GloBI cross-predator flower photo/ALAN result.

Requires PRIVATE original plant photo-cell×year numeric source (do not commit).
Positive global GloBI source links are NOT local effective pollination or
documented predator occupancy. Site-level records remain private.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import t as student_t

EXPECTED_SHA='32000250071106f81f203ca8c40e2faa67d02152e06afd442d6f93df5f69d233'
GROUPS=('bats','geckos','centipedes','spiders','birds','mantises')
SCOPES=('all','non_inaturalist_trophich')
RELATIONS=('flower_either','db_pollinates','visit')

def check(ok,message):
    if not bool(ok):raise ValueError(message)

def load(path,expected_sha=EXPECTED_SHA):
    if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected_sha:
        raise ValueError('Source identity table SHA256 mismatch')
    d=pd.read_csv(path,low_memory=False,dtype={'target_id':str,'land_id':str})
    required={'target_id','year','region5','land_class','land_id','viirs','population_density',
              'photo_attempts','n_species','temperature','precipitation','elevation'}
    counts={f'count_{t}_{g}_{r}' for t in SCOPES for g in GROUPS for r in RELATIONS}
    check(not (required|counts)-set(d),'Source schema mismatch')
    check(d.target_id.is_unique,'Duplicated original photo site-year source')
    check(set(d.year)=={2019,2024} and set(d.land_class)=={'island','mainland'},
          'Unexpected source year or land class')
    for col in counts:
        check(d[col].ge(0).all() and d[col].le(d.n_species).all(),
              'GloBI-positive recorded species exceed photographed accepted species')
    d['island']=d.land_class.eq('island').astype(int)
    d['region_year']=d.year.astype(str)+'|'+d.region5.astype(str)
    return d

def fit(d,response):
    s=d.loc[d.photo_attempts.ge(3)&d.n_species.ge(2)].copy()
    support=s.groupby('region_year').land_class.nunique()
    s=s.loc[s.region_year.isin(support.index[support.eq(2)])].copy()
    pos=s.loc[s[response].gt(0)]
    isl=pos.loc[pos.island.eq(1)];main=pos.loc[pos.island.eq(0)]
    counts={'positive_island_cells':len(isl),'positive_mainland_cells':len(main),
            'positive_island_geographic_blocks':int(isl.region5.nunique()),
            'source_island_site_years':int(s.island.sum()),
            'source_mainland_site_years':int((1-s.island).sum()),
            'source_geographic_blocks':int(s.region5.nunique())}
    if not(len(isl)>=20 and len(main)>=40 and isl.region5.nunique()>=10 and s.region5.nunique()>=15):
        return {'estimable':False,'reason':'predeclared_source_positive_photo_geography_gate','support':counts}
    y=100*s[response].to_numpy(float)/s.n_species.to_numpy(float)
    a=s.copy()
    a['loglight']=np.log1p(s.viirs.clip(lower=0))
    a['island_light']=s.island*a.loglight
    a['logpop']=np.log1p(s.population_density.clip(lower=0))
    a['logphotos']=np.log1p(s.photo_attempts)
    a['logspecies']=np.log1p(s.n_species)
    a['lograin']=np.log1p(s.precipitation.clip(lower=0))
    a['elevation_km']=s.elevation/1000
    cols=['loglight','island','island_light','logpop','logphotos','logspecies',
          'temperature','lograin','elevation_km']
    X=a[cols].to_numpy(float)
    check(np.isfinite(X).all() and np.isfinite(y).all(),'Nonfinite source model inputs')
    group=a.groupby(['region_year','land_class'],observed=True).size()
    keys=pd.MultiIndex.from_frame(a[['region_year','land_class']])
    w=group.groupby(level=0).min().reindex(a.region_year).to_numpy(float)/group.reindex(keys).to_numpy(float)
    check(np.isfinite(w).all() and (w>0).all(),'Invalid year/geography land-class weights')
    code,regions=pd.factorize(a.region_year,sort=True)
    mat=np.column_stack([X,y])
    totals=np.bincount(code,weights=w)
    for j in range(mat.shape[1]):
        mean=np.bincount(code,weights=w*mat[:,j])/totals
        mat[:,j]=(mat[:,j]-mean[code])*np.sqrt(w)
    X,Y=mat[:,:-1],mat[:,-1]
    beta,_,rank,_=np.linalg.lstsq(X,Y,rcond=None)
    if rank!=len(cols):return {'estimable':False,'reason':'rank_deficient_year_region_model','support':counts}
    residual=Y-X@beta
    geography,groups=pd.factorize(a.region5,sort=True)
    scores=X*residual[:,None]
    sums=np.column_stack([np.bincount(geography,weights=scores[:,i],minlength=len(groups)) for i in range(len(cols))])
    bread=np.linalg.inv(X.T@X)
    cov=bread@(sums.T@sums)@bread*len(groups)/(len(groups)-1)*(len(Y)-1)/(len(Y)-len(cols))
    j=cols.index('island_light')
    se=float(np.sqrt(cov[j,j])) if cov[j,j]>0 else None
    tcrit=float(student_t.ppf(.975,len(groups)-1))
    return {'estimable':True,'beta_island_minus_mainland_pp_per_log1pVIIRS':float(beta[j]),
            'nominal_95CI':[float(beta[j]-tcrit*se),float(beta[j]+tcrit*se)] if se else None,
            'nominal_p':float(2*student_t.sf(abs(beta[j]/se),len(groups)-1)) if se else None,
            'source_site_years':len(a),'spatial_blocks':len(groups),'support':counts,
            'response':'share of photographed accepted plants with positive moth AND other guild flower records SOMEWHERE in GloBI'}

def run(private_source,expected_sha=EXPECTED_SHA):
    d=load(private_source,expected_sha)
    return {'schema':'alan_GloBI_v011_crossguild_photo_light_private_source_replay_v1',
       'private_source_SHA256':expected_sha,'years':[2019,2024],
       'GloBI_global_visitor_positive_not_local_pollination':True,
       'not_predator_abundance_or_causal_ALAN':True,
       'models':[{'source_scope':scope,'animal_guild':group,'GloBI_relation':rel,
           'test':fit(d,f'count_{scope}_{group}_{rel}')}
           for scope in SCOPES for group in GROUPS for rel in RELATIONS]}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--private-source',type=Path,required=True)
    p.add_argument('--expected-sha256',default=EXPECTED_SHA)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    check(not a.output.exists(),'Refuse overwrite')
    result=run(a.private_source,a.expected_sha256)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'estimated':[(x['source_scope'],x['animal_guild'],x['GloBI_relation']) for x in result['models'] if x['test']['estimable']]}))
if __name__=='__main__':main()
