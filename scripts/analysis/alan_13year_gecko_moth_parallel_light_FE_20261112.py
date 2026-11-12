#!/usr/bin/env python3
"""Parallel positive-only GBIF gecko and moth photo ratios against previous-year PCNL.

Same GSHHG island/year/0.05-degree cells, original SHA-pinned 13-year sources.
No animal population abundance, actual predator-prey capture, pollination or
species absence is established by a community photo record.
"""
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.sparse import coo_matrix,csc_matrix,hstack
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsmr
from scipy.stats import t as student_t

from alan_13year_moth_gecko_repeat_site_FE_20261112 import read_pair,repeated_sites
from alan_13year_moth_gecko_same_island_source_ALAN_20261112 import cohort

def fit(a,outcome='gecko',two_way=False,weight='island_equal'):
    a=a.copy();n=len(a)
    if n<100 or a.island.nunique()<15:
        return {'estimable':False,'reason':'too_few_original_island_photo_sites'}
    a['response']=100*a.gecko_share if outcome=='gecko' else a.moth_share
    predictors=['loglight','log_insect','log_reptile']+([] if two_way else ['latitude'])
    fixed=['site','islandyear'] if two_way else ['islandyear']
    if weight=='island_equal':
        w=1/a.groupby('island').island.transform('size').to_numpy(float)
    elif weight=='cell_equal':w=np.ones(n)
    elif weight=='islandyear_equal':
        w=1/a.groupby('islandyear').islandyear.transform('size').to_numpy(float)
    else:raise ValueError('Unknown source weight')
    sw=np.sqrt(w)
    dims=[];codes=[]
    for col in fixed:
        c,keys=pd.factorize(a[col],sort=True)
        dims.append(len(keys));codes.append(c)
    Z=hstack([csc_matrix((sw,(np.arange(n),c)),shape=(n,k))
              for c,k in zip(codes,dims)],format='csc')
    V=a[predictors+['response']].to_numpy(float)*sw[:,None]
    if not np.isfinite(V).all():raise ValueError('Nonfinite GBIF photo source')
    for j in range(V.shape[1]):
        sol,stop,*_=lsmr(Z,V[:,j],atol=1e-11,btol=1e-11,maxiter=8000)
        if stop==7:raise RuntimeError('FE sparse projection not converged')
        V[:,j]-=Z@sol
    X,Y=V[:,:-1],V[:,-1]
    beta,_,rank,_=np.linalg.lstsq(X,Y,rcond=1e-11)
    if rank!=len(predictors):return {'estimable':False,'reason':'rank_deficient_FE'}
    res=Y-X@beta
    cluster,keys=pd.factorize(a.island,sort=True);G=len(keys)
    score=np.column_stack([np.bincount(cluster,weights=X[:,j]*res,minlength=G)
                           for j in range(len(predictors))])
    bread=np.linalg.inv(X.T@X)
    if two_way:
        ns,ny=dims;s,y=codes;N=ns+ny
        graph=coo_matrix((np.ones(2*n),(np.r_[s,ns+y],np.r_[ns+y,s])),
                         shape=(N,N)).tocsr()
        components=connected_components(graph,directed=False)[0]
        fe_rank=N-components
    else:
        components=1;fe_rank=dims[0]
    df=n-fe_rank-len(predictors)
    if df<=0:raise ValueError('No source FE residual degrees of freedom')
    cov=bread@score.T@score@bread*(G/(G-1))*((n-1)/df)
    se=math.sqrt(max(0,float(cov[0,0])))
    tcrit=float(student_t.ppf(.975,G-1))
    return {'estimable':True,'cell_years':n,'sites':int(a.site.nunique()),
      'independent_islands':G,'islandyears':int(a.islandyear.nunique()),
      'response':outcome+'_photographic_source_share_pp',
      'site_FE':two_way,'weight':weight,
      'slope_light':float(beta[0]),
      'CI95':[float(beta[0]-tcrit*se),float(beta[0]+tcrit*se)],
      'nominal_p':float(2*student_t.sf(abs(beta[0]/se),G-1)) if se else None,
      'FE_connected_components':int(components),'residual_df':df}

def run(insects,reptiles):
    a=read_pair(insects,reptiles)
    spatial=cohort(a)
    if (len(spatial),spatial.island.nunique())!=(3424,50):
        raise ValueError('Frozen co-observed spatial cohort drift')
    repeated=repeated_sites(spatial,.02)
    if (len(repeated),repeated.site.nunique(),repeated.island.nunique())!=(1979,454,30):
        raise ValueError('Frozen actual repeated source cohort drift')
    models={}
    for guild in ('gecko','moth'):
        models[guild]={
          'same_islandyear_spatial_50islands':fit(spatial,guild),
          'same_islandyear_spatial_in_repeat_site_30islands':fit(repeated,guild),
          'same_exact_site_temporal_30islands':fit(repeated,guild,True),
          'same_exact_site_temporal_cell_equal':fit(repeated,guild,True,'cell_equal'),
        }
    return {
      'schema':'ALAN_GBIF_13yr_gecko_moth_parallel_source_photo_light_vs_site_FE_v1',
      'source_SHA_insects':'64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852',
      'source_SHA_gecko_scink':'03aa855c9e7691855bfee4baa099f0d3a31db8d53066ec36a95afda48c3ba339',
      'years':list(range(2013,2026)),
      'actual_same_original_island_cell_prior_light_source_support':15255,
      'spatial_50island_common_source_cell_years':3424,
      'repeated_30island_true_site_years':1979,
      'results':models,
      'limitations':{
        'photo_source_fraction_not_animal_population_abundance':True,
        'GloBI_global_flower_use_not_GBIF_local_pollination':True,
        'zero_guild_photo_record_not_true_species_absence':True,
        'identical_photo_cell_not_proven_trophic_contact':True,
        'original_PCNL_prior_calendar_year_not_ground_lamp_illuminance':True,
        'spatial_full_cohort_50islands_and_temporal_repeat_30islands_differ':True,
        'nominal_exploratory_p_values_unadjusted':True,
        'private_source_species_observer_and_coordinates_not_exported':True}}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--insects',type=Path,required=True)
    p.add_argument('--reptiles',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():raise ValueError('Refuse source result overwrite')
    result=run(a.insects,a.reptiles)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+'\n')
    print('ALAN_GECKO_AND_MOTH_PARALLEL_SOURCE_LIGHT',json.dumps(result['results']['gecko']))
