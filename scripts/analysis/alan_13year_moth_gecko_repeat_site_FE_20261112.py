#!/usr/bin/env python3
"""2013–2025 source-authenticated GBIF photo cell/year repeat-site ALAN analysis.

Independent source-paired photographed insect families and photographed
geckos/skinks are observational evidence, not field predation or pollination.
No 0.05° site identity, coordinates or species×island source rows are exported.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csc_matrix, hstack
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsmr
from scipy.stats import t as student_t

from alan_13year_moth_gecko_same_island_source_ALAN_20261112 import (
    source, SHA_INSECT, SHA_REPTILE, GROUP, cohort
)
PREDICTORS=('loglight','gecko_share','cross','log_insect','log_reptile')

def require(ok, why):
    if not bool(ok): raise ValueError(why)

def read_pair(insects: Path, reptiles: Path):
    i=source(insects,SHA_INSECT,set(GROUP+[
        'moth_family_records','butterfly_family_records','all_four_family_records',
        'prior_log1p_PCNL','PCNL_year','center_inside_same_original_GSHHG_island']))
    r=source(reptiles,SHA_REPTILE,set(GROUP+[
        'gecko_source_photo_count','scink_source_photo_count',
        'original_total_reptile_photo_source_events','log1p_prior_PCNL',
        'prior_PCNL_year','center_inside_same_original_GSHHG_island']))
    require(len(i)==228361 and len(r)==34056,'Source event-year grain drift')
    a=i.merge(r,on=GROUP,how='inner',validate='one_to_one',suffixes=('_insect','_reptile'))
    require(len(a)==17944,'Source island×year×cell overlap changed')
    a=a.loc[a.center_inside_same_original_GSHHG_island_insect.eq(True)&
            a.center_inside_same_original_GSHHG_island_reptile.eq(True)&
            a.prior_log1p_PCNL.notna()&a.log1p_prior_PCNL.notna()].copy()
    require(len(a)==15255,'Original PCNL geography coverage changed')
    require(np.allclose(a.prior_log1p_PCNL,a.log1p_prior_PCNL,atol=1e-10),
            'Independent PCNL source values disagree')
    require((a.PCNL_year.eq(a.year-1)&a.prior_PCNL_year.eq(a.year-1)).all(),
            'PCNL exposure is not before the photo year')
    require(a.all_four_family_records.eq(
            a.moth_family_records+a.butterfly_family_records).all(),
            'Moth/butterfly source count changed')
    require(a.original_total_reptile_photo_source_events.eq(
            a.gecko_source_photo_count+a.scink_source_photo_count).all(),
            'Gecko/skink photo source count changed')
    a['islandyear']=a.island.astype(str)+'|'+a.year.astype(str)
    a['site']=a.island.astype(str)+'|'+a.latbin.astype(str)+'|'+a.lonbin.astype(str)
    a['loglight']=a.prior_log1p_PCNL.astype(float)
    a['moth_share']=100*a.moth_family_records/a.all_four_family_records
    a['gecko_share']=a.gecko_source_photo_count/a.original_total_reptile_photo_source_events
    a['cross']=a.loglight*a.gecko_share
    a['log_insect']=np.log1p(a.all_four_family_records)
    a['log_reptile']=np.log1p(a.original_total_reptile_photo_source_events)
    a['latitude']=a.latitude_center_insect
    return a

def repeated_sites(a,span=.02):
    # Site identity is the SAME physical GSHHG island and 0.05° photo-cell
    years=a.groupby('site').year.nunique()
    x=a.loc[a.site.isin(years.index[years.ge(2)])].copy()
    ranges=x.groupby('site').loglight.agg(lambda z:z.max()-z.min())
    x=x.loc[x.site.isin(ranges.index[ranges.gt(span)])].copy()
    return x

def fit(a,*,weight='island_equal',two_way=True,include_effort=True):
    a=a.copy()
    if len(a)<100 or a.island.nunique()<15:
        return {'estimable':False,'reason':'too_few_independent_island_siteyears'}
    cols=list(PREDICTORS if include_effort else PREDICTORS[:3])
    if not two_way:cols.append('latitude')
    fe_names=('site','islandyear') if two_way else ('islandyear',)
    if weight=='island_equal':
        w=1/a.groupby('island').island.transform('size').to_numpy(float)
    elif weight=='cell_equal':w=np.ones(len(a))
    elif weight=='islandyear_equal':
        w=1/a.groupby('islandyear').islandyear.transform('size').to_numpy(float)
    else:raise ValueError('Unknown source weight')
    require(np.isfinite(w).all() and (w>0).all(),'Nonpositive photograph source weights')
    sw=np.sqrt(w)
    n=len(a)
    variables=a[cols+['moth_share']].to_numpy(float)*sw[:,None]
    require(np.isfinite(variables).all(),'Nonfinite source photo outcome')
    matrices=[];fe_dims={};FE=[]
    for name in fe_names:
        codes,keys=pd.factorize(a[name],sort=True)
        require((codes>=0).all(),'Missing fixed effect')
        matrices.append(csc_matrix((sw,(np.arange(n),codes)),
                        shape=(n,len(keys))))
        fe_dims[name]=len(keys);FE.append(codes)
    Z=hstack(matrices,format='csc')
    v=np.empty_like(variables)
    for j in range(variables.shape[1]):
        coeff,status,*_=lsmr(Z,variables[:,j],atol=1e-11,btol=1e-11,maxiter=8000)
        require(status!=7,'Sparse FE projection did not converge')
        v[:,j]=variables[:,j]-Z@coeff
    X,Y=v[:,:-1],v[:,-1]
    beta,_,rank,_=np.linalg.lstsq(X,Y,rcond=1e-11)
    if rank!=len(cols):return {'estimable':False,'reason':'fixed_effect_collinearity'}
    if two_way:
        ns,ny=fe_dims['site'],fe_dims['islandyear']
        site,iy=FE
        graph=coo_matrix((np.ones(2*n),(
            np.r_[site,ns+iy],np.r_[ns+iy,site])),shape=(ns+ny,ns+ny)).tocsr()
        components,_=connected_components(graph,directed=False)
        rankFE=ns+ny-components
    else:
        components=1;rankFE=sum(fe_dims.values())
    df=n-rankFE-len(cols)
    require(df>0,'Insufficient residual degrees of freedom')
    err=Y-X@beta
    cluster,ck=pd.factorize(a.island,sort=True)
    G=len(ck)
    score=np.column_stack([np.bincount(cluster,
          weights=X[:,j]*err,minlength=G) for j in range(len(cols))])
    bread=np.linalg.inv(X.T@X)
    cov=bread@score.T@score@bread*(G/(G-1))*((n-1)/df)
    tcrit=float(student_t.ppf(.975,G-1))
    effects={}
    for j,key in enumerate(cols):
        se=math.sqrt(max(0,float(cov[j,j])))
        effects[key]={'beta':float(beta[j]),'SE_cluster_island':se,
                      'CI95':[float(beta[j]-tcrit*se),float(beta[j]+tcrit*se)],
                      'nominal_p':float(2*student_t.sf(abs(beta[j]/se),G-1)) if se else None}
    return {'estimable':True,'cell_years':n,'sites':int(a.site.nunique()),
      'independent_islands':G,'island_years':int(a.islandyear.nunique()),
      'source_photo_years':int(a.year.nunique()),
      'fixed_effects':fe_dims,'source_connected_FE_components':int(components),
      'effective_residual_df':df,'weighting':weight,'coefficients':effects}

def counts(a):
    return {'cell_years':int(len(a)),'sites':int(a.site.nunique()),
            'independent_islands':int(a.island.nunique()),
            'island_years':int(a.islandyear.nunique()),
            'years':sorted(map(int,a.year.unique()))}

def run(insects,reptiles):
    a=read_pair(insects,reptiles)
    primary=cohort(a)
    require((len(primary),primary.island.nunique(),
            primary.islandyear.nunique())==(3424,50,209),
            'Original spatial 13-year model not exactly recovered')
    x=repeated_sites(primary,.02)
    require((len(x),x.site.nunique(),x.island.nunique())==(1979,454,30),
            'True repeated same-site source cohort changed')
    x0=repeated_sites(primary,0)
    high=repeated_sites(primary,.1)
    results={
      'spatial_on_repeated_cohort':fit(x,two_way=False),
      'primary_same_site_and_islandyear_FE':fit(x),
      'same_site_cell_equal':fit(x,weight='cell_equal'),
      'same_site_islandyear_equal':fit(x,weight='islandyear_equal'),
      'same_site_any_light_variation':fit(x0),
      'same_site_0p1plus_light_variation':fit(high),
      'same_site_without_photo_effort_covariates':fit(x,include_effort=False)
    }
    # Photo opportunity gates apply BEFORE fitting, and still require repeat
    # years of the SAME source cell after restriction.
    strict={}
    for name,column,mincount in (
      ('min3_reptile_photos','original_total_reptile_photo_source_events',3),
      ('min20_insect_photos','all_four_family_records',20)):
        xx=x.loc[x[column].ge(mincount)].copy()
        n=xx.groupby('site').year.nunique()
        xx=xx.loc[xx.site.isin(n.index[n.ge(2)])].copy()
        strict[name]={'support':counts(xx),'fit':fit(xx)}
    # Post-hoc descriptive influence check, not 30 independent tests.
    leave=[]
    for island in sorted(x.island.unique()):
        z=fit(x.loc[x.island.ne(island)])
        if z['estimable']:
            leave.append(z['coefficients']['cross']['beta'])
    report={
      'schema':'alan_true_island_same_grid_multi_year_gecko_moth_priorPCNL_FE_v1',
      'run_method':'direct SHA256 frozen full original GBIF two-arm 2013–2025 photo source replay',
      'original_source_SHA256':{'four_insect_families':SHA_INSECT,'geckos_and_skinks':SHA_REPTILE},
      'site_definition':'GSHHG island ID plus 0.05-degree grid; >=2 photo years and prior-light span >0.02',
      'observation_years':list(range(2013,2026)),
      'source_support':{
        'actual_same_island_cell_and_year_insect_reptile_source_matches':17944,
        'with_valid_equal_past_year_light':15255,
        'original_spatial_cohort':counts(primary),
        'same_site_any_light_variation':counts(x0),
        'same_site_prior_loglight_span_gt0p02':counts(x),
        'same_site_prior_loglight_span_gt0p1':counts(high),
        'observed_reptile_photo_fraction_exact_zero':float(x.gecko_share.eq(0).mean()),
        'observed_reptile_photo_fraction_exact_one':float(x.gecko_share.eq(1).mean()),
        'median_reptile_photo_opportunities':float(x.original_total_reptile_photo_source_events.median()),
        'median_repeated_site_loglight_range':float(x.groupby('site').loglight.agg(lambda v:v.max()-v.min()).median())
      },
      'same_site_models':results,
      'strict_photo_opportunity_gates':strict,
      'leave_one_island_out_descriptive':{
        'n_overlapping_omissions':len(leave),
        'negative_beta_omissions':sum(v<0 for v in leave),
        'beta_range':[float(min(leave)),float(max(leave))],
        'not_independent_replications':True
      },
      'interpretation_boundaries':{
        'photos_not_individual_animal_population_abundance':True,
        'no_photo_record_NOT_true_moth_or_gecko_absence':True,
        'no_actual_gecko_moth_predation_observed':True,
        'no_moth_pollen_carriage_or_effective_pollination_observed':True,
        'source_photo_observation_effort_not_standardized':True,
        'prior_corrected_PCNL_not_local_lamp_illuminance':True,
        'statistical_interaction_not_mechanistic_predation':True,
        'original_positive_README_PL64_result_untouched':True,
        'no_private_exact_plant_photo_sites_or_observation_rows_exported':True
      }}
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--insects',type=Path,required=True)
    ap.add_argument('--reptiles',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():raise ValueError('Do not overwrite prior result')
    r=run(args.insects,args.reptiles)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('ALAN_TRUE_REPEATED_SOURCE_ISLAND_MOTH_GECKO',json.dumps(
      r['same_site_models']['primary_same_site_and_islandyear_FE'],sort_keys=True))
