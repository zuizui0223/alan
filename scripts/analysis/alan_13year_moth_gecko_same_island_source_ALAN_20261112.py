#!/usr/bin/env python3
"""Actual 2013–2025 GBIF moth-butterfly / gecko-skink same-island photos × prior PCNL.

Positive documented four-family insect PHOTO observations and gecko/skink
PHOTO observations; no inference of observed predation, pollination,
population density or biological absence from unrecorded cases.
No individual observation or 0.05° source rows are exported.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import t as student_t

SHA_INSECT="64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852"
SHA_REPTILE="03aa855c9e7691855bfee4baa099f0d3a31db8d53066ec36a95afda48c3ba339"
GROUP=['island','year','latbin','lonbin']
def require(x,reason):
    if not x: raise ValueError(reason)

def source(path,sha,required):
    require(hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha,"Source CSV.gz SHA mismatch")
    d=pd.read_csv(path,low_memory=False)
    require(required<=set(d),f"Source schema drift {required-set(d)}")
    require(not d.duplicated(GROUP).any(),"Original island/year/cell key duplication")
    return d

def cohort(s,insect_min=5,reptile_min=1):
    d=s.loc[s.all_four_family_records.ge(insect_min) &
      s.moth_family_records.gt(0)&s.butterfly_family_records.gt(0)&
      s.original_total_reptile_photo_source_events.ge(reptile_min)].copy()
    counts=d.groupby('islandyear')[['gecko_source_photo_count','scink_source_photo_count']].sum()
    both=counts.index[counts.gecko_source_photo_count.gt(0)&counts.scink_source_photo_count.gt(0)]
    d=d.loc[d.islandyear.isin(both)].copy()
    g=d.groupby('islandyear').loglight.agg(['size','min','max'])
    supported=g.index[g['size'].ge(2)&(g['max']-g['min']).ge(.02)]
    return d.loc[d.islandyear.isin(supported)].copy()

def fit(d,weight='island_equal'):
    require(len(d)>=100 and d.island.nunique()>=20,"Source common support below floor")
    a=d.copy()
    a['share']=100*a.moth_family_records/a.all_four_family_records
    a['gecko_fraction']=a.gecko_source_photo_count/a.original_total_reptile_photo_source_events
    a['cross']=a.loglight*a.gecko_fraction
    a['log_insect']=np.log1p(a.all_four_family_records)
    a['log_reptile']=np.log1p(a.original_total_reptile_photo_source_events)
    a['latitude']=a.latitude_center_insect-a.latitude_center_insect.mean()
    cols=['loglight','gecko_fraction','cross','log_insect','log_reptile','latitude']
    r,g=pd.factorize(a.islandyear,sort=True)
    c,clusters=pd.factorize(a.island,sort=True)
    if weight=='island_equal':
        w=1/a.groupby('island').island.transform('size').to_numpy(float)
    elif weight=='cell_equal':
        w=np.ones(len(a))
    elif weight=='islandyear_equal':
        w=1/a.groupby('islandyear').islandyear.transform('size').to_numpy(float)
    else:raise ValueError(weight)
    z=a[cols+['share']].to_numpy(float)
    require(np.isfinite(z).all() and np.isfinite(w).all(),"Nonfinite original source")
    totals=np.bincount(r,weights=w,minlength=len(g))
    for j in range(z.shape[1]):
        means=np.bincount(r,weights=w*z[:,j],minlength=len(g))/totals
        z[:,j]=(z[:,j]-means[r])*np.sqrt(w)
    X,y=z[:,:-1],z[:,-1]
    beta,_,rank,_=np.linalg.lstsq(X,y,rcond=1e-11)
    require(rank==len(cols),"FE design not identifiable")
    residual=y-X@beta
    scores=np.column_stack([np.bincount(c,weights=X[:,j]*residual,minlength=len(clusters))
                            for j in range(len(cols))])
    bread=np.linalg.inv(X.T@X)
    G=len(clusters)
    V=bread@scores.T@scores@bread*(G/(G-1))*((len(a)-1)/(len(a)-len(cols)))
    crit=float(student_t.ppf(.975,G-1))
    effect={}
    for j,k in enumerate(cols):
        se=math.sqrt(max(0,V[j,j]))
        effect[k]={'beta':float(beta[j]),
          'ci95':[float(beta[j]-crit*se),float(beta[j]+crit*se)],
          'p_nominal':float(2*student_t.sf(abs(beta[j]/se),G-1)) if se else None}
    return {'photo_cell_years':len(a),'independent_islands':G,
       'island_year_fixed_effects':len(g),
       'source_moth_and_butterfly_photo_records':int(a.all_four_family_records.sum()),
       'source_gecko_and_skink_photo_records':int(a.original_total_reptile_photo_source_events.sum()),
       'weighting':weight,'coefficient':effect}

def run(insect_file,reptile_file):
    i=source(insect_file,SHA_INSECT,set(GROUP+['moth_family_records','butterfly_family_records','all_four_family_records','prior_log1p_PCNL','PCNL_year','center_inside_same_original_GSHHG_island']))
    r=source(reptile_file,SHA_REPTILE,set(GROUP+['gecko_source_photo_count','scink_source_photo_count','original_total_reptile_photo_source_events','log1p_prior_PCNL','prior_PCNL_year','center_inside_same_original_GSHHG_island']))
    require(len(i)==228361 and len(r)==34056,"Source row population drift")
    d=i.merge(r,on=GROUP,how='inner',validate='one_to_one',suffixes=('_insect','_reptile'))
    require(len(d)==17944,"Original co-observation identities changed")
    d=d.loc[d.center_inside_same_original_GSHHG_island_insect.eq(True) &
            d.center_inside_same_original_GSHHG_island_reptile.eq(True) &
            d.prior_log1p_PCNL.notna() & d.log1p_prior_PCNL.notna()].copy()
    require(len(d)==15255,"Valid same-island source population drift")
    require(np.allclose(d.prior_log1p_PCNL,d.log1p_prior_PCNL,atol=1e-10),
        "Two original PCNL source exposures disagree at the same cell/year")
    require(d.PCNL_year.eq(d.year-1).all() and d.prior_PCNL_year.eq(d.year-1).all(),
        "Source uses future or non-prior satellite exposure")
    require(d.all_four_family_records.eq(d.moth_family_records+d.butterfly_family_records).all(),
        "Insect source counts are not conserved")
    require(d.original_total_reptile_photo_source_events.eq(
        d.gecko_source_photo_count+d.scink_source_photo_count).all(),
        "Reptile source counts are not conserved")
    d['islandyear']=d.island.astype(str)+'|'+d.year.astype(str)
    d['loglight']=d.prior_log1p_PCNL
    main=cohort(d)
    strict=cohort(d,reptile_min=3)
    high=cohort(d,insect_min=20)
    res={'primary_island_equal':fit(main),'strict_reptile_3plus':fit(strict),
         'strict_insect_20plus':fit(high),
         'primary_cell_equal':fit(main,'cell_equal'),
         'primary_islandyear_equal':fit(main,'islandyear_equal')}
    return {'schema':'alan_true_GBIF_moth_gecko_light_13years_v1',
      'source_insect_sha256':SHA_INSECT,'source_reptile_sha256':SHA_REPTILE,
      'source_years':list(range(2013,2026)),
      'real_exact_geo_year_cell_matches':17944,
      'source_valid_same_cell_prior_PCNL_matches':15255,
      'source_positive_only_insect_and_lizard_guilds':True,
      'model':'Within actual GSHHG island×year fixed effects, island-cluster CR1 CI; response documented moth-photo fraction of moth+butterfly photo observations, predictor gecko-photo fraction of gecko+skink observations × prior-year log1p PCNL; controls site photo opportunities, latitude; >=2 cells and >=.02 year-prior light span; true same 0.05 degree island/year source cell IDs',
      'model_estimates':res,
      'biological_limits':{
       'gecko_photo_fraction_not_density_or_predation':True,
       'moth_photo_fraction_not_pollen_delivery_or_population':True,
       'unrecorded_GBIF_species_not_species_absence':True,
       'same_cell_not_same_flower_or_predator_prey_event':True,
       'observational_photo_sampling_not_independent_trap_effort':True,
       'historical_prior_PCNL_not_local_lamp_lux':True,
       'no_original_sensitive_site_coordinates_exported':True}}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--insects',type=Path,required=True)
    parser.add_argument('--reptiles',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    if a.output.exists():raise ValueError("Do not overwrite previous source result")
    out=run(a.insects,a.reptiles)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print('ACTUAL_GBIF_13YEAR_CROSSGUILD',json.dumps(out['model_estimates']['primary_island_equal']))
