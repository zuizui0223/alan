#!/usr/bin/env python3
"""Same 0.05-degree reptile PHOTO cell repeated over original 2013-2025 years.

A direct sensitivity to fixed urban/site observation preference for gecko
versus Scincidae source records. Do NOT assume that within-cell variation
in GBIF photo share reflects population recruitment or preference: source
gecko/skink zero is only no photo of that category, and no observer/time
effort control is available.

Historical PCNL is always observation year-1 (2012-24). Two-way
site identity + calendar-year weighted fixed effects, real island-cluster
uncertainty when original site/year/contrast support exists; otherwise
fail closed and report exactly why.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
from scripts.analysis.alan_gecko_skink_recorded_light_niche_FE_20261110 import (
    load_original,support_coobserved)
SCHEMA="ALAN_gecko_vs_scink_original_repeated_source_gridcells_temporal_PCNl_FE_v1"

def require(ok,why):
    if not bool(ok):raise ValueError("ALAN_GECKO_REPEATED_SOURCE: "+str(why))

def time_panel(d,min_per_site_years=2,min_light_span=.02):
    x,source_support=support_coobserved(d,min_each=1,min_cell_photo_events=1)
    x["site_id"]=(x.island.astype(str)+"|"+
                  x.latbin.astype(str)+"|"+x.lonbin.astype(str))
    site=x.groupby("site_id",sort=True).agg(
        source_observed_years=("year","nunique"),
        source_real_observation_rows=("year","size"),
        PCNL_min=("log1p_prior_PCNL","min"),
        PCNL_max=("log1p_prior_PCNL","max"),
        source_real_gecko_photo_count=("gecko_source_photo_count","sum"),
        source_real_skink_photo_count=("scink_source_photo_count","sum"))
    good=site.index[site.source_observed_years.ge(min_per_site_years)&
               (site.PCNL_max-site.PCNL_min).ge(min_light_span)&
               site.source_real_gecko_photo_count.ge(1)&
               site.source_real_skink_photo_count.ge(1)]
    selected=x.loc[x.site_id.isin(good)].copy()
    support={"original_full_gecko_skink_source":source_support,
        "original_repeated_same_grid_sites_2plus_years_with_within_site_PCNL_span":len(good),
        "real_repeated_siteyear_source_rows":len(selected),
        "real_repeated_temporal_islands":int(selected.island.nunique()),
        "real_years_observed":sorted(selected.year.unique().astype(int).tolist()),
        "observed_years_per_site_quantiles":(
            {str(v):float(site.loc[good,"source_observed_years"].quantile(v))
             for v in (.1,.5,.9)} if len(good) else {}),
        "year_minus_one_pixel_differed_by_at_least":min_light_span,
        "repeated_sites_have_at_least_one_photo_of_each_reptile_category_over_all_observed_years":True,
        "original_photo_grid_not_actual_identically_surveyed_road_or_streetlamp":True}
    return selected,support

def repeated_weighted_demean(A,site,year,weights,max_iter=400,tol=2e-10):
    z=np.asarray(A,dtype=float).copy()
    codes=[pd.factorize(site,sort=True)[0],
           pd.factorize(year,sort=True)[0]]
    for iteration in range(max_iter):
        prior=z.copy()
        for c in codes:
            denom=np.bincount(c,weights=weights)
            require((denom>0).all(),"Unsupported original site or year FE")
            means=np.vstack([np.bincount(c,weights=weights*z[:,j])/denom
                             for j in range(z.shape[1])]).T
            z-=means[c]
        if np.max(np.abs(z-prior)) < tol:
            break
    else:
        return None,None
    means_after=[]
    for c in codes:
        den=np.bincount(c,weights=weights)
        means_after.append(max(float(np.max(np.abs(np.bincount(
            c,weights=weights*z[:,j])/den))) for j in range(z.shape[1])))
    require(max(means_after)<1e-7,"Original site and year effects not fully removed")
    return z,iteration+1

def fit_panel(data,support,mode="island_equal"):
    n=len(data);g=int(data.island.nunique())
    years=int(data.year.nunique());sites=int(data.site_id.nunique())
    if n<120 or g<25 or sites<60 or years<5 or n<sites+years+12:
        return {"status":"unidentifiable_repeated_cell_year_gecko_observation_light",
                "source_support":support,"N":n,"independent_islands":g,
                "source_photo_sites":sites,"observation_years":years}
    x=data.copy()
    x["light"]=x.log1p_prior_PCNL
    x["n_events"]=x.log_original_reptile_record_count
    xx=x[["light","n_events"]].to_numpy(float)
    yy=x.gecko_photo_record_percent.to_numpy(float)
    w= (1/x.groupby("island").island.transform("size").to_numpy(float)
        if mode=="island_equal" else np.ones(n,float))
    require(np.isfinite(w).all() and w.min()>0,"Original cluster source weights invalid")
    all_data=np.column_stack([yy,xx])
    demean,iterations=repeated_weighted_demean(
        all_data,x.site_id,x.year,w)
    if demean is None:
        return {"status":"original_unbalanced_siteyear_FE_not_converged",
                "source_support":support}
    X=demean[:,1:];y=demean[:,0]
    k=X.shape[1]
    bread_inv=(X*w[:,None]).T@X
    rank=int(np.linalg.matrix_rank(bread_inv))
    if rank<k or np.linalg.cond(bread_inv)>1e9:
        return {"status":"near_collinear_within_site_light_and_photo_effort",
                "source_support":support,"source_fe_matrix_rank":rank}
    bread=np.linalg.inv(bread_inv)
    beta=bread@((X*w[:,None]).T@y)
    res=y-X@beta
    groups=pd.factorize(x.island,sort=True)[0]
    scores=np.zeros((g,k))
    np.add.at(scores,groups,X*(w*res)[:,None])
    cov=bread@(scores.T@scores)@bread
    df=n-sites-years+1-k
    require(df>0,"Not enough true source within-site and year information")
    cov*=g/(g-1)*(n-1)/df
    se=np.sqrt(np.maximum(0,np.diag(cov)))
    cutoff=float(t.ppf(.975,g-1))
    estimates={}
    for j,name in enumerate(("source_prior_light","source_photo_event_count")):
        estimates[name]={"beta":float(beta[j]),
           "original_independent_island_cluster_SE":float(se[j]),
           "original_island_cluster_nominal_95CI":
             [float(beta[j]-cutoff*se[j]),float(beta[j]+cutoff*se[j])],
           "nominal_exploratory_two_sided_p":
              float(2*t.sf(abs(beta[j]/se[j]),g-1)) if se[j]>0 else None}
    return {"status":"estimated_observed_source_photo_composition_only",
            "source_independent_islands":g,"source_repeated_grid_sites":sites,
            "source_site_year_photos":n,"year_fixed_effects":years,
            "source_site_fixed_effects":sites,
            "fit_2way_FE_projection_iterations":iterations,
            "source_island_equal_weight":mode=="island_equal",
            "source_estimates":estimates,
            "source_photo_not_population_recruitment_or_moth_predation":True,
            "year_minus_one_satellite_not_local_ground_lux":True}

def run(source,receipt,out):
    data,_=load_original(source,receipt)
    x,s=time_panel(data)
    result={"schema":SCHEMA,
       "source_original_13year_GBIF_2013_2025":"10.15468/dl.8yrbe7",
       "source_original_PCNl_2012_2024":"10.5281/zenodo.17013414",
       "true_source_repeated_site_year_support":s,
       "island_equal_temporal_site_and_year_FE":fit_panel(x,s,"island_equal"),
       "siteyear_equal_temporal_site_and_year_FE":fit_panel(x,s,"site_equal"),
       "source_photo_record_share_not_reptile_true_occupancy":True,
       "no_actual_gecko_consumed_moth_or_pollen_effect":True,
       "not_a_causal_longitudinal_gecko_light_treatment":True}
    out=Path(out)
    require(not out.exists(),"No overwriting real 13year repeat-source model")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_GECKO_SCINK_GENUINE_REPEATED_SITE_LIGHT",
          json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("source","receipt","out"):
        p.add_argument("--"+key,required=True,type=Path)
    a=p.parse_args();run(a.source,a.receipt,a.out)
