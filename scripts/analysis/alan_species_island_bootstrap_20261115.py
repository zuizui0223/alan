#!/usr/bin/env python3
"""Original 703-species photographic ALAN first-stage island block bootstrap.

Re-fit first-stage source sufficient statistics when whole ISLANDS (not plant
species) are resampled. A species appears in many of the same photographed
cells as other species, so inverse-variance independent-effect meta-regression
does NOT supply valid overall uncertainty here.

This is an observational *photo nonrecord* analysis without population/climate
covariates; bootstrap inference is conditional on that incomplete model.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    load_sites,load_labels,original_accepted_species,require
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits
from scripts.analysis.alan_species_effortadjusted_firststage_20261115 import (
    prepare_support,site_arrays,fit_one,island_fold
)

MIN_ISLANDS=30
SEED=20261115
N_BOOT=600

def island_sufficient_stats(site,site_ids,strata,universe):
    s=site.loc[site.iy.isin(strata)].copy()
    sid=s.sid.to_numpy(np.int64)
    y=np.isin(sid,site_ids,assume_unique=True).astype(float)
    yc=y-s.assign(y=y).groupby("iy",observed=True).y.transform("mean").to_numpy(float)
    x=s.xc.to_numpy(float)
    e=s.ec.to_numpy(float)
    code=pd.Categorical(s.island,categories=universe).codes
    require((code>=0).all(),"Unknown source island")
    g=len(universe)
    n=np.bincount(code,minlength=g)
    weight=1/n[code]
    def collect(a):
        return np.bincount(code,weights=a,minlength=g)
    a=np.column_stack([
        collect(weight*x*x),
        collect(weight*x*e),
        collect(weight*e*e),
    ])
    b=np.column_stack([
        collect(weight*x*yc),
        collect(weight*e*yc),
    ])
    return a,b,int((n>0).sum())

def solve_multispecies(A,B,mult):
    """Same two-regressor source original slope; multiply island score blocks."""
    aa=np.einsum("sik,i->sk",A,mult,optimize=True)
    bb=np.einsum("sik,i->sk",B,mult,optimize=True)
    det=aa[:,0]*aa[:,2]-aa[:,1]**2
    valid=np.isfinite(det)&(det>1.e-13)
    val=np.full(len(det),np.nan)
    val[valid]=100*(bb[valid,0]*aa[valid,2]-bb[valid,1]*aa[valid,1])/det[valid]
    return val

def group_contrast(beta,aut):
    is_auto=aut.astype(bool)
    a=beta[is_auto];b=beta[~is_auto]
    aa=a[np.isfinite(a)];bb=b[np.isfinite(b)]
    if len(aa)<max(15,int(.90*len(a))) or len(bb)<max(15,int(.90*len(b))):
        return np.nan
    return float(aa.mean()-bb.mean())

def bootstrap_islands(A,B,autonomy,high_quality,folds,n_boot=N_BOOT,seed=SEED):
    n,G,_=A.shape
    require(B.shape==(n,G,2) and len(autonomy)==n and len(high_quality)==n,
            "Inconsistent source panel")
    orig=solve_multispecies(A,B,np.ones(G))
    original=group_contrast(orig,autonomy)
    m0=(folds==0).astype(float);m1=(folds==1).astype(float)
    test0=solve_multispecies(A,B,m0)
    test1=solve_multispecies(A,B,m1)
    fold0=group_contrast(test0,autonomy)
    fold1=group_contrast(test1,autonomy)
    rng=np.random.default_rng(seed)
    contrast=np.full(n_boot,np.nan)
    contrast_high=np.full(n_boot,np.nan)
    for j in range(n_boot):
        counts=np.bincount(rng.integers(0,G,G),minlength=G).astype(float)
        bet=solve_multispecies(A,B,counts)
        contrast[j]=group_contrast(bet,autonomy)
        contrast_high[j]=group_contrast(bet[high_quality],autonomy[high_quality])
    valid=contrast[np.isfinite(contrast)]
    high=contrast_high[np.isfinite(contrast_high)]
    require(len(valid)>=.90*n_boot,"Too few island bootstrap replicates qualified")
    return {
      "n_species_in_fixed_30island_cohort":n,
      "n_islands_in_common_gradient_universe":G,
      "autonomous_species":int(autonomy.sum()),
      "absent_species":int((~autonomy).sum()),
      "high_only_species":int(high_quality.sum()),
      "original_ge30_island_autonomous_minus_absent_pp":original,
      "fold0_same_species_cohort_diff_pp":fold0,
      "fold1_same_species_cohort_diff_pp":fold1,
      "bootstrap_unit":"entire genuine GSHHG island from the common 243-island source",
      "bootstrap_method":"Refit both species-stage light and effort coefficients from weighted original island score blocks, then compare equal-species-mean slopes",
      "n_draws_requested":n_boot,
      "n_valid_draws":int(len(valid)),
      "bootstrap_diff_95percentile_ci_pp":[float(np.quantile(valid,.025)),float(np.quantile(valid,.975))],
      "bootstrap_diff_median_pp":float(np.median(valid)),
      "bootstrap_draw_negative_fraction":float((valid<0).mean()),
      "high_quality_original_diff_pp":group_contrast(orig[high_quality],autonomy[high_quality]),
      "high_quality_bootstrap_n_valid":int(len(high)),
      "high_quality_bootstrap_diff_95percentile_ci_pp":(
           [float(np.quantile(high,.025)),float(np.quantile(high,.975))]
           if len(high) else None),
      "island_bootstrap_not_independent_geographic_region_or_climate_control":True
    }

def run(sites_path,globi_path,groups,ledger,out):
    require(not out.exists(),"Refuse source result overwrite")
    sites=load_sites(sites_path)
    flora=original_accepted_species(groups,sites,load_labels(globi_path))
    ledger_hm=direct_plant_traits(ledger)
    ledger_h=direct_plant_traits(ledger,high_only=True)
    taxa=pd.DataFrame({"species":flora.species.drop_duplicates()})
    taxa["autonomy"]=taxa.species.map(ledger_hm.autonomous_selfing_capacity)
    taxa["high_quality"]=taxa.species.map(ledger_h.autonomous_selfing_capacity).isin(
        ("autonomous","absent"))
    site,positive,groupspec,selected,coverage=prepare_support(
        sites,flora,taxa,span=.02)
    require((coverage["eligible_islands"],coverage["strict_selfing_species"])==
            (243,703),"Original 703 species support drift")
    site,source=site_arrays(site)
    universe=np.sort(site.island.astype(str).unique())
    folded=np.array([island_fold(x) for x in universe],dtype=int)
    trait=taxa.loc[taxa.species.isin(selected)].set_index("species")
    pos=dict((s,g.sid.to_numpy(np.int64))
             for s,g in positive.groupby("species",observed=True))
    inform=dict((s,g.iy.to_numpy(int))
             for s,g in groupspec.groupby("species",observed=True))
    A=[];B=[];status=[];aut=[];hq=[]
    for j,sp in enumerate(sorted(selected)):
        aa,bb,ncl=island_sufficient_stats(site,pos[sp],inform[sp],universe)
        if ncl<MIN_ISLANDS:continue
        base=fit_one(site,inform[sp],pos[sp])
        require(base["ok"] and base["n_islands"]==ncl,
                "Within-stratum original model differs")
        beta=solve_multispecies(aa[None,...],bb[None,...],np.ones(len(universe)))[0]
        require(np.isfinite(beta) and abs(beta-base["beta_pp"])<1.e-7,
                "Bootstrap sufficient-statistic slope differs from source main")
        A.append(aa);B.append(bb)
        status.append(ncl)
        aut.append(str(trait.loc[sp,"autonomy"])=="autonomous")
        hq.append(bool(trait.loc[sp,"high_quality"]))
        if len(A)%50==0:print("ORIGINAL_ISLAND_CLUSTER_STATS",len(A),flush=True)
    require(len(A)==224,"Frozen independent-island support drift; stop")
    A=np.stack(A)
    B=np.stack(B)
    result=bootstrap_islands(A,B,np.array(aut,dtype=bool),
                             np.array(hq,dtype=bool),folded)
    require(abs(result["original_ge30_island_autonomous_minus_absent_pp"]-
                (0.6600646913763487-0.6714059365470156))<1.e-7,
            "Prior original-source 30+ island result cannot be reconstructed")
    report={
      "schema":"ALAN_autonomous_selfing_30island_original_species_slopes_island_bootstrap_v1",
      "status":"REAL_ORIGINAL_SOURCES_ISLAND_BLOCK_REPLICATION_EXPLORATORY",
      "source":{"n_original_islands":562,"n_original_site_years":77647,
                "eligible_gradient_islands":243,"eligible_photo_site_years":72472,
                "strict_selfing_matched":703,
                "source_sites_and_botanical_traits_are_SHA_VERIFIED":True},
      "block_bootstrap":result,
      "limits":{
        "photographic_nonrecord_not_true_absence":True,
        "no_independent_population_and_climate_controls":True,
        "species_with_30plus_islands_are_geographic_prevalence_selected":True,
        "no_full_703_species_meta_effect_claim":True,
        "within_island_year_observation_effort_only":True,
        "shared_island_resampling_does_not_resolve_spatial_neighbor_dependence":True,
        "no_true_pollen_deposition_or_reproductive_output":True,
        "no_species_names_or_private_source_location_exported":True
      }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_SPECIES_ISLAND_BLOCK_BOOTSTRAP",json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("sites","globi","groups","ledger","out"):
        p.add_argument("--"+a,required=True,type=Path)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.ledger,a.out)
