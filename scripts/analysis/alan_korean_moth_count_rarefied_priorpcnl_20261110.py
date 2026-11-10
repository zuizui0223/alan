#!/usr/bin/env python3
"""Original Korean island moth-trap count and rarefied-richness vs prior PCNL.

Reanalyse independently published 2007-2025 Korean island light-trap sample
abundances (377 *survey site-year* samples on 30 islands; 1179 moth species)
with source-MD5-verified PCNL strictly year-minus-one at the ACTUAL trap site.
Not insect-density, light-source choice, nocturnal floral visitation, or
functional pollination. Light-trap catchability is inherently light dependent.
No trapping-night effort available in the original release.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import t
import statsmodels.formula.api as smf

from scripts.acquisition.alan_korean_island_moth_real_source_preflight_20261110 import (
    SOURCES,get_original)
from scripts.analysis.alan_korean_moth_prior_year_pcnl_17year_20261110 import (
    SCHEMA as SOURCE_SCHEMA,require,checksum)

SCHEMA="ALAN_real_30island_moth_counts_rarefaction_v_past_PCNl_observational_v1"

def expected_rarefaction(counts,m=20):
    """Expected species in m individuals; exact finite-population hypergeometric."""
    x=np.asarray(counts,dtype=float)
    if not (np.isfinite(x).all() and (x>=0).all() and
            np.all(x==np.floor(x))):
        raise ValueError("Original species abundance not nonnegative whole counts")
    N=int(x.sum())
    if N<m:return np.nan
    observed=x[x>0]
    available=N-observed
    miss=np.zeros(len(observed),dtype=float)
    ok=available>=m
    # log C(N-n_i,m)/C(N,m); m! cancels, preserving numerical stability.
    miss[ok]=np.exp(gammaln(available[ok]+1)-
                    gammaln(available[ok]-m+1)-
                    gammaln(N+1)+gammaln(N-m+1))
    return float(np.sum(1-np.clip(miss,0,1)))

def source_frame(pcnl_table,pcnl_receipt):
    pcnl_info=json.loads(Path(pcnl_receipt).read_text())
    require(pcnl_info["schema"]==SOURCE_SCHEMA and
            pcnl_info["original_survey_site_years"]==377 and
            pcnl_info["original_distinct_islands"]==30 and
            checksum(pcnl_table)==pcnl_info["joined_full_source_CSV_SHA256"],
            "Original 17-year yearminus1 source was not SHA verified")
    light=pd.read_csv(pcnl_table)
    need={"island","sites","year","raster_year","prior_year_log1p_PCNL",
          "prior_year_PCNL_native","Elevation"}
    require(need.issubset(light.columns) and
            (light.raster_year==light.year-1).all(),
            "Historical year leakage or lost site-year source")
    raw=pd.read_excel(io.BytesIO(get_original(*SOURCES["occurrences"])))
    species=[c for c in raw.columns if str(c).startswith("sp")]
    require(len(raw)==377 and len(species)==1179 and
            not raw.duplicated(["island","year","sites"]).any(),
            "Original moth matrix changed")
    counts=raw[species].apply(pd.to_numeric,errors="coerce")
    require(counts.notna().all().all() and (counts>=0).all().all() and
            np.all(counts.to_numpy()==np.floor(counts.to_numpy())),
            "Original moth counts cannot be invented/imputed")
    x=raw[["island","year","sites"]].copy()
    x["original_moths_captured"]=counts.sum(axis=1).astype(int)
    x["observed_moth_species"]=counts.gt(0).sum(axis=1).astype(int)
    x["rarefied20_expected_species"]=[
       expected_rarefaction(row,m=20) for row in counts.to_numpy()]
    x["rarefied10_expected_species"]=[
       expected_rarefaction(row,m=10) for row in counts.to_numpy()]
    hill=[]
    for row in counts.to_numpy():
        values=row[row>0].astype(float);N=values.sum()
        if N<=0:hill.append(np.nan)
        else:
            proportions=values/N
            hill.append(float(np.exp(-np.sum(proportions*np.log(proportions)))))
    x["hill1_effective_species"]=hill
    require(int(x.original_moths_captured.sum())==30442 and
            x.island.nunique()==30,"Original 30,442 real moth individuals not reproduced")
    d=x.merge(light[["island","sites","year","raster_year",
                     "prior_year_log1p_PCNL","prior_year_PCNL_native",
                     "Elevation"]],
              on=["island","sites","year"],how="left",validate="one_to_one",indicator=True)
    require(d._merge.eq("both").all() and len(d)==377,
            "Moth source site-year not actually linked to original PCNL cells")
    islands=pd.read_excel(io.BytesIO(get_original(*SOURCES["islands"])))
    require(len(islands)==30 and islands.island.is_unique and
            {"area","distance","elevation"}.issubset(islands.columns),
            "No author-documented 30 source islands environmental data")
    d=d.merge(islands,on="island",how="left",validate="many_to_one")
    require(d.area.gt(0).all() and d.distance.ge(0).all(),
            "Missing original source island area and isolation")
    d["log1p_moths_captured"]=np.log1p(d.original_moths_captured)
    d["log1p_species_richness"]=np.log1p(d.observed_moth_species)
    d["elevation_km"]=pd.to_numeric(d.Elevation,errors="coerce")/1000.
    d["island_year"]=d.island.astype(str)+"|"+d.year.astype(int).astype(str)
    return d,pcnl_info

def valid_within_island_year(d,outcome,light_span_min=.02):
    z=d.loc[d.prior_year_log1p_PCNL.notna() &
            d.elevation_km.notna() & d[outcome].notna()].copy()
    stats=z.groupby("island_year").agg(
        distinct_sites=("sites","nunique"),
        n_different_light=("prior_year_log1p_PCNL","nunique"),
        light_span=("prior_year_log1p_PCNL",lambda t:float(t.max()-t.min())))
    good=stats.index[(stats.distinct_sites>=2)&
                     (stats.n_different_light>=2)&
                     (stats.light_span>=light_span_min)]
    x=z.loc[z.island_year.isin(good)].copy()
    return x,{"source_rows_original":len(d),
              "true_prior_PCNL_measured":int(d.prior_year_log1p_PCNL.notna().sum()),
              "eligible_outcome_rows":len(z),
              "within_island_year_contrast_groups":len(good),
              "within_island_year_contrast_islands":int(x.island.nunique()),
              "within_island_year_contrast_site_years":len(x),
              "fixed_minimum_spatial_loglight_span":light_span_min,
              "actual_moth_survey_years":sorted(x.year.unique().astype(int).tolist())}

def fit(d, outcome,weight="island_equal",covariate_elevation=True):
    x,support=valid_within_island_year(d,outcome)
    G=int(x.island.nunique())
    if G<10 or len(x)<max(60,support["within_island_year_contrast_groups"]+8):
        return {"status":"not_estimated_insufficient_real_islands_or_within_year_sites",
                "outcome":outcome,"support":support}
    if weight=="island_equal":
        weights=1/x.groupby("island").island.transform("size").to_numpy(float)
    elif weight=="siteyear_equal":
        weights=np.ones(len(x))
    elif weight=="island_year_equal":
        weights=1/x.groupby("island_year").island_year.transform("size").to_numpy(float)
    else:raise ValueError("Unknown source weighting")
    formula=(outcome+" ~ prior_year_log1p_PCNL"+
             (" + elevation_km" if covariate_elevation else "")+
             " + C(island_year)")
    model=smf.wls(formula,x,weights=weights).fit(
        cov_type="cluster",cov_kwds={"groups":x.island,"use_correction":True})
    coef="prior_year_log1p_PCNL"
    b=float(model.params[coef]);se=float(model.bse[coef])
    require(np.isfinite([b,se]).all() and se>0,"Moth source fit singular")
    q=float(t.ppf(.975,G-1))
    return {"status":"estimated_exploratory_original_field_data",
            "outcome":outcome,"source_unit":"one sampled site-year; 30 true islands maximum",
            "within_island_year_FEs":True,
            "weights":weight,"formula":formula,
            "n_original_islands_clustered":G,"support":support,
            "beta_per_log1p_past_year_PCNL":b,
            "island_cluster_SE":se,
            "95CI_nominal_cluster_islands":[b-q*se,b+q*se],
            "nominal_unadjusted_two_sided_p":float(2*t.sf(abs(b/se),G-1)),
            "model_R_squared_descriptive":float(model.rsquared),
            "not_independent_species_replicates":True,
            "published_trap_nights_unavailable":True,
            "light_trap_catch_not_local_population_density_or_pollination":True}

def run(prior_source,prior_receipt,out):
    d,source=source_frame(prior_source,prior_receipt)
    primary="rarefied20_expected_species"
    outcomes=(primary,"log1p_moths_captured","log1p_species_richness",
              "rarefied10_expected_species","hill1_effective_species")
    models={}
    for name in outcomes:
        models[name]=fit(d,name,"island_equal",True)
    models["primary_site_year_equal"]=fit(d,primary,"siteyear_equal",True)
    models["primary_island_year_equal"]=fit(d,primary,"island_year_equal",True)
    models["primary_no_elevation"]=fit(d,primary,"island_equal",False)
    source_support={
        "all_true_source_moth_samples":len(d),
        "all_true_islands":int(d.island.nunique()),
        "all_original_observed_moth_individuals":int(d.original_moths_captured.sum()),
        "all_source_moth_species":1179,
        "source_samples_with_valid_prior_PCNl":int(d.prior_year_log1p_PCNL.notna().sum()),
        "source_islands_with_valid_prior_PCNl":int(
            d.loc[d.prior_year_log1p_PCNL.notna(),"island"].nunique()),
        "source_surveys_with_20_or_more_original_moths":int(
            d.original_moths_captured.ge(20).sum()),
        "source_surveys_with_10_or_more_original_moths":int(
            d.original_moths_captured.ge(10).sum()),
        "island_year_groups":int(d.island_year.nunique()),
        "site_ID_year_repeats":int(d.groupby(["island","sites"]).year.nunique().ge(2).sum()),
        "original_sampling_years":sorted(d.year.unique().astype(int).tolist()),
        "no_trap_night_count_or_date_in_source":True}
    report={
        "schema":SCHEMA,
        "source_DOI_real_moth_island_data":"10.5281/zenodo.21438943",
        "source_DOI_year_prior_PCNL":"10.5281/zenodo.17013414",
        "original_PCNL_site_year_table_sha256":checksum(prior_source),
        "source_quality":source_support,
        "model_specifications":models,
        "rarefaction_sampling_of_individuals_not_standardized_trap_nights":True,
        "no_moth_density_or_actual_pollen_delivery":True,
        "no_evidence_for_true_night_pollinator_gain_or_loss_from_light_traps_alone":True,
        "no_posthoc_significance_based_selection":True,
        "no_GloBI_plant_moth_claims_used":True
    }
    out=Path(out)
    require(not out.exists(),"Do not overwrite signed real outcome model")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_TRUE_KOREAN_ISLAND_MOTH_ABUNDANCE_PRIOR_ALAN_RESULT",
          json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--receipt",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.receipt,a.out)
