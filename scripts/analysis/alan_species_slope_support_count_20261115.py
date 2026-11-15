#!/usr/bin/env python3
"""Count species-level historical ALAN support from genuine photo-flora.
No fitted light effects, species×site coordinates, or restricted traits exported.
Photo nonrecord must never be interpreted as true plant absence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, LIGHT, load_sites, load_labels, original_accepted_species, require
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits
from scripts.analysis.alan_real_daynight_reproductive_PCNL_562_20261113 import experimental_daynight

def island_fold(island):
    return hashlib.sha256(str(island).encode("utf-8")).digest()[0] & 1

def count_at_threshold(sites, flora, taxa, threshold):
    """Island×year light contrast, then species×survey-opportunity support."""
    assert not sites.duplicated(KEY).any()
    assert not flora.duplicated(KEY+["species"]).any()
    g=(sites.groupby(["island","year"],observed=True)[LIGHT]
       .agg(n_sites="size",light_min="min",light_max="max",
            median_light="median").reset_index())
    g=g.loc[g.n_sites.ge(6)&(g.light_max-g.light_min).ge(threshold)].copy()
    cand=sites.merge(g[["island","year","median_light"]],
                     on=["island","year"],how="inner",validate="many_to_one")
    cand["hi"]=(cand[LIGHT]>cand.median_light).astype("int8")
    halves=(cand.groupby(["island","year","hi"],observed=True).size()
            .unstack("hi",fill_value=0))
    for k in (0,1):
        if k not in halves:halves[k]=0
    ok=halves.loc[halves[0].ge(2)&halves[1].ge(2)].reset_index()[
        ["island","year"]]
    cand=cand.merge(ok,on=["island","year"],validate="many_to_one")
    cand["fold"]=cand.island.map(island_fold).astype("int8")
    obs=flora.merge(cand[KEY+["hi","fold"]],on=KEY,
                    how="inner",validate="many_to_one")
    obs["lo"]=1-obs.hi
    totals=(cand.groupby(["island","year"],observed=True).size()
            .rename("surveyed_sites").reset_index())
    obs_group=(obs.groupby(["species","island","year","fold"],observed=True)
               .agg(n=("species","size"),lo=("lo","sum"),hi=("hi","sum"))
               .reset_index())
    obs_group=obs_group.merge(totals,on=["island","year"],
                              validate="many_to_one")
    # Presence and nonrecord must vary within an otherwise light-varying
    # observational stratum. Never claim that nonrecord equals true absence.
    informative=obs_group.loc[obs_group.n.lt(obs_group.surveyed_sites)].copy()
    agg=(informative.groupby("species",observed=True)
         .agg(n=("n","sum"),strata=("n","size"),islands=("island","nunique"),
              years=("year","nunique"),lo=("lo","sum"),hi=("hi","sum"))
         .reset_index())
    agg["screen"]=agg.n.ge(12)&agg.strata.ge(2)
    agg["core"]=(agg.n.ge(30)&agg.strata.ge(6)&
                 agg.islands.ge(2)&agg.years.ge(2))
    agg["two_halves"]=agg.core&agg.lo.ge(5)&agg.hi.ge(5)
    fold=(informative.groupby(["species","fold"],observed=True)
          .agg(n=("n","sum"),strata=("n","size")).reset_index())
    fold["supported"]=fold.n.ge(10)&fold.strata.ge(2)
    held=fold.groupby("species",observed=True).supported.sum().rename(
        "supported_folds").reset_index()
    agg=agg.merge(held,on="species",how="left",validate="one_to_one")
    agg["holdout"]=agg.supported_folds.ge(2)
    agg["core_holdout"]=agg.core&agg.holdout
    agg["strict"]=agg.two_halves&agg.holdout
    agg=taxa.merge(agg,on="species",how="left",validate="one_to_one")
    indicator=["screen","core","two_halves","holdout","core_holdout","strict"]
    for c in indicator:agg[c]=agg[c].eq(True)
    classes={
      "all_observed":agg.species.notna(),
      "selfing_known_high_medium":agg.selfing_hm.isin(["autonomous","absent"]),
      "selfing_autonomous":agg.selfing_hm.eq("autonomous"),
      "selfing_absent":agg.selfing_hm.eq("absent"),
      "selfing_known_high_only":agg.selfing_high.isin(["autonomous","absent"]),
      "experimental_diel":agg.diel.isin(["night","day"]),
      "experimental_night_stronger":agg.diel.eq("night"),
      "experimental_day_stronger":agg.diel.eq("day"),
      "selfing_and_experimental_diel":(
           agg.selfing_hm.isin(["autonomous","absent"])&
           agg.diel.isin(["night","day"])),
      "selfing_absent_and_night_stronger":(
           agg.selfing_hm.eq("absent")&agg.diel.eq("night")),
      "selfing_autonomous_and_night_stronger":(
           agg.selfing_hm.eq("autonomous")&agg.diel.eq("night"))
    }
    counts={label:{"matched_species":int(mask.sum()),
                   **{c:int((mask&agg[c]).sum()) for c in indicator}}
            for label,mask in classes.items()}
    return {"exposure_span_min_log1p":threshold,
            "eligible_photo_site_years":int(len(cand)),
            "eligible_islands":int(cand.island.nunique()),
            "eligible_island_year_strata":int(len(ok)),
            "eligible_species_site_years":int(len(obs)),
            "distinct_species_in_eligible_strata":int(obs.species.nunique()),
            "species_with_informative_island_year":int(informative.species.nunique()),
            "by_trait_evidence":counts}

def run(sites_path,globi_path,groups,ledger,study,study_traits,out):
    assert not out.exists()
    sites=load_sites(sites_path)
    flora=original_accepted_species(groups,sites,load_labels(globi_path))
    require((len(sites),sites.island.nunique(),len(flora),
             flora.species.nunique())==(77647,562,2012925,14330),
            "Original source frame has changed")
    hm=direct_plant_traits(ledger,high_only=False)
    high=direct_plant_traits(ledger,high_only=True)
    diel,_=experimental_daynight(study,study_traits)
    require((len(diel),int(diel.night_stronger.sum()),
             int(diel.day_stronger.sum()))==(134,75,59),
            "Original independent day/night source labels drifted")
    require(diel.species.is_unique,"Day/night duplicate plant")
    taxa=pd.DataFrame({"species":flora.species.drop_duplicates()})
    taxa["selfing_hm"]=taxa.species.map(hm.autonomous_selfing_capacity)
    taxa["selfing_high"]=taxa.species.map(high.autonomous_selfing_capacity)
    diel_map=pd.Series(np.where(diel.night_stronger,"night","day"),
                       index=diel.species)
    taxa["diel"]=taxa.species.map(diel_map)
    source={
      "islands":562,"site_years":len(sites),
      "source_species_site_years":len(flora),"accepted_species":len(taxa),
      "direct_selfing_global_high_medium":int(
         hm.autonomous_selfing_capacity.isin(["autonomous","absent"]).sum()),
      "direct_selfing_photographed_high_medium":int(
         taxa.selfing_hm.isin(["autonomous","absent"]).sum()),
      "direct_selfing_photographed_high_only":int(
         taxa.selfing_high.isin(["autonomous","absent"]).sum()),
      "direct_selfing_photographed_excluding_diel_experimental_species":int(
         (taxa.selfing_hm.isin(["autonomous","absent"])&taxa.diel.isna()).sum()),
      "daynight_experiment_species_photographed":int(taxa.diel.notna().sum()),
      "matched_selfing_and_diel":int(
         (taxa.selfing_hm.isin(["autonomous","absent"])&taxa.diel.notna()).sum())
    }
    gates={str(t):count_at_threshold(sites,flora,taxa,t) for t in (.02,.10)}
    result={
      "schema":"ALAN_species_slope_source_support_counts_v1",
      "status":"SOURCE_CENSUS_ONLY_SPECIES_SLOPES_NOT_YET_FITTED",
      "source":source,"light_gradient_gates":gates,
      "eligibility_rules":{
        "stratum":"same island × original observation year",
        "photo_cells_per_stratum":6,
        "photographed_cells_per_light_half":2,
        "screen":"12 positive species-cell years and 2 within-stratum informative groups",
        "core":"30 positive species-cell years; 6 informative island×years; 2 islands; 2 observed years",
        "two_halves":"core and >=5 positive photo cells in both light strata",
        "holdout":"whole islands hash-split into two folds; each >=10 positive cells and 2 informative island×years",
        "strict":"core and two-halves and holdout"
      },
      "limitations":{
        "present_absent_is_photographed_present_vs_not_recorded":True,
        "no_original_population_trend_or_seed_output":True,
        "island_photo_archive_only_not_global_GBIF_coverage":True,
        "no_worldpop_climate_or_observer_adjustment_in_eligibility_count":True,
        "island_hash_folds_are_not_spatially_separated_oceanic_regions":True,
        "no_species_slope_or_meta_regression_fitted":True,
        "night_benefit_is_not_proven_moth_pollination":True,
        "original_trait_by_site_not_redistributed":True
      }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_SLOPE_SUPPORT_FINAL_COUNTS",
          json.dumps({"source":source,
              "gates":{k:{"islands":v["eligible_islands"],
                          "site_years":v["eligible_photo_site_years"],
                          "trait_counts":v["by_trait_evidence"]}
                  for k,v in gates.items()}},ensure_ascii=False),
          flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("sites","globi","groups","ledger","study","study_traits","out"):
        p.add_argument("--"+a,required=True,type=Path)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.ledger,a.study,a.study_traits,a.out)
