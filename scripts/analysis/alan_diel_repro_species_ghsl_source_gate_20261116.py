#!/usr/bin/env python3
"""Directly observed day/night reproductive *function* versus ALAN photo slopes.

External Kendall/Nicholson original 2025 day-vs-night exclusion reproduction
is NOT ALAN manipulation, not a locally matched population, not moth-only
pollination and not intrinsic selfing. Reconstruct author study species from
pinned original 1652 treatment rows and 142 independent floral trait rows.

A globally observed original 562-island species photo source and native
WorldPop + WorldClim + independent JRC GHSL population/urban controls are
reused unchanged. Fit one per-species photographic source slope with same
island-year strata and urban built-up controls as PR #84. Quantify source
support BEFORE choosing whether an island bootstrap is identifiable.

No source species×location, original site-level reproduction or original
endangered coordinates are exported.
"""
from __future__ import annotations
import argparse, json, numpy as np, pandas as pd
from pathlib import Path
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    load_sites,load_labels,original_accepted_species,require
)
from scripts.analysis.alan_real_daynight_reproductive_PCNL_562_20261113 import experimental_daynight
from scripts.analysis.alan_species_worldpop_conditional_20261115 import (
    within_year_island_photo_design,candidate_species_photo_support
)
from scripts.analysis.alan_species_worldclim_population_20261115 import matching_climate_gradient
from scripts.analysis.alan_species_full_worldpop_climate_20261116 import source_annual_population
from scripts.analysis.alan_species_ghsl_builtup_urban_control_20261116 import (
    authentic_ghsl_source,fit_source_species
)
from scripts.analysis.alan_species_effortadjusted_firststage_20261115 import percentile

SEED=20261116
BOOT=400
MIN_BOOT_SPECIES_PER_GUILD=8
MIN_BOOT_ISLANDS=10

def source_diel_labels(study,traits,photographed_species):
    d,original=experimental_daynight(study,traits)
    require((len(d),int(d.night_stronger.sum()),int(d.day_stronger.sum()))==
            (134,75,59),"Original author diel reproduction source drift")
    labels=d.loc[d.species.isin(photographed_species)].copy()
    require(len(labels)==52 and int(labels.night_stronger.sum())==27 and
            int(labels.day_stronger.sum())==25,
            "Original fully accepted WCVP photo species×diel experiment match drift")
    require(labels.species.is_unique,"Duplicated experimental plant taxon")
    labels["guild"]=np.where(labels.night_stronger,"night_higher","day_higher")
    labels["genus"]=labels.species.str.split().str[0]
    labels=labels.set_index("species")
    # Independence checks: external empirical comparison across original source
    # plant families and original publications, not local pollen delivery.
    mixed_families=(labels.groupby("plant_family",observed=True).guild.nunique())
    mixed_genera=(labels.groupby("genus",observed=True).guild.nunique())
    support={
      "original_independent_daynight_original_reproductive_species":len(d),
      "exact_botanical_photo_species_matched":len(labels),
      "night_higher_species":int(labels.guild.eq("night_higher").sum()),
      "day_higher_species":int(labels.guild.eq("day_higher").sum()),
      "observed_original_study_families":int(labels.plant_family.nunique()),
      "observed_original_genera":int(labels.genus.nunique()),
      "families_with_both_night_higher_and_day_higher":int(mixed_families.eq(2).sum()),
      "genera_with_both_night_higher_and_day_higher":int(mixed_genera.eq(2).sum()),
      "source_species_with_available_perhour_exclusion":int(labels.perhour_available.sum()),
      "perhour_source_class_discordant_from_unadjusted":int(
          ((labels.guild.eq("night_higher")&labels.perhour_day_stronger)|
           (labels.guild.eq("day_higher")&labels.perhour_night_stronger)).sum()),
      "reproductive_published_source_not_ALAN_light_treatment":True,
      "no_local_species_population_pollinator_identity":True
    }
    return labels,original,support

def summary(fit):
    d=pd.DataFrame(fit)
    cols=("L_baseline_pp","L_built2015_pp","L_built2020_pp",
          "built2015_pp_per_log_percent_built",
          "built2020_pp_per_log_percent_built")
    groups={"all":d,
        "published_night_repro_higher":d.loc[d.guild.eq("night_higher")],
        "published_day_repro_higher":d.loc[d.guild.eq("day_higher")]}
    out={}
    for key,g in groups.items():
        out[key]={"n_plants":len(g),
           "source_island_ge10":int(g.n_islands.ge(10).sum()),
           "source_island_ge20":int(g.n_islands.ge(20).sum()),
           "source_island_ge30":int(g.n_islands.ge(30).sum()),
           "published_source_families":int(g.plant_family.nunique()),
           **{col:percentile(g[col]) for col in cols}}
    for col in cols:
        g1=groups["published_night_repro_higher"][col]
        g2=groups["published_day_repro_higher"][col]
        out[col+"_night_minus_day"]=float(g1.mean()-g2.mean()) if len(g1) and len(g2) else None
    return out

def bootstrap_diel(genuine_models,universe):
    chosen=[v for v in genuine_models if v["q"]["n_islands"]>=MIN_BOOT_ISLANDS]
    night=sum(v["guild"]=="night_higher" for v in chosen)
    day=len(chosen)-night
    result={"strict_source_support":"AT_LEAST_10_ISLANDS_PER_SPECIES",
            "n_source_eligible":len(chosen),
            "n_night_higher":night,"n_day_higher":day,
            "bootstrap_species_per_guild_minimum":MIN_BOOT_SPECIES_PER_GUILD}
    if min(night,day)<MIN_BOOT_SPECIES_PER_GUILD:
        result["status"]="NOT_ESTIMABLE_DIRECT_DIEL_SOURCE_INDEPENDENT_SPECIES_SUPPORT"
        return result
    U=list(universe);idx={x:i for i,x in enumerate(U)}
    n=len(chosen);g=len(U)
    A=np.zeros((n,g,6,6),float)
    B=np.zeros((n,g,6),float)
    is_night=np.zeros(n,bool)
    for i,row in enumerate(chosen):
        is_night[i]=(row["guild"]=="night_higher")
        j=np.array([idx[z] for z in row["islands"]])
        A[i,j,:,:]=row["A"][:,:6,:6]
        B[i,j,:]=row["B"][:,:6]
    def solve(mult):
        aa=np.einsum("g,ngij->nij",mult,A,optimize=True)
        bb=np.einsum("g,ngi->ni",mult,B,optimize=True)
        return np.linalg.solve(aa,bb[...,None])[...,0]
    source=solve(np.ones(g))
    target=np.asarray([v["q"]["L_built2015_pp"] for v in chosen],float)
    require(np.max(abs(100*source[:,0]-target))<1e-7,
            "Actual source original species coefficients cannot reconstruct from island score blocks")
    rng=np.random.default_rng(SEED)
    slopes=np.full(BOOT,np.nan)
    for i in range(BOOT):
        w=np.bincount(rng.integers(0,g,g),minlength=g)
        try:
            draw=solve(w.astype(float))*100
        except np.linalg.LinAlgError:continue
        if not np.isfinite(draw).all():continue
        slopes[i]=draw[is_night,0].mean()-draw[~is_night,0].mean()
    vals=slopes[np.isfinite(slopes)]
    if len(vals)<.8*BOOT:
        result.update(status="NOT_ESTIMABLE_SOURCE_ISLAND_BOOTSTRAP_RANK",
                     n_valid_boot=len(vals),n_requested=BOOT)
        return result
    result.update(status="EXPLORATORY_SOURCE_DIEL_REPRO_FUNCTION_LIGHT_SLOPE_COMPARATOR",
                  n_islands_resampled=g,
                  observed_night_minus_day_light_slope_pp=float(
                      target[is_night].mean()-target[~is_night].mean()),
                  n_bootstrap_valid=len(vals),n_bootstrap_requested=BOOT,
                  ci95_island_percentile_pp=[float(x) for x in np.quantile(vals,[.025,.975])],
                  original_published_study_group_not_local_pollinator_service=True,
                  no_causal_ALAN_or_moth_function_claim=True)
    return result

def run(site_file,globi,groups,annual_worldpop,climate_root,ghsl_source,
        study,study_traits,out):
    require(not out.exists(),"Refuse source overwrite")
    original=load_sites(site_file)
    flora=original_accepted_species(groups,original,load_labels(globi))
    labels,origin,meta=source_diel_labels(
        study,study_traits,set(flora.species.unique()))
    pops,pmeta=source_annual_population(original,annual_worldpop)
    require(pmeta["source_population_valid_sites"]==69525,
            "Full authenticated WorldPop source missing")
    photo,photo_meta=within_year_island_photo_design(pops)
    require(len(photo)==64966 and photo.island.nunique()==202,
            "Prior frozen 703-species photo support changed")
    climate,climate_meta=matching_climate_gradient(photo,climate_root)
    require(len(climate)==64944,"Frozen prior WorldClim source drift")
    built,b_meta=authentic_ghsl_source(ghsl_source,climate)
    pos,iy,rawsupport=candidate_species_photo_support(
        built,flora,labels.guild)
    data=[];blocks=[];failure={}
    for name in sorted(labels.index):
        if name not in pos or name not in iy:
            why="NO_INFORMATIVE_ISLAND_YEAR_FOR_PHOTOGRAPHED_PLANT"
            failure[why]=failure.get(why,0)+1
            continue
        ans,reason=fit_source_species(built,pos[name],iy[name])
        if reason:
            failure[reason]=failure.get(reason,0)+1
            continue
        q,islands,A,B=ans
        guild=str(labels.loc[name,"guild"])
        family=str(labels.loc[name,"plant_family"])
        data.append({**q,"guild":guild,"plant_family":family})
        blocks.append({"q":q,"guild":guild,"islands":islands,"A":A,"B":B})
    group_results=summary(data) if data else {"status":"NOT_ESTIMABLE_NO_DIRECT_DIEL_STUDY_SPECIES"}
    model_support={
       **rawsupport,
       "source_direct_diel_documented_plants":len(labels),
       "estimated_plants":len(data),
       "night_higher_estimated":sum(x["guild"]=="night_higher" for x in data),
       "day_higher_estimated":sum(x["guild"]=="day_higher" for x in data),
       "nonestimable_source_reasons":failure}
    report={
      "schema":"ALAN_direct_published_night_day_reproduction_guild_vs_photoALAN_urban_controls_20261116_v1",
      "status":"SOURCE_AUTHENTIC_PUBLISHED_DIEL_SPECIES_FEASIBILITY_AND_EXPLORATORY_SLOPE_DIAGNOSTIC",
      "source":{
         "original_plant_species_site_years":2012925,
         "original_ALAN_photo_sites":77647,
         "actual_worldpop_photo_sites":69525,
         "source_predefined_autonomous_selfing_703_cohort_not_used_as_selection":True,
         "source_original_diel_author_data_pinned_by_Git_blob":True,
         "prior_2015_GHSL_raster_sha256":b_meta["source_ghsl_2015_raster_sha256"],
         "photo_source_after_population_climate_built":len(built),
         "photo_source_islands":int(built.island.nunique())
      },
      "actual_published_daynight_exclosure_source":meta,
      "original_source_eligible_support":model_support,
      "descriptive_species_light_slope_moderation":group_results,
      "original_whole_island_bootstrap":bootstrap_diel(blocks,built.island.unique()),
      "limitations":{
        "original_external_day_night_fruitseed_function_not_locally_matched_to_ALAN":True,
        "original_function_not_proof_of_moth_effective_pollination":True,
        "no_light_experimental_reproduction_comparison":True,
        "photo_nonrecord_is_not_true_species_absence":True,
        "source_species_guild_heterogeneous_plant_lineages_and_papers":True,
        "genus_and_study_dependence_not_fully_adjusted":True,
        "2015_builtup_only_proxy_urban_habitat_not_full_landcover":True,
        "previous_year_PCNL_not_local_nighttime_lux":True,
        "no_raw_protected_species_site_or_original_reproduction_locality_exported":True,
        "not_a_confirmatory_meta_analysis":True
      }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ORIGINAL_DIEL_FUNCTION_GHSL_SOURCE",json.dumps({
        "diel":meta,"support":model_support,"groups":group_results,
        "bootstrap":report["original_whole_island_bootstrap"]
    },ensure_ascii=False,allow_nan=False),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("sites","globi","groups","annual_worldpop",
                 "climate_root","ghsl_source","study","study_traits","out"):
        p.add_argument("--"+flag,type=Path,required=True)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.annual_worldpop,a.climate_root,
        a.ghsl_source,a.study,a.study_traits,a.out)
