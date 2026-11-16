#!/usr/bin/env python3
"""Source-complete island ALAN species response with original annual WorldPop.

Join all original ALAN 77,647 photo-site-year rows by their checksum-locked
row position to ORIGINAL five-year, source-SHA-verified WorldPop rasters. There
is no dependency on another photographic frame for population measurement.

Validate the previous 23,018 WorldPop exact-pixel overlap estimates AGAINST
independently re-extracted full-coverage population, not replacing genuine
NoData or relabeling survey years.

Recompute fixed preselected 703 measured-autonomous-selfing species photo
record indicator slopes using same-stratum island×observed-year FE, original
photo effort, previous-calendar-year PCNL, source-year WorldPop and original
WorldClim normals. Repeat whole-island block bootstrap with species-stage
refits. These remain correlational photograph-record analyses.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, load_sites, load_labels, original_accepted_species, require
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import direct_plant_traits
from scripts.analysis.alan_species_effortadjusted_firststage_20261115 import prepare_support, percentile
from scripts.analysis.alan_species_worldpop_conditional_20261115 import (
    within_year_island_photo_design, candidate_species_photo_support, fit_species
)
from scripts.analysis.alan_species_worldclim_population_20261115 import (
    matching_climate_gradient, fit_five, output_summary, boot_islands
)
from scripts.analysis.alan_mois_worldpop_exact_pixel_support_20261115 import join as prior_mois_pixel_join
from scripts.analysis.alan_worldpop_full_77647_sites_20261116 import (
    RASTERS, file_sha256
)

ORIGINAL_N=77647
BASE_SOURCE_SHA="440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c"
N_BOOT=400

def source_annual_population(all_sites,annual_dir):
    require(len(all_sites)==ORIGINAL_N and all_sites.island.nunique()==562,
            "Source population and source flora frame drift")
    all_sites=all_sites.reset_index(drop=True).copy()
    all_sites["site_id"]=np.arange(len(all_sites),dtype=np.int32)
    original_years=all_sites.groupby("year",observed=True).size().to_dict()
    data=[];annual=[]
    for year in sorted(RASTERS):
        path=Path(annual_dir)/f"year{year}"
        receipt=path/"source_receipt.json"
        values=path/f"worldpop_R2024B_UC_{year}_original_ALAN_photo_sites.csv.gz"
        require(receipt.is_file() and values.is_file(),
                "Missing original full-sampled annual WorldPop source")
        rec=json.loads(receipt.read_text())
        reference_sha, reference_bytes, reference_n=RASTERS[year]
        require(rec["schema"]=="ALAN_original_full_photo_island_WorldPop_G2_R2024B_UC_v1"
            and rec["same_year_only"] and rec["population_year"]==year
            and rec["original_raster_bytes"]==reference_bytes
            and rec["raster_metadata"]["source_original_sha256"]==reference_sha
            and rec["photo_site_source_SHA256"]==BASE_SOURCE_SHA
            and rec["original_alansites"]==reference_n==original_years[year],
            "Wrong original source/raster/year receipt")
        require(file_sha256(values)==rec["new_population_source_values_sha256"],
                "Original annual source population values SHA changed")
        v=pd.read_csv(values,dtype={"site_id":np.int32,"population_status":str})
        require(len(v)==reference_n and v.site_id.is_unique and
                v.original_observation_year.eq(year).all() and
                v.worldpop_year.eq(year).all(),
                "Annual original population one-to-one key/year mismatch")
        require(v.loc[v.population_status.ne("OK"),
                      "population_density_per_km2"].isna().all() and
                v.loc[v.population_status.eq("OK"),
                      "population_density_per_km2"].ge(0).all(),
                "Original population NODATA or true-zero conflated")
        data.append(v[["site_id","population_status","population_density_per_km2"]])
        annual.append({
            "year":year,"all_photo_source_site_years":reference_n,
            "source_population_valid":int(v.population_status.eq("OK").sum()),
            "source_population_nodata":int(v.population_status.eq("NODATA").sum()),
            "source_outside_raster":int(v.population_status.eq("OUTSIDE_RASTER").sum()),
            "source_true_zero_population":int(
                 (v.population_status.eq("OK")&
                  v.population_density_per_km2.eq(0)).sum()),
            "source_raster_sha256":reference_sha,
            "source_site_population_csv_sha256":rec["new_population_source_values_sha256"]
        })
    merged_values=pd.concat(data,ignore_index=True)
    require(len(merged_values)==len(all_sites) and merged_values.site_id.is_unique,
            "At least one original ALAN photographic site/year missing")
    s=all_sites.merge(merged_values,on="site_id",how="left",validate="one_to_one")
    require(s.population_status.notna().all(),"Raw year extraction failed to cover one or more source sites")
    s["pop_density"]=s.population_density_per_km2
    status=s.population_status.value_counts().to_dict()
    status_by_selection=s.groupby("population_status").agg(
        source_site_years=("site_id","size"),
        distinct_islands=("island","nunique"),
        median_preexisting_light=("pcnl_1992_to_previous_year_mean_log1p","median"),
        median_plant_observation_effort=("all_accepted_observation_uuids","median"),
    )
    strata={
      str(k):{"site_years":int(v.source_site_years),
              "islands":int(v.distinct_islands),
              "median_historical_logPCNL":float(v.median_preexisting_light),
              "median_original_plant_photo_effort":float(v.median_plant_observation_effort)}
      for k,v in status_by_selection.iterrows()
    }
    return s,{"annual_original_worldpop_source":annual,
              "WorldPop_original_NODATA_selection_bias_source_diagnostic":strata,
              "original_alan_source_photo_site_years":len(s),
              "original_alan_islands":int(s.island.nunique()),
              "source_population_status_counts":{k:int(v) for k,v in status.items()},
              "source_population_valid_sites":int(s.pop_density.notna().sum()),
              "source_population_valid_islands":int(
                  s.loc[s.pop_density.notna(),"island"].nunique()),
              "source_population_true_zeros":int(s.pop_density.eq(0).sum()),
              "no_source_photo_site_dropped_before_true_NODATA_gate":True,
              "all_years_exact_original_WorldPop_pixels_and_SHA":True,
              "population_not_interpolated_from_MOIS_photo_sites":True}

def audit_original_overlap(s,previous_mois_sources):
    ref, subset=prior_mois_pixel_join(
        s.drop(columns=["pop_density"]),previous_mois_sources,return_data=True)
    require(ref["matched_usable_site_years"]==23018,
            "Previous WorldPop-MOIS photographic overlap changed")
    known=subset.loc[subset.pop_density.notna(),
                     ["site_id","pop_density"]].rename(
                     columns={"pop_density":"pop_density_previous"})
    assert known.site_id.is_unique
    check=s.merge(known,on="site_id",how="inner",validate="one_to_one")
    require(len(check)==23018,"Missing original pixel-matched photo site for validation")
    require(check.pop_density.notna().all(),"Full source lost population at previously measured original pixels")
    old=check.pop_density_previous.to_numpy(float)
    new=check.pop_density.to_numpy(float)
    abs_err=np.abs(old-new)
    rel_err=abs_err/np.maximum(1,abs(old))
    require(np.all(np.isclose(old,new,rtol=1e-5,atol=1e-4)),
            "Official full-raster population differs from prior same-year same-pixel original")
    s=s.copy()
    s["previous_MOIS_photo_overlap"]=s.site_id.isin(known.site_id)
    cover=s.groupby("previous_MOIS_photo_overlap").agg(
         n_sites=("site_id","size"),
         n_valid_population=("pop_density","count"),
         n_islands=("island","nunique"),
         median_population=("pop_density","median"),
         median_historical_logPCNL=("pcnl_1992_to_previous_year_mean_log1p","median"),
         median_plant_photo_observations=("all_accepted_observation_uuids","median"))
    separated={}
    for flag,row in cover.iterrows():
        separated["prior_MOIS_photo_overlap" if flag else "new_direct_ALAN_photo_sites"]={
             k:(int(row[k]) if k in ("n_sites","n_valid_population","n_islands")
                else float(row[k])) for k in cover.columns}
    return s,{"previous_exact_pixel_source_matched_sites":len(check),
              "maximum_absolute_population_difference_per_km2":float(abs_err.max()),
              "median_absolute_population_difference_per_km2":float(np.median(abs_err)),
              "maximum_relative_population_difference_at_least_one_person_per_km2":
                  float(rel_err.max()),
              "prior_and_new_population_original_raster_values_reconciled":True,
              "new_original_source_sites_not_selected_by_MOIS_photo_opportunity":True,
              "old_photo_overlap_vs_new_direct_sites":separated}

def run(site_file,globi,groups,ledger,annual_worldpop,
        prior_mois_sources,climate_root,out):
    require(not out.exists(),"Original full-coverage result overwrite refused")
    original=load_sites(site_file)
    flora=original_accepted_species(groups,original,load_labels(globi))
    trait=direct_plant_traits(ledger)
    plants=pd.DataFrame({"species":flora.species.drop_duplicates()})
    plants["autonomy"]=plants.species.map(trait.autonomous_selfing_capacity)
    _,_,_,fixed703,s0=prepare_support(original,flora,plants,span=.02)
    require((s0["strict_selfing_species"],len(fixed703))==(703,703),
            "Frozen original 703 direct-selfing trait cohort drift")
    original_labels=plants.set_index("species").autonomy.loc[sorted(fixed703)]
    population,source_receipt=source_annual_population(original,annual_worldpop)
    population,overlap=audit_original_overlap(population,prior_mois_sources)
    base,baseinfo=within_year_island_photo_design(population)
    climate,climinfo=matching_climate_gradient(base,climate_root)
    ids,strata,coverage=candidate_species_photo_support(
        climate,flora,original_labels)
    fits=[];blocks=[];fail={}
    for j,sp in enumerate(sorted(fixed703)):
        if sp not in ids or sp not in strata:
            fail["NO_SOURCE_INFORMATIVE_ISLAND_YEAR"]=fail.get(
                "NO_SOURCE_INFORMATIVE_ISLAND_YEAR",0)+1
            continue
        res,why=fit_five(climate,ids[sp],strata[sp])
        if why is not None:
            fail[why]=fail.get(why,0)+1
            continue
        q,islands,aa,bb=res
        old,old_reason=fit_species(climate,ids[sp],strata[sp])
        require(old_reason is None,"Previously supported light+effort same-source baseline became unavailable")
        oldfit,*_=old
        q["beta_L_effort_only_pp"]=oldfit["beta_L_without_population_pp"]
        aut=str(original_labels[sp])
        fits.append({**q,"autonomy":aut})
        blocks.append({"q":q,"autonomy":aut,"islands":islands,
                       "AA":aa,"BB":bb})
        if (j+1)%100==0:
            print("ALAN_ORIGINAL_FULL_WORLDPOP_SPECIES_FIT_PROGRESS",
                  j+1,"of",len(fixed703),flush=True)
    require(len(fits)>0,"No original population-climate species source estimable")
    # Predeclare ≥20 independent islands as primary common-support
    # robustness gate: avoids near-singular 5-predictor fits in very sparse
    # island resampling. All 703 source results retained descriptively.
    primary_islands=20
    blocks_for_boot=[r for r in blocks if r["q"]["n_islands"]>=primary_islands]
    bootstrap=boot_islands(blocks_for_boot,climate.island.unique())
    bootstrap["predeclared_min_independent_islands_for_bootstrap"]=primary_islands
    combined=output_summary(fits)
    groups=pd.DataFrame(fits)
    combined["same_sample_light_effort_only_autonomous_minus_absent_pp"]=float(
         groups.loc[groups.autonomy.eq("autonomous"),"beta_L_effort_only_pp"].mean()-
         groups.loc[groups.autonomy.eq("absent"),"beta_L_effort_only_pp"].mean())
    combined["same_sample_light_effort_only_species_slope_pp"]=percentile(
         groups.beta_L_effort_only_pp)
    report={
       "schema":"ALAN_703_autonomy_species_FULL_annual_WorldPop_WorldClim_original_sites_v1",
       "status":"ORIGINAL_SOURCE_FULL_SITE_WORLDPOP_PLUS_STATIC_WORLDCLIM_EXPLORATORY",
       "source":{
           "original_alan_source_sites":len(original),
           "original_alan_islands":int(original.island.nunique()),
           "original_plant_species_site_years":len(flora),
           "full_703_direct_trait_species":703,
           "worldpop_model":"Global2 R2024B UC v1 same observation calendar year original 1km",
           "prior_NTL":"1992 through observation-year-minus-one original Li-Cao PCNL mean log1p",
           "worldclim_model":"BIO1 + log1p(BIO12) WorldClim 2.1 10-minute 1970–2000 normals"
       },
       "worldpop_full_direct_source":source_receipt,
       "independent_prior_23k_pixel_original_replay":overlap,
       "complete_worldpop_light_cohort":baseinfo,
       "climate_matched_source":climinfo,
       "species_support":{
          **coverage,"original_fixed703":703,
          "species_fitted":len(fits),
          "autonomous_fitted":int(sum(x["autonomy"]=="autonomous" for x in fits)),
          "autonomy_absent_fitted":int(sum(x["autonomy"]=="absent" for x in fits)),
          "nonfit_source_reasons":fail,
          "original_species_with_20plus_island_support":len(blocks_for_boot),
          "20plus_autonomous":int(sum(x["autonomy"]=="autonomous" for x in blocks_for_boot)),
          "20plus_absent":int(sum(x["autonomy"]=="absent" for x in blocks_for_boot))
       },
       "same_cohort_model_comparison":combined,
       "whole_island_bootstrap":bootstrap,
       "limitations":{
          "source_population_direct_not_MOIS_photo_overlap":True,
          "true_nodata_never_imputed_to_zero":True,
          "photo_species_nonrecord_not_true_species_absence":True,
          "no_species_population_trends_or_reproduction":True,
          "no_independently_measured_floral_illumination":True,
          "WorldClim_static_not_observation_year_weather":True,
          "no_urban_land_cover_microhabitat_adjustment":True,
          "no_observer_preference_or_independent_detection_model":True,
          "no_phylogenetic_meta_regression_estimated":True,
          "same_703_species_original_trait_ascertainment_selection":True,
          "exploratory_preselected_support_filters_no_causal_claim":True,
          "no_restricted_species_trait_by_location_exported":True
       }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,
            indent=2,allow_nan=False)+"\n")
    print("ALAN_FULL_WORLDPOP_77647_SPECIES_RESULTS",json.dumps({
        "population":source_receipt,
        "original_23k_replay":overlap,
        "matched_photo_design":baseinfo,
        "climate_design":{
            "population_worldclim_photo_site_years_after_new_light_gate":
                climinfo["population_worldclim_photo_site_years_after_new_light_gate"],
            "population_worldclim_supported_islands":
                climinfo["population_worldclim_supported_islands"],
            "population_worldclim_island_years":
                climinfo["population_worldclim_island_years"],
            "within_island_year_global_residual_predictor_corr":
                climinfo["within_island_year_global_residual_predictor_corr"]},
        "species_support":report["species_support"],
        "conditional_model":report["same_cohort_model_comparison"],
        "bootstrap":bootstrap,
    },ensure_ascii=False,allow_nan=False),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("sites","globi","groups","ledger","annual_worldpop",
                 "prior_mois_sources","climate_root","out"):
        p.add_argument("--"+flag,required=True,type=Path)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.ledger,a.annual_worldpop,
        a.prior_mois_sources,a.climate_root,a.out)
