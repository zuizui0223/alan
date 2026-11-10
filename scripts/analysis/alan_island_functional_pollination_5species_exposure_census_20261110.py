#!/usr/bin/env python3
"""Audited macro opportunity for INDEPENDENTLY VERIFIED island non-bee pollination.

Original source: 93,209,351 iNaturalist photos already collapsed to
WCVP accepted plants / GSHHG island photo site-years; historic PCNL strictly
before original photo observation years, original Kew WCVP botanical-country
native source, independently measured 2026-08 GBIF documented gecko source.
This is an EXPOSURE and source-support census, NOT plant decline or
demonstrated evolutionary transitions to night/vertebrate pollination.

NO row for a location of a critically endangered plant is released.
Only aggregated species×ARCHIPELAGO summaries and island counts are saved.
A literature-confirmed archipelago endemic != verified wild individual.
Source Gekkonidae photo occurrence != local predator activity or occupancy.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_island_endemic_night_moth_pair_support_20261110 import (
    sources,SEEDS,SITE_SHA)
from scripts.analysis.alan_wcvp_regional_origin_pcnl_20261110 import SCHEMA as WCVP_SCHEMA

KEY=["island","latbin","lonbin","year"]
SCHEMA="ALAN_20261110_independent_island_endemic_functional_pollination_PCNl_gecko_exposure_gate_v1"
NATIVE_SOURCE_SHA="2f30984f627425c7dff96d7502c30bb1d28d8291835961fc059e8956204bc457"
GECKO_SOURCE_SHA="c2e0c823ee9fcf61182b920724e488d9130b1d1ca9d53266511c2da8873983c9"
HEMIDACTYLUS_SOURCE_SHA="91d09507f84e062632feed8aef18f09c1fbf7797c812ed3ab147b5e32851b9e2"
LIGHT="pcnl_1992_to_previous_year_mean_log1p"

# Scientific evidence tier is verified independently of GBIF and GloBI:
# DO NOT call any of these an independently proven ancestral evolutionary switch.
# Note: "lizard pollination" in Canarias can be Gallotia, not gecko.
CONFIRMED={
    "Ipomoea habeliana":{"guild":"effective_nocturnal_hawkmoth",
       "archipelago":"Galapagos","evidence":"doi:10.1111/j.1095-8339.2009.00963.x",
       "evolutionary_ancestral_shift_proven":False},
    "Canarina canariensis":{"guild":"effective_opportunistic_birds",
       "archipelago":"Canary","evidence":"PMID:18931830; PMID:21900611",
       "evolutionary_ancestral_shift_proven":False},
    "Digitalis canariensis":{"guild":"effective_opportunistic_birds",
       "archipelago":"Canary","evidence":"doi:10.3732/ajb.0800055",
       "evolutionary_ancestral_shift_proven":False},
    "Roussea simplex":{"guild":"effective_native_day_gecko",
       "archipelago":"Mauritius","evidence":"doi:10.1086/593050; Kew POWO",
       "evolutionary_ancestral_shift_proven":False,
       "day_gecko_potential_bird_function_additional_studies":True},
    "Ruizia blackburniana":{"guild":"effective_native_day_gecko",
       "archipelago":"Mauritius","evidence":"PMID:17262697; accepted synonym of Trochetia blackburniana",
       "evolutionary_ancestral_shift_proven":False},
}
# Only source taxonomic location filters (coarse, roomy): reject exotic island
# photographs outside the independently documented native archipelago; even
# IN-band photos may be gardens and cannot be declared wild.
REGION_BOX={
    "Galapagos":{"lat_min":-2.5,"lat_max":1.8,"lon_min":-92.5,"lon_max":-89},
    "Canary":{"lat_min":27.0,"lat_max":30.5,"lon_min":-19.2,"lon_max":-12.4},
    "Mauritius":{"lat_min":-21.0,"lat_max":-19.6,"lon_min":57.0,"lon_max":58.7},
}

def require(ok,reason):
    if not bool(ok):raise ValueError("ALAN_ENDMIC_EXPOSURE_REAL: "+str(reason))

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def source_external(native,gecko,genus):
    require(sha(native)==NATIVE_SOURCE_SHA and
            sha(gecko)==GECKO_SOURCE_SHA and
            sha(genus)==HEMIDACTYLUS_SOURCE_SHA,
            "Original native botanical-country or independent 2026-08 gecko source SHA changed")
    n=pd.read_csv(native,dtype={"island":str},low_memory=False)
    g=pd.read_csv(gecko,dtype={"island":str})
    h=pd.read_csv(genus,dtype={"island":str})
    require(len(n)==77647 and n.island.nunique()==562 and
            not n.duplicated(KEY).any() and
            "recorded_species__native_at_TDWG3" in n.columns and
            LIGHT in n.columns,
            "Original Kew regional native site-year/photo observation cohort changed")
    require({"island","year","latbin","lonbin","actual_source_taxon_family",
             "actual_GBIF_human_observation_records"}.issubset(g.columns) and
            not g.duplicated(KEY+["actual_source_taxon_family"]).any(),
            "Original gecko source family is not verified by island site year")
    require({"island","year","latbin","lonbin","original_genus",
             "original_source_records"}.issubset(h.columns) and
            not h.duplicated(KEY+["original_genus"]).any(),
            "Original gecko source genus is not verified by island site year")
    return n,g,h

def archipelago_filter(d,archipelago):
    if d.empty:return d.copy()
    bounds=REGION_BOX[archipelago]
    lon=(d.lonbin.astype(float)+.5)*.05-180
    lat=(d.latbin.astype(float)+.5)*.05-90
    return d.loc[lon.between(bounds["lon_min"],bounds["lon_max"]) &
                 lat.between(bounds["lat_min"],bounds["lat_max"])].copy()

def visit_support(original,gecko,genus):
    if original.empty:
        return {"photo_sites_with_gecko_record_same_cell_year":0,
                "original_islands_with_gecko_earlier_observation_years":0,
                "original_islands_with_Hemidactylus_genus_earlier_photo":0,
                "photo_sites_with_gecko_record_earlier_same_island":0,
                "photo_sites_with_Hemidactylus_genus_earlier_same_island":0}
    d=gecko.loc[gecko.actual_source_taxon_family.ne("Scincidae")].copy()
    d=d.groupby(KEY,sort=False).actual_GBIF_human_observation_records.sum().rename(
       "same_cell_gecko_photodocumentation").reset_index()
    j=original[KEY].merge(d,on=KEY,how="left",validate="one_to_one")
    w=original[["island","year"]].copy().drop_duplicates()
    first=gecko.loc[gecko.actual_source_taxon_family.ne("Scincidae")].groupby(
        "island").year.min()
    first_hemi=genus.loc[genus.original_genus.eq("Hemidactylus")].groupby(
        "island").year.min()
    w["earlier_gecko_documented"]=(
        w.island.map(first).lt(w.year).fillna(False))
    w["earlier_hemi_documented"]=(
        w.island.map(first_hemi).lt(w.year).fillna(False))
    v=original[["island","year"]].merge(w,on=["island","year"],
                                      how="left",validate="many_to_one")
    require(len(v)==len(original),"Original site-year shape changed by gecko source")
    return {
      "photo_sites_with_gecko_record_same_cell_year":int(j.same_cell_gecko_photodocumentation.fillna(0).gt(0).sum()),
      "original_islands_with_gecko_earlier_observation_years":int(w.loc[w.earlier_gecko_documented,"island"].nunique()),
      "original_islands_with_Hemidactylus_genus_earlier_photo":int(w.loc[w.earlier_hemi_documented,"island"].nunique()),
      "photo_sites_with_gecko_record_earlier_same_island":int(v.earlier_gecko_documented.sum()),
      "photo_sites_with_Hemidactylus_genus_earlier_same_island":int(v.earlier_hemi_documented.sum()),
      "source_gecko_2013_2025_DWCA_not_prior_to_2013":True,
      "zero_photographs_not_proof_of_local_gecko_absence":True,
      "Hemidactylus_genus_not_a_verified_invasive_population":True,
    }

def census(study,traits,pcnl,groups,native,gecko,genus):
    d,macro,valid,studied,author_trait=sources(study,traits,pcnl,groups)
    n,g,h=source_external(native,gecko,genus)
    require(set(CONFIRMED).issubset(set(SEEDS)),
            "High-evidence published plant candidate not pinned in PR30 source register")
    require(sha(pcnl)==SITE_SHA,"Source 33-year PCNL annual data drift")
    site=valid[KEY+[LIGHT]].copy()
    require(not site.duplicated(KEY).any(),"Real 85k site/year source duplicated")
    actual=d.loc[d.species.isin(CONFIRMED)].copy()
    actual["in_signed_original_source_photo_cohort"]=True
    actual=actual.merge(site,on=KEY,how="left",validate="many_to_one")
    rows=[]
    counts_by_guild={}
    for species,meta in CONFIRMED.items():
        q=actual.loc[actual.species.eq(species)].copy()
        sel=archipelago_filter(q,meta["archipelago"])
        lit=sel.loc[sel[LIGHT].notna()].copy()
        source_unique_macro=lit.merge(macro[KEY],
                   on=KEY,how="inner",validate="many_to_one")
        baseline=n.loc[n["recorded_species__native_at_TDWG3"].gt(0),
                       KEY+[LIGHT,"recorded_species__native_at_TDWG3"]].copy()
        lit_native=lit.merge(baseline,on=KEY,how="inner",validate="many_to_one",
                              suffixes=("_candidate","_native"))
        # Native comparator can contain all TDWG3-range native plants,
        # not proven island natives or wild source individuals.
        native_island_years=baseline.loc[baseline.island.isin(lit_native.island.unique())]
        pairs=[]
        for (island,year),c in lit_native.groupby(["island","year"],sort=True):
            ref=native_island_years.loc[
                native_island_years.island.eq(island)&
                native_island_years.year.eq(year)]
            if len(ref)<2 or ref[LIGHT].nunique()<2:continue
            weight=ref["recorded_species__native_at_TDWG3"].to_numpy(float)
            v=float(np.average(ref[LIGHT].to_numpy(float),weights=weight))
            cov=float(np.mean(c[f"{LIGHT}_candidate"].to_numpy(float)))
            pairs.append({"island":str(island),"year":int(year),
                          "difference":cov-v,
                          "candidate_source_sites":len(c)})
        evidence=visit_support(lit,g,h)
        row={
            "accepted_plant_name":species,"evidence_guild":meta["guild"],
            "literature_verified_native_archipelago":meta["archipelago"],
            "source_reference":meta["evidence"],
            "independent_ancestral_pollinator_switch_proven":False,
            "original_photo_IDs_all_GSHHG_islands":int(q.photo_ids.sum()),
            "original_photo_observation_UUIDs_all_GSHHG_islands":int(q.observation_uuids.sum()),
            "all_island_polygon_units_observed":int(q.island.nunique()),
            "all_island_photo_grid_years":len(q),
            "same_verified_native_archipelago_photo_grid_years":len(sel),
            "same_verified_native_archipelago_distinct_GSHHG_island_units":
                int(sel.island.nunique()),
            "same_verified_native_archipelago_photo_IDs":int(sel.photo_ids.sum()),
            "candidate_original_photo_sites_with_actual_prior_1992to2024_PCNL":len(lit),
            "candidate_original_islands_with_actual_prior_1992to2024_PCNL":
                int(lit.island.nunique()),
            "candidate_median_historical_prior_log1p_PCNL":
                float(lit[LIGHT].median()) if len(lit) else None,
            "candidate_prior_loglight_quantiles":{
                str(z):float(lit[LIGHT].quantile(z)) for z in (.1,.5,.9)
            } if len(lit) else {},
            "candidate_sites_in_562_original_within_island_fullflora_ALAN_cohort":
                len(source_unique_macro),
            "candidate_islands_in_562_original_ALAN_cohort":
                int(source_unique_macro.island.nunique()),
            "candidate_sites_co_located_with_any_TDWG3_region_native_observed_flora":
                len(lit_native),
            "source_candidate_island_years_with_comparable_native_regional_reference":
                len(pairs),
            "source_comparable_reference_islands":
                len(set(k["island"] for k in pairs)),
            "source_candidate_minus_TDWG3_regional_native_photo_reference_logPCNL_median":
                float(np.median([p["difference"] for p in pairs])) if pairs else None,
            "reference_comparison_is_not_true_island_native_or_source_effort_equal":True,
            "original_gecko_reptile_photo_documentation":evidence,
            "original_photo_points_cannot_assert_wild_native_individuals":True,
        }
        rows.append(row)
        counts_by_guild.setdefault(meta["guild"],{"species":0,"native_archipelago_sites_with_light":0,
                                                   "source_native_archipelago_island_units":set()})
        a=counts_by_guild[meta["guild"]]
        a["species"]+=1
        a["native_archipelago_sites_with_light"]+=len(lit)
        a["source_native_archipelago_island_units"].update(lit.island.unique())
    cohorts={k:{"literature_confirmed_species":v["species"],
                "observed_archipelago_source_photo_grid_years_with_historic_light":
                  v["native_archipelago_sites_with_light"],
                "native_archipelago_GSHHG_island_units_with_any_supported_photo":
                  len(v["source_native_archipelago_island_units"])}
             for k,v in counts_by_guild.items()}
    report={
      "schema":SCHEMA,
      "original_93m_iNaturalist_photo_registry_count":93209351,
      "original_WCVP_accepted_plant_siteyear_photo_replay_source":True,
      "original_33year_real_prior_PCNL_source_SHA256":SITE_SHA,
      "original_regional_TDWG3_native_proxy_source_SHA256":NATIVE_SOURCE_SHA,
      "original_frozen_2026_08_gecko_source_SHA256":GECKO_SOURCE_SHA,
      "source_gecko_genus_original_SHA256":HEMIDACTYLUS_SOURCE_SHA,
      "original_signed_study_effective_moth_plant_only_Ipomoea":True,
      "independent_evolutionary_ancestral_pollinator_shifts_verified":0,
      "literature_verified_functional_pollination_species_n":len(CONFIRMED),
      "original_original_photo_site_year_population":85318,
      "original_within_island_fullflora_exposure_cohort":77647,
      "original_within_island_fullflora_GSHHG_islands":562,
      "source_verifiable_photo_coverage_by_pollinator_evidence":cohorts,
      "source_candidate_species_exposure":rows,
      "counterfactual_TDWG3_native_reference_not_true_island_endemic_population":True,
      "GloBI_visitation_tier_not_used_as_effective_pollination_evidence":True,
      "no_plant_fitness_population_insect_abundance_or_predation_measured":True,
      "no_verified_Hemidactylus_alien_status_per_island":True,
      "no_critically_endangered_species_precise_coordinates_released":True,
      "no_effects_regressed_when_archipelagos_do_not_offer_paired_guilds":True}
    return report

def run(study,traits,pcnl,groups,native,gecko,genus,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite verified original source macro census")
    d=census(study,traits,pcnl,groups,native,gecko,genus)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(d,sort_keys=True,indent=2,default=str)+"\n")
    print("ALAN_TRUE_ISLAND_SPECIALIZED_POLLLINATION_LIGHT_GECKO_EXPOSURE",
       json.dumps({"cohorts":d["source_verifiable_photo_coverage_by_pollinator_evidence"],
          "candidate_species":d["source_candidate_species_exposure"],
          "ancestral_shifts_verified":d["independent_evolutionary_ancestral_pollinator_shifts_verified"]},
       sort_keys=True),flush=True)
    return d

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("study","traits","pcnl","groups","native","gecko","genus","out"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    run(a.study,a.traits,a.pcnl,a.groups,a.native,a.gecko,a.genus,a.out)
