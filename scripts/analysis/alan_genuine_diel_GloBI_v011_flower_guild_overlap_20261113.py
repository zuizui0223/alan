#!/usr/bin/env python3
"""True independent GloBI v0.11 source GUILD x actual published day/night
fruit/seed reproductive function: source-support and positivity audit.

Inputs are authenticated ORIGINAL sources:
  * 301 actual paired same-field-study×plant×site×year×season fruit/seed
    reproductive exclusion outcomes with MD5-authenticated previous-year
    1992–2021 PCNL; 118 studies and 123 source plant names.
  * ALL 24,453,954 original GloBI v0.11 interpreted interaction records,
    compressed SOURCE-verified moth/butterfly plant-positive flags; original
    candidate flower interactions 64,316, NOT observed nocturnal pollen.
  * Independent Kendall/Nicholson original botanical bloom-time source.

CRITICAL: GloBI v0.11 moth/butterfly flower positive evidence is from any
location and time, and has ZERO original qualified flower-interaction rows
with BOTH geocoordinates and event date. No local ecological guild or temporal
matching is identified. Missing positive guild source NEVER means true
nocturnal/diurnal pollinator absence. Butterfly and non-butterfly
Lepidoptera are broad documented clades, not individually nocturnal/diurnal
behavior.

No species-location/study original rows or restricted raw traits exported.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, LIGHT, require
)
from scripts.analysis.alan_flower_multitraits_reproductive_diel_transfer_20261113 import (
    original_traits
)

FLAGS_SHA="cf750192b4d641a2d904f4f74d5b87d41fc97e04822a67d513d9bad74ece8e01"
SOURCE_SHA="670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3"
AUTHOR_TRAITS_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
COHORT=(301,118,123)
TIERS={
   "all_source_v011":("all_moth_any_family","all_butterfly_any_family"),
   "non_inat_trophich_source":(
       "non_inat_trophich_moth_any_family",
       "non_inat_trophich_butterfly_any_family"),
   "non_inat_visits_source":(
       "non_inat_visits_moth_any_family",
       "non_inat_visits_butterfly_any_family"),
   "non_inat_trophich_insect_species_named":(
       "non_inat_trophich_moth_species_insect",
       "non_inat_trophich_butterfly_species_insect")
}
MIN_GROUP_SPECIES=20
MIN_GROUP_STUDIES=25


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for c in iter(lambda:stream.read(4<<20),b""):h.update(c)
    return h.hexdigest()


def original_v011_flags(csv:Path,receipt:Path):
    require(sha256(csv)==FLAGS_SHA,
            "Real GloBI full v0.11 positive flowering interaction taxon flags drift")
    meta=json.loads(Path(receipt).read_text())
    require(meta["status"]=="FULL_GLOBI_V011_SOURCE_COMPLETE" and
            meta["sha256"]==SOURCE_SHA and
            meta["all_interaction_rows"]==24453954 and
            meta["qualified_flower_interaction_rows"]==64316 and
            meta["plant_names_with_any_positive_claim"]==4547,
            "Not the original fully scanned source family-positive GloBI v0.11")
    require(meta["georeferenced_and_dated_rows"]=={"moth":0,"butterfly":0}
            and meta["flower_visits_not_effective_nocturnal_pollination"]
            and meta["missing_GloBI_not_pollinator_absence"],
            "Original positive guild event time+geo support changed")
    flags=pd.read_csv(csv,compression="gzip",dtype={"plant_name":str})
    require(len(flags)==4547 and flags.plant_name.is_unique,
            "Original positive taxon name flags duplicates or incomplete")
    for a,b in TIERS.values():
        require(a in flags.columns and b in flags.columns,
                "Cannot infer a pollinator family from source column absence")
        for col in (a,b):
            require(flags[col].isin([0,1]).all(),"GloBI source positive evidence invalid")
    return flags,meta


def summary_cells(x,flag):
    y=x.loc[flag]
    positive_names=y.species.unique()
    s=y[["species","night_share"]].groupby("species",observed=True).night_share.median()
    bloom=y.loc[y.bloom_period_simple.isin(("night","day"))]
    return {
       "original_true_paired_field_study_site_year_metrics":int(len(y)),
       "distinct_original_published_studies":int(y.study_key.nunique()),
       "distinct_original_plant_species":int(y.species.nunique()),
       "distinct_source_plant_families":int(y.plant_family.nunique()),
       "night_bloom_source_species_documented":int(bloom.loc[
           bloom.bloom_period_simple.eq("night"),"species"].nunique()),
       "day_bloom_source_species_documented":int(bloom.loc[
           bloom.bloom_period_simple.eq("day"),"species"].nunique()),
       "original_flowering_unknown_or_both_species":int(
           y.loc[~y.bloom_period_simple.isin(("night","day")),"species"].nunique()),
       "source_reproductive_night_share_median_across_original_pairs":(
           float(y.night_share.median()) if len(y) else None),
       "source_reproductive_night_share_median_across_source_plant_species":(
           float(s.median()) if len(s) else None),
       "historic_prior_light_zero_source_pairs":int(y[LIGHT].eq(0).sum()),
       "historic_prior_light_positive_source_pairs":int(y[LIGHT].gt(0).sum()),
       "same_site_time_moth_pollination_effectiveness_NOT_ascertained":True,
    }


def compare_positive_guild_documentation(x,scope,a,b):
    m=x[a].eq(1)
    d=x[b].eq(1)
    groups={
      "moth_positive_source_any_daytime_guild_unknown":m,
      "butterfly_positive_source_any_nighttime_guild_unknown":d,
      "both_moth_and_butterfly_positive":m&d,
      "moth_positive_butterfly_NOT_documented":m&~d,
      "butterfly_positive_moth_NOT_documented":d&~m,
      "no_positive_moth_butterfly_claim_UNKNOWN_BIOLOGICALLY":~(m|d)
    }
    tab={name:summary_cells(x,q) for name,q in groups.items()}
    only_moth=tab["moth_positive_butterfly_NOT_documented"]
    only_bfly=tab["butterfly_positive_moth_NOT_documented"]
    a_st=only_moth["distinct_original_published_studies"]
    b_st=only_bfly["distinct_original_published_studies"]
    aa=only_moth["distinct_original_plant_species"]
    bb=only_bfly["distinct_original_plant_species"]
    if (min(a_st,b_st)>=MIN_GROUP_STUDIES and
        min(aa,bb)>=MIN_GROUP_SPECIES):
        raise ValueError("New substantial independent source overlap! Must preregister model prior to inference.")
    return {
      "source_scope":scope,
      "observed_positive_guild_documentation_counts":tab,
      "can_infer_real_moth_vs_butterfly_original_reproductive_contrast":{
         "status":"NOT_ESTIMABLE_INDEPENDENT_SOURCE_PLANT_AND_PAPER_GUILD_SUPPORT",
         "minimum_required_unique_plant_species_PER_class":MIN_GROUP_SPECIES,
         "minimum_required_independent_studies_PER_class":MIN_GROUP_STUDIES,
         "source_moth_only_documented_plant_species":aa,
         "source_butterfly_only_documented_plant_species":bb,
         "source_moth_only_original_studies":a_st,
         "source_butterfly_only_original_studies":b_st,
         "not_a_p_value_and_no_negative_guild_assigned":True
      },
      "biological_limits":{
         "GloBI_source_positive_not_effective_pollination":True,
         "GloBI_source_moth_animal_may_be_day_active":True,
         "GloBI_source_butterfly_animal_may_be_night_active":True,
         "GloBI_butterfly_NOT_documented_does_NOT_mean_no_daytime_pollinators":True,
         "GloBI_moth_NOT_documented_does_NOT_mean_no_nighttime_pollinators":True,
         "source_plant_name_matching_exact_only_no_genus_inferred":True,
         "no_original_same_experimental_site_year_animal_event_matches":True
      }
    }


def run(joined:Path,receipt:Path,author_trait:Path,flags:Path,
        globi_receipt:Path,out:Path):
    if out.exists():raise FileExistsError("Scientific original aggregate output exists")
    source,original_receipt=load_source(joined,receipt)
    require((len(source),source.study_key.nunique(),source.species.nunique())==COHORT,
            "Original field experiment study population drift")
    flora=original_traits(author_trait)
    require(flora.species.is_unique,"Independent author study plant trait duplicated")
    flora=flora[["species","bloom_period_simple","plant_family"]].copy()
    flora.bloom_period_simple=flora.bloom_period_simple.replace(
        {"NA":"unknown","":"unknown"})
    flora.plant_family=flora.plant_family.replace(
        {"NA":"unknown","":"unknown"})
    x=source.merge(flora,on="species",how="left",validate="many_to_one")
    x.bloom_period_simple=x.bloom_period_simple.fillna("unknown")
    x.plant_family=x.plant_family.fillna("unknown")
    orig_positive,metadata=original_v011_flags(flags,globi_receipt)
    x=x.merge(orig_positive,left_on="species",right_on="plant_name",how="left",
              validate="many_to_one")
    for a,b in TIERS.values():
        x[a]=x[a].fillna(0).astype(int)
        x[b]=x[b].fillna(0).astype(int)
    require(len(x)==COHORT[0],"GloBI positive-only identity join changed field study denominators")
    source_hit=int(x.plant_name.notna().sum())
    source_hit_species=int(x.loc[x.plant_name.notna(),"species"].nunique())
    require((source_hit,source_hit_species)==(63,23),
            "Exact original 301 reproduction species and GloBI confirmed positive match drift")
    output={
        "schema":"ALAN_original301_reproductive_nightshare_GloBI_v011_true_visitor_source_support_v1",
        "status":"REAL_301_EXPERIMENT_GLOBI_24M_SOURCE_GUILD_OVERLAP_AUDITED_NOT_CAUSAL",
        "source":{
            "published_paired_reproduction_PCNLPrior_SHA256":original_receipt["joined_file_sha256"],
            "original_paired_day_night_reproductive_site_year_metrics":COHORT[0],
            "original_published_independent_studies":COHORT[1],
            "original_source_plant_species":COHORT[2],
            "experiment_outcome":"night/(night+day) genuine fruit/seed production at actual study original field site and year",
            "prior_satellite_PCNL_raster_years":[1992,2021],
            "original_full_GloBI_v011_all_interactions":24453954,
            "GloBI_v011_flower_interactions":64316,
            "GloBI_v011_plant_positive_names":4547,
            "GloBI_original_species_guild_flags_SHA256":FLAGS_SHA,
            "GloBI_source_SHA256":SOURCE_SHA,
            "real_source_plant_species_exactly_matched_visitor_positive_GloBI":source_hit_species,
            "real_source_reproductive_pairs_with_any_positive_animal_guild_GloBI":source_hit,
            "original_GloBI_source_georeferenced_AND_event_dated_moth_flower_rows":
                metadata["georeferenced_and_dated_rows"]["moth"],
            "original_GloBI_source_georeferenced_AND_event_dated_butterfly_flower_rows":
                metadata["georeferenced_and_dated_rows"]["butterfly"],
        },
        "guild_positive_source_tier_support":{
            name:compare_positive_guild_documentation(x,name,*columns)
            for name,columns in TIERS.items()},
        "strict_functional_conclusion":{
            "biological_guild_availability_at_reproduction_FIELD_SITE_YEAR_IDENTIFIED":False,
            "moth_vs_daytime_reproductive_service_contrast_identified":False,
            "functional_redundancy_of_daytime_pollinators_demonstrated":False,
            "positive_visitation_source_is_not_local_effective_pollen":True,
            "no_biological_absence_from_missing_source":True,
            "GloBI_and_experiment_study_can_share_published_sources_not_strictly_independent_measurement":True,
            "actual_day_night_function_source_is_EXPERIMENTALLY_PAIRED":True,
            "external_GloBI_guild_species_only_and_not_same_study_site_or_year":True,
            "outcome_and_orig_prior_light_not_co_labeled_for_positive_source_selection":True,
            "ALAN_not_experimentally_manipulated":True,
            "guild_specific_local_predator_shift_not_measured":True,
            "no_original_taxon_by_site_or_study_identifiers_exported":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_VERIFIED_GLOBAL_EXPERIMENTAL_DIEL_GLOBI_GUILD_SUPPORT",
          json.dumps({
            "source":output["source"],
            "positive_evidence":{
                scope:{
                  name:{"species":r["distinct_original_plant_species"],
                        "papers":r["distinct_original_published_studies"],
                        "pairs":r["original_true_paired_field_study_site_year_metrics"],
                        "species_nightshare_median":
                            r["source_reproductive_night_share_median_across_source_plant_species"]}
                  for name,r in obj["observed_positive_guild_documentation_counts"].items()
                  if name in ("moth_positive_source_any_daytime_guild_unknown",
                              "butterfly_positive_source_any_nighttime_guild_unknown",
                              "both_moth_and_butterfly_positive",
                              "moth_positive_butterfly_NOT_documented",
                              "butterfly_positive_moth_NOT_documented")}
                for scope,obj in output["guild_positive_source_tier_support"].items()
            }
          },ensure_ascii=False),flush=True)
    return output


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("joined","receipt","author_trait","globi_flags","globi_receipt","out"):
        p.add_argument("--"+arg,required=True,type=Path)
    s=p.parse_args()
    run(s.joined,s.receipt,s.author_trait,s.globi_flags,s.globi_receipt,s.out)
