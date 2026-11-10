#!/usr/bin/env python3
"""Audit real exclusion-experiment evidence against the ALAN WCVP accepted flora.

External independent functional evidence:
Kendall & Nicholson (Ecology Letters, 2025, DOI:10.1111/ele.70036)
Original authors' public GitHub `liamkendall/bumble-bat`:
 - data/DvN_LiteratureDataClean_2024-09-13.csv
 - data/SppTraitData_clean.csv

Crucial limitation: these studies exclude pollinators by TIME OF DAY, not by
moth family. A nocturnal pollen contribution may involve other insects, birds,
bats and plant selfing; do not assign moth specificity from day/night trials.
Original raw source contains 145 species names; published final meta-analysis
reports 139 species after analysis preparation. No raw treatment effect
differences are silently equated to effective nocturnal dependence.

Compute strict WCVP accepted-species exact-name overlap, original 562 GSHHG
island presence and within-island ALAN geographic support. Output metadata
only, NO unreviewed binary moth-dependence labels or photographic RGB.
"""
from __future__ import annotations
import argparse, hashlib, json, sqlite3
from pathlib import Path
import pandas as pd
import numpy as np

TRAIT_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
STUDY_BLOB="39a82d966afcdde670603add1c881ace8791ffe9"
ALAN_ACCEPTED_SHA="198057897839bfa0657a3c3c3d9ce2ed4904041c350e05d5c61ced0dd5f5cd3d"
ALAN_PCNL_SHA="440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c"
KEY=["island","latbin","lonbin","year"]
SCHEMA="ALAN_diel_exclusion_experiment_plant_function_source_overlap_v1"


def require(ok,why):
    if not ok:raise ValueError("ALAN_TRUE_DIEL_EXPERIMENT_SOURCE: "+str(why))


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):h.update(part)
    return h.hexdigest()


def git_blob_sha(path):
    source=Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(source)}\0".encode()+source).hexdigest()


def normalize_name(series):
    return series.fillna("").astype(str).str.split().str.join(" ")


def experimental_records(trait_file,study_file):
    require(git_blob_sha(trait_file)==TRAIT_BLOB,
            "Original 2025 experimental plant-trait evidence Git blob changed")
    require(git_blob_sha(study_file)==STUDY_BLOB,
            "Original 2025 experimental day/night exclusion study data Git blob changed")
    trait=pd.read_csv(trait_file,dtype=str,keep_default_na=False)
    raw=pd.read_csv(study_file,dtype=str,keep_default_na=False)
    require(len(trait)==142 and len(raw)==1652,
            "Original published supplementary source row count changed")
    needed={"plant_species","plant_family","lifespan","breeding_system"}
    require(needed.issubset(trait.columns),
            "Independent experimental plant trait source lacks core fields")
    assert {"study_ID","plant_species","treatment_condition",
            "treatment_effectiveness_metric","effectiveness_value",
            "sample_size"}.issubset(raw.columns)
    trait["species"]=normalize_name(trait.plant_species)
    raw["species"]=normalize_name(raw.plant_species)
    require(trait.species.ne("").all() and raw.species.ne("").all(),
            "An experimental source plant name is unresolved")
    require(trait.species.is_unique,"Original trait species records not unique")
    # One experimental study×species×benefit metric has genuine day AND night
    # treatment arms; it is an evidence opportunity, NOT an effect direction.
    group_cols=["study_ID","species","treatment_effectiveness_metric"]
    groups=raw.groupby(group_cols,observed=True).treatment_condition.agg(
        lambda v:tuple(sorted(set(v))))
    both=groups.map(lambda z:{"day pollination","night pollination"}<=z)
    require(bool(both.all()) and len(groups)==262,
            "Original 262 paired day/night study×plant×metric groups incomplete")
    study_count=raw.groupby("species",observed=True).study_ID.nunique().rename(
        "n_experimental_studies")
    names=pd.DataFrame({"species":sorted(set(raw.species)|set(trait.species))})
    names=names.merge(trait[["species","lifespan","breeding_system",
                            "plant_family","bloom_period_simple","life_form"]],
                      on="species",how="left",validate="one_to_one")
    names=names.merge(study_count,on="species",how="left",validate="one_to_one")
    names["has_day_night_exclusion_study"]=names.n_experimental_studies.notna()
    names["has_trait_record"]=names.plant_family.notna()
    count=raw.treatment_condition.value_counts().to_dict()
    evidence={
       "published_meta_analysis_DOI":"10.1111/ele.70036",
       "published_final_meta_analysis_species":139,
       "original_clean_trait_rows":int(len(trait)),
       "original_clean_study_rows":int(len(raw)),
       "original_distinct_study_plant_names":int(raw.species.nunique()),
       "original_distinct_trait_plant_names":int(trait.species.nunique()),
       "original_raw_studies":int(raw.study_ID.nunique()),
       "raw_study_species_benefit_metric_day_night_pairs":int(both.sum()),
       "original_treatment_arms":{k:int(v) for k,v in count.items()},
       "moth_pollinator_taxa_not_identified_by_exclusion":True,
       "night_vs_day_treatments_not_a_binary_moth_dependence_label":True,
       "original_phenotype_trait_file_git_blob":TRAIT_BLOB,
       "original_day_night_study_csv_git_blob":STUDY_BLOB
    }
    return names,evidence


def source_match(accepted_csv,sites_csv,groups,experiment):
    require(sha256(accepted_csv)==ALAN_ACCEPTED_SHA,
            "Original 14,330 accepted WCVP plant source prevalence data changed")
    require(sha256(sites_csv)==ALAN_PCNL_SHA,
            "Original 33-year PCNL true photo site/year source changed")
    accepted=pd.read_csv(accepted_csv)
    sites=pd.read_csv(sites_csv)
    require(len(accepted)==14330 and accepted.species.is_unique,
            "Original 14,330 photo-observed WCVP accepted species changed")
    require(len(sites)==85318,"Signed original PCNL site-year denominator changed")
    accepted["species"]=normalize_name(accepted.species)
    sites["island"]=sites.island.astype(str)
    sites=sites.loc[sites.historic_ntl_complete.astype(str).str.lower().eq("true") &
          sites.pcnl_years_observed.eq(sites.pcnl_years_expected)].copy()
    s=sites.groupby(["island","year"],observed=True)[
        "pcnl_1992_to_previous_year_mean_log1p"].agg(["size","min","max"])
    keep=set(s.loc[s["size"].ge(2)&(s["max"]-s["min"]).ge(.02)].reset_index().island)
    sites=sites.loc[sites.island.isin(keep)].copy()
    require(len(sites)==77647 and sites.island.nunique()==562,
            "Exact 562 original island study reference changed")
    experiment["species"]=normalize_name(experiment.species)
    joined=experiment.merge(accepted,on="species",how="left",
                            validate="one_to_one",indicator=True)
    joined["source_matched_WCVP_accepted"] = joined._merge.eq("both")
    joined=joined.drop(columns="_merge")
    matched_names=sorted(set(joined.loc[joined.source_matched_WCVP_accepted,"species"]))
    acc=[]
    for group in range(16):
        path=Path(groups)/f"allflora-island-group-{group:02d}.sqlite"
        require(path.is_file(),f"Missing original 93m photo source group {group}")
        with sqlite3.connect(path) as con:
            for start in range(0,len(matched_names),60):
                names=matched_names[start:start+60]
                if not names:continue
                placeholders=",".join("?" for _ in names)
                q=f"""SELECT island,latbin,lonbin,year,species,
                       moth_positive_photo_ids,butterfly_positive_photo_ids
                       FROM flora WHERE species IN ({placeholders})"""
                acc.append(pd.read_sql_query(q,con,params=names))
    found=pd.concat(acc,ignore_index=True) if acc else pd.DataFrame(
        columns=KEY+["species","moth_positive_photo_ids","butterfly_positive_photo_ids"])
    found["island"]=found.island.astype(str)
    found=found.merge(sites[KEY],on=KEY,how="inner",validate="many_to_one")
    found["gloBI_moth_positive"]=found.moth_positive_photo_ids.gt(0)
    found["gloBI_butterfly_positive"]=found.butterfly_positive_photo_ids.gt(0)
    duplicate=found.duplicated(KEY+["species"],keep=False)
    if duplicate.any():
        n=found.loc[duplicate].groupby(KEY+["species"],observed=True)[
            ["gloBI_moth_positive","gloBI_butterfly_positive"]].nunique()
        require(n.to_numpy().max()<=1,"Original photo groups disagree on GloBI positive claim")
    found=found.drop_duplicates(KEY+["species"])
    support=found.groupby("species",observed=True).agg(
        n_islands=("island","nunique"),
        n_source_site_years=("year","size"),
        n_GloBI_moth_positive=("gloBI_moth_positive","sum"),
        n_GloBI_butterfly_positive=("gloBI_butterfly_positive","sum")).reset_index()
    return joined.merge(support,on="species",how="left",validate="one_to_one"),int(found.island.nunique())


def run(study,traits,accepted,sites,groups,out):
    candidates,meta=experimental_records(traits,study)
    joined,n_distinct_islands=source_match(accepted,sites,groups,candidates)
    mask=joined.source_matched_WCVP_accepted & joined.has_day_night_exclusion_study
    matched=joined.loc[mask]
    matched=matched.loc[matched.n_source_site_years.notna()]
    result={
        "schema":SCHEMA,"experiment_source":meta,
        "ALAN_photo_original_source_IDs":93209351,
        "ALAN_reference_islands":562,
        "ALAN_recorded_accepted_plant_species":14330,
        "experimental_raw_names_total":len(joined),
        "experimental_raw_study_names_exact_in_WCVP_accepted":int(mask.sum()),
        "experimental_raw_study_names_with_genuine_ALAN_island_photo_sites":len(matched),
        "experimental_raw_study_matched_distinct_islands":int(n_distinct_islands),
        "experimental_raw_study_matched_species_island_incidence_sum":int(
            pd.to_numeric(matched.n_islands,errors="coerce").sum()),
        "experimental_source_photo_site_species_units":int(
            matched.n_source_site_years.fillna(0).sum()),
        "matched_source_trait_rows":int(matched.has_trait_record.sum()),
        "matched_GloBI_moth_positive_source_species":int(
            matched.n_GloBI_moth_positive.fillna(0).gt(0).sum()),
        "matched_GloBI_butterfly_positive_source_species":int(
            matched.n_GloBI_butterfly_positive.fillna(0).gt(0).sum()),
        "no_distinct_moth_function_declared":True,
        "not_a_fitted_ALAN_response_if_support_limited":True,
        "night_pollination_exclusion_cannot_assign_moth_identity":True,
        "published_139_species_not_identical_to_145_raw_species_entries":True,
        "exact_WCVP_name_matches_only_synonyms_not_auto_reconciled":True
    }
    folder=Path(out)
    require(not folder.exists(),"Do not overwrite an independent experimental source")
    folder.mkdir(parents=True)
    fields=["species","has_day_night_exclusion_study","has_trait_record",
            "source_matched_WCVP_accepted","n_experimental_studies",
            "n_islands","n_source_site_years",
            "n_GloBI_moth_positive","n_GloBI_butterfly_positive"]
    joined[fields].to_csv(folder/"diel_experiment_accepted_plant_source_overlap.csv",
                          index=False)
    (folder/"diel_exclusion_experiment_ALAN_source_overlap.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ALAN_DIEL_139PUBLISHED_FUNCTION_SOURCE_MATCH",
          json.dumps(result,sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    for k in ("study","traits","accepted","sites","groups","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    run(a.study,a.traits,a.accepted,a.sites,a.groups,a.out)


if __name__=="__main__":main()
