#!/usr/bin/env python3
"""Real author-published diel-exclusion reproductive comparisons x prior PCNL on 562 islands.

This is source-reconciled GLOBAL ECOLOGICAL OBSERVATION, not causal ALAN effects,
moth-only function, or site-specific experimentally measured pollination.

No missing GloBI or experiment record is treated as a true biological absence.
Original experiment treatment comparisons are made ONLY within the exact
study, species, outcome metric, measurement units, and original day/night arms.
The author study is not reanalyzed as a significance meta-analysis.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, LIGHT, load_sites, load_labels, original_accepted_species,
    site_outcome_model, species_model, require
)

STUDY_BLOB="39a82d966afcdde670603add1c881ace8791ffe9"
TRAIT_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
REPRO_METRICS=("fruit set","seed set")
MAIN_STUDY_PAIRS=262
ACCEPTED_SOURCE=14330
N_SOURCE_ROWS=2012925
N_SITES=77647
N_ISLANDS=562


def blob_hash(path):
    a=Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(a)}\0".encode()+a).hexdigest()


def canon(s):
    return s.fillna("").astype(str).str.split().str.join(" ")


def experimental_daynight(study_path:Path, trait_path:Path):
    require(blob_hash(study_path)==STUDY_BLOB,"Original author study blob drift")
    require(blob_hash(trait_path)==TRAIT_BLOB,"Original author trait blob drift")
    s=pd.read_csv(study_path,dtype=str,keep_default_na=False)
    t=pd.read_csv(trait_path,dtype=str,keep_default_na=False)
    require((len(s),len(t))==(1652,142),"Original author source row counts changed")
    required={"study_ID","plant_species","treatment_effectiveness_metric",
              "treatment_condition","effectiveness_value","effectiveness_units",
              "day_pollination_hrs","night_pollination_hrs"}
    require(required.issubset(s.columns),"Author treatment columns changed")
    require({"plant_species","bloom_period_simple","breeding_system"}.issubset(t.columns),
            "Author trait columns changed")
    s["species"]=canon(s.plant_species)
    t["species"]=canon(t.plant_species)
    require(t.species.is_unique,"Non-unique author plant trait")
    grouped=s.groupby(["study_ID","species","treatment_effectiveness_metric"],
                      observed=True,sort=False)
    n_compare=0
    rows=[]
    for (study,species,metric),g in grouped:
        night=g[g.treatment_condition.eq("night pollination")]
        day=g[g.treatment_condition.eq("day pollination")]
        if night.empty or day.empty:continue
        n_compare+=1
        # Keep ORIGINAL comparable response units within paired study metric.
        n_units=set(night.effectiveness_units)-{"","NA"}
        d_units=set(day.effectiveness_units)-{"","NA"}
        if len(n_units)>1 or len(d_units)>1 or (n_units and d_units and n_units!=d_units):
            continue
        def group_mean(x):
            vals=pd.to_numeric(x.effectiveness_value,errors="coerce")
            if vals.isna().any() or (vals<0).any():return None
            weights=pd.to_numeric(x.sample_size,errors="coerce") if "sample_size" in x else None
            if weights is not None and weights.notna().all() and (weights>0).all():
                return float(np.average(vals.to_numpy(float),weights=weights.to_numpy(float)))
            return float(vals.mean())
        a,b=group_mean(night),group_mean(day)
        if a is None or b is None or a+b<=0:continue
        day_hours=pd.to_numeric(g.day_pollination_hrs,errors="coerce")
        night_hours=pd.to_numeric(g.night_pollination_hrs,errors="coerce")
        hour_ok=(day_hours.notna().all() and night_hours.notna().all()
                 and (day_hours>0).all() and (night_hours>0).all())
        day_hr=float(day_hours.median()) if hour_ok else None
        night_hr=float(night_hours.median()) if hour_ok else None
        per_hour=(a/night_hr)/((a/night_hr)+(b/day_hr)) if hour_ok else None
        rows.append({"study":study,"species":species,"metric":metric,
                     "night_share":a/(a+b),"night_share_per_hour":per_hour})
    require(n_compare==MAIN_STUDY_PAIRS,"Original 262 experimental day/night pair units changed")
    raw=pd.DataFrame(rows)
    require(len(raw)>=190,"Too few valid original comparable experimental endpoints")

    # Primary is DIRECT reproductive outcome, not floral visitation, color,
    # body-pollen carriage or a source-labelled moth guild.
    repro=raw[raw.metric.isin(REPRO_METRICS)].copy()
    require(repro.study.nunique()>=75 and repro.species.nunique()>=100,
            "Insufficient reproductive day/night paired evidence")
    by_study=(repro.groupby(["study","species"],observed=True)
                    .agg(night_share=("night_share","mean"),
                         night_share_per_hour=("night_share_per_hour","mean"),
                         metrics=("metric","nunique")).reset_index())
    effect=(by_study.groupby("species",observed=True)
                  .agg(n_studies=("study","nunique"),
                       raw_reproductive_pairs=("metrics","sum"),
                       night_share=("night_share","median"),
                       night_share_per_hour=("night_share_per_hour","median"))
                  .reset_index())
    effect["night_stronger"]=effect.night_share.gt(.5+1.e-12)
    effect["day_stronger"]=effect.night_share.lt(.5-1.e-12)
    effect["exact_tie"]=~(effect.night_stronger|effect.day_stronger)
    effect["perhour_night_stronger"]=effect.night_share_per_hour.gt(.5+1.e-12)
    effect["perhour_day_stronger"]=effect.night_share_per_hour.lt(.5-1.e-12)
    effect["perhour_available"]=effect.night_share_per_hour.notna()
    # Study-trait flowering time and breeding system describe the published
    # study species, NOT local island populations nor autonomous selfing.
    trait=t[["species","bloom_period_simple","breeding_system","plant_family"]].copy()
    trait["bloom_period_simple"]=trait.bloom_period_simple.replace(
        {"NA":pd.NA,"":pd.NA}).astype("string")
    effect=effect.merge(trait,on="species",how="left",validate="one_to_one")
    effect["night_bloom"]=effect.bloom_period_simple.eq("night").fillna(False)
    effect["day_bloom"]=effect.bloom_period_simple.eq("day").fillna(False)
    effect["both_bloom"]=effect.bloom_period_simple.eq("both").fillna(False)
    require((effect[["night_stronger","day_stronger","exact_tie"]].sum(axis=1)==1).all(),
            "Experimental reproductive ranking not a partition")
    cohort_counts={
      "original_source_study_rows":len(s),
      "original_day_night_study_species_metric_pairs":n_compare,
      "comparable_numeric_pairs_all_metrics":len(raw),
      "comparable_reproductive_pairs":len(repro),
      "reproductive_species_with_study_comparisons":len(effect),
      "reproductive_independent_studies":int(by_study.study.nunique()),
      "reproductive_nocturnal_higher_species":int(effect.night_stronger.sum()),
      "reproductive_diurnal_higher_species":int(effect.day_stronger.sum()),
      "reproductive_equal_species":int(effect.exact_tie.sum()),
      "reproductive_perhour_supported_species":int(effect.perhour_available.sum()),
      "original_authors_trait_names":len(t),
      "source_treatment_effectiveness_metrics":sorted(set(raw.metric)),
      "interpretation":"Day/night removal-experiment reproductive endpoints; not identified moth species or species-site pollination",
    }
    return effect,cohort_counts


def source_regression(d,sites,flag,*,min_positive_islands=30,min_positive_sites=150):
    pos=d.loc[d[flag],KEY].drop_duplicates()
    n_islands=int(pos.island.nunique());n_sites=int(len(pos))
    sample={"islands_with_positive_taxa":n_islands,
            "site_years_with_positive_taxa":n_sites,
            "island_years_with_positive_taxa":int(pos[["island","year"]].drop_duplicates().shape[0]),
            "positive_species_site_years":int(d[flag].sum())}
    if n_islands<min_positive_islands or n_sites<min_positive_sites:
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_ISLAND_SUPPORT","source_support":sample}
    k=d.groupby(KEY,observed=True)[flag].sum().reset_index(name="n")
    s=sites.merge(k,on=KEY,how="left",validate="one_to_one")
    require(s.n.notna().all(),"Signed photo-site source incomplete")
    require(s.n.le(s.all_accepted_recorded_species).all(),"Outcome exceeds total source documented plants")
    fit=site_outcome_model(s,s.n.to_numpy(float))
    return {"status":"EXPLORATORY_DOCUMENTED_REPRODUCTIVE_SOURCE_COMPOSITION",
            "fit":fit, "source_support":sample}


def signed_contrast(d,sites,a,b,min_positive_islands=25):
    def count(flag):
        return d.loc[d[flag]].groupby(KEY,observed=True).size().reindex(
            pd.MultiIndex.from_frame(sites[KEY]),fill_value=0).to_numpy(int)
    aa,bb=count(a),count(b)
    ia=int(sites.loc[aa>0,"island"].nunique())
    ib=int(sites.loc[bb>0,"island"].nunique())
    if min(ia,ib)<min_positive_islands:
        return {"status":"NOT_ESTIMABLE_SIGNED_CONTRAST_INDEPENDENT_ISLAND_SUPPORT",
                "night_positive_islands":ia,"day_positive_islands":ib}
    f=site_outcome_model(sites,aa-bb)
    return {"status":"EXPLORATORY_PAIRED_SOURCE_LABEL_CONTRAST",
            "model":f,"night_positive_islands":ia,"day_positive_islands":ib,
            "note":"Same denominator and model; exact difference of slopes, not independent comparisons of p values"}


def calculate(sites_path,globisource,groups,study,traits,output):
    if output.exists():raise FileExistsError(output)
    sites=load_sites(sites_path)
    require(len(sites)==N_SITES and sites.island.nunique()==N_ISLANDS,"562 island source changed")
    labels=load_labels(globisource)
    d=original_accepted_species(groups,sites,labels)
    require(len(d)==N_SOURCE_ROWS and d.species.nunique()==ACCEPTED_SOURCE,"Original photographic plant records changed")
    exposure,cohort=experimental_daynight(study,traits)
    m=d.merge(exposure,on="species",how="left",validate="many_to_one")
    require(len(m)==len(d),"Changed original plant observation denominators")
    flags=("night_stronger","day_stronger","exact_tie","perhour_night_stronger",
           "perhour_day_stronger","perhour_available","night_bloom","day_bloom","both_bloom")
    for k in flags:m[k]=m[k].fillna(False).astype(bool)
    m["reproductive_experiment_documented"]=m.night_share.notna()
    m["bloom_documented"]=m.night_bloom|m.day_bloom|m.both_bloom
    m["night_stronger_moth_source"]=m.night_stronger&m.moth_original
    m["day_stronger_moth_source"]=m.day_stronger&m.moth_original

    support={}
    for k in ["night_stronger","day_stronger","exact_tie","reproductive_experiment_documented",
              "night_bloom","day_bloom","both_bloom","bloom_documented",
              "night_stronger_moth_source","day_stronger_moth_source",
              "perhour_night_stronger","perhour_day_stronger","perhour_available"]:
        support[k]={"experimental_species":int(exposure[k].sum()) if k in exposure else None,
                    "photo_observed_accepted_species":int(m.loc[m[k],"species"].nunique())}
    candidate_models={}
    for k in support:candidate_models[k]=source_regression(m,sites,k)
    direct=signed_contrast(m,sites,"night_stronger","day_stronger")
    perhour=signed_contrast(m,sites,"perhour_night_stronger","perhour_day_stronger")
    blossom=signed_contrast(m,sites,"night_bloom","day_bloom")
    # Cross-product support: do we have independent island×year source
    # sites with both experimentally night- and day-favored species?
    k=["island","year","latbin","lonbin"]
    s1=set(map(tuple,m.loc[m.night_stronger,k].drop_duplicates().to_numpy()))
    s2=set(map(tuple,m.loc[m.day_stronger,k].drop_duplicates().to_numpy()))
    overlap=s1 & s2
    same_genus=(exposure.assign(genus=exposure.species.str.split().str[0])
                .groupby("genus").agg(night=("night_stronger","sum"),
                                      day=("day_stronger","sum")))
    result={
     "schema":"ALAN_562islands_original93m_daynight_reproductive_exclusion_vs_priorPCNL_20261113_v1",
     "status":"ACTUAL_AUTHOR_PAIR_AND_ORIGINAL_ISLAND_PCNL_SOURCE_EXECUTED_EXPLORATORY",
     "source":{
       "study":"Kendall and Nicholson 2025 original public author data, DOI 10.1111/ele.70036",
       "author_study_original_gitblob":STUDY_BLOB,
       "author_trait_original_gitblob":TRAIT_BLOB,
       "photo_source":"original signed 93,209,351 photo ID index, 16 original accepted species site SQLite groups",
       "observed_site_years":len(sites),"islands":int(sites.island.nunique()),
       "observed_accepted_species":int(d.species.nunique()),
       "photo_species_site_years":len(d),
       "nightlight":"previous-year-only mean log1p Li-Cao corrected PCNL (1992 through photo-year-1)",
     },
     "study_source":cohort,
     "observational_label_source_support":support,
     "models":candidate_models,
     "direct_within_identical_photo_site_response_contrasts":{
       "night_vs_day_experimental_reproductive_strength":direct,
       "night_vs_day_perhour_sensitivity":perhour,
       "night_vs_day_bloom_period_source":blossom,
     },
     "within_source_photo_cell_year_support":{
       "photo_site_years_containing_both_night_and_day_stronger_species":len(overlap),
       "distinct_islands_with_both_source_groups_same_photo_cell_year":len(set(t[0] for t in overlap)),
       "observed_plant_genera_with_night_and_day_experimental_species":int(((same_genus.night>0)&(same_genus.day>0)).sum()),
     },
     "no_causal_or_moth_claim":{
       "night_day_exclusion_not_specific_to_moths":True,
       "species_mean_outside_experiment_geography_not_local_pollen_deposition":True,
       "old_GloBI_moth_positive_record_is_not_effective_pollination":True,
       "GloBI_absence_not_species_absence":True,
       "documented_cohort_is_study_selected_and_cannot_estimate_global_night_pollination_prevalence":True,
       "past_mean_pcNL_is_not_measured_33year_reproductive_trend":True,
       "individual_study_effects_not_independently_inferentially_tested":True,
       "multiple_unadjusted_posthoc_comparisons":True,
       "no_study_species_photo_location_rows_exported":True,
     }
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ACTUAL_DIEL_EXCLUSION_562_ISLAND_PCNL_RESULTS",json.dumps({
       "study":cohort,"groups":support,"contrast":result["direct_within_identical_photo_site_response_contrasts"],
       "main_models":{k:candidate_models[k] for k in ["night_stronger","day_stronger","reproductive_experiment_documented","night_bloom","day_bloom"]},
       "overlap":result["within_source_photo_cell_year_support"],
    },ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("sites","globi","groups","study","traits","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    calculate(a.sites,a.globi,a.groups,a.study,a.traits,a.out)
