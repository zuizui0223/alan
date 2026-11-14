#!/usr/bin/env python3
"""Verify primary published fruit/seed effect table and source guild in Marques.

Original publication: Marques et al. 2015, Anais da Academia Brasileira de
Ciencias, DOI as printed 10.590/0001-3765201520140340. Original paper:
https://www.scielo.br/j/aabc/a/cSzbHtNkrVKjms8RK76zHMh/?lang=en

The original paper contains field records of bat Glossophaga soricina
and hummingbirds at the same research forest area and publishes actual
day/night/open reproductive fruit/seed treatments. However, its METHODS
does not separately establish the visitor observation year. Original 2010
study-year tag is from a secondary source, not an independently certified
original contemporaneous visitor-event date.

The values from published Table I are descriptive source facts, NOT a
new ecological discovery. This is a second, independent PUBLISHED PAPER
with direct reproductive evidence and field visitor taxon identities,
not a second complete same-year visitor-population join.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

ORIGINAL_BLOB = "39a82d966afcdde670603add1c881ace8791ffe9"
STUDY_ID = "Marques_2015"
SOURCE_DOI_AS_PRINTED = "10.590/0001-3765201520140340"
SOURCE_URL = "https://www.scielo.br/j/aabc/a/cSzbHtNkrVKjms8RK76zHMh/?lang=en"
PUB_FRUIT_COUNTS = {
    "day pollination": (81,808),
    "night pollination": (647,788),
    "open pollination": (327,420),
}
PUB_SEEDS_PER_FRUIT = {
    "day pollination": 3.8,
    "night pollination": 139.6,
    "open pollination": 73.9,
}


def require(x,msg):
    if not bool(x):
        raise ValueError("ALAN_MARQUES_PRIMARY_EFFECT_TABLE: "+str(msg))


def gitblob(path):
    b=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()


def verify(original_author,output):
    output=Path(output)
    require(not output.exists(),"Cannot overwrite previous source report")
    require(gitblob(original_author)==ORIGINAL_BLOB,
            "Original 1652 author treatment Git object changed")
    df=pd.read_csv(original_author,dtype=str,keep_default_na=False)
    require(len(df)==1652 and df.study_ID.nunique()==140,
            "Original global author study source changed")
    a=df.loc[df.study_ID.eq(STUDY_ID)].copy()
    require(len(a)==6 and a.plant_species.nunique()==1,
            "Original published primary table 2 metrics × 3 arms source changed")
    require(set(a.treatment_effectiveness_metric)=={"fruit set","seed set"},
            "Original fruit/seed source missing")
    require(set(a.treatment_condition)==set(PUB_FRUIT_COUNTS),
            "Original day/night/open treatment source changed")
    require(a.year_start.eq("2010").all() and
            a.Pollinator_observations_included.eq("y").all(),
            "Original publication metadata year/visitor indication changed")
    require(a.coordinates_lat.nunique()==1 and a.coordinates_lon.nunique()==1,
            "Original source group changed site geography")
    checks={}
    for t,(fruits,flowers) in PUB_FRUIT_COUNTS.items():
        b=a.loc[a.treatment_effectiveness_metric.eq("fruit set") &
                a.treatment_condition.eq(t)]
        require(len(b)==1,"Original primary fruit-set arm absent")
        raw=float(b.iloc[0].effectiveness_value)
        orig=fruits/flowers
        require(int(b.iloc[0].sample_size)==flowers and abs(raw-orig)<1e-6,
                "Published original Table I and global source fruit treatment disagree")
        checks[t]={"published_fruits":fruits,"published_flowers":flowers,
                   "published_fruitset_fraction":float(orig),
                   "global_author_matched_fraction":raw}
    seed={}
    for t,expected in PUB_SEEDS_PER_FRUIT.items():
        b=a.loc[a.treatment_effectiveness_metric.eq("seed set") &
                a.treatment_condition.eq(t)]
        require(len(b)==1 and int(b.iloc[0].sample_size)==50 and
                abs(float(b.iloc[0].effectiveness_value)-expected)<1e-8,
                "Published Table I seeds-per-fruit original value differs")
        seed[t]=expected
    v=checks["night pollination"]["published_fruitset_fraction"]
    day=checks["day pollination"]["published_fruitset_fraction"]
    result={
        "schema":"alan_primary_Marques2015_original_global_reproduction_taxon_source_v1",
        "status":"AUTHENTIC_PUBLIC_PRIMARY_FRUIT_SEED_TABLE_VERIFIED_STUDY_FIELD_VISITORS_YEAR_NOT_PROVEN",
        "provenance":{
            "original_1652_treatment_git_blob":ORIGINAL_BLOB,
            "author_article_DOI_as_printed_not_crossref_asserted":SOURCE_DOI_AS_PRINTED,
            "original_published_primary_HTML":SOURCE_URL,
            "original_global_reproductive_study_key":STUDY_ID,
            "global_compiled_study_year_not_primary_independently_checked":2010,
            "source_global_reproductive_metric_pairs":2,
            "original_source_same_unique_botanical_taxon":True,
            "study_level_field_visitor_observation_flag":True,
        },
        "original_primary_table_I_fruit_set":checks,
        "original_primary_table_I_seeds_per_fruit":seed,
        "functional_source_descriptive":{
            "night_divided_by_sum_night_day_original_fruit_fraction":
                float(v/(v+day)),
            "published_reproductive_fruit_set_night_minus_day":float(v-day),
            "night_reproductive_fruit_set_exceeds_day":bool(v>day),
            "primary_article_65h_field_visitor_observations":True,
            "primary_article_documented_bat_visitor":"Glossophaga soricina",
            "primary_article_documented_hummingbird_visitors":[
                "Thalurania glaucopis","Phaethornis pretrei"],
            "published_study_area_visitor_and_trial_locations_shared":True,
            "visitor_field_year_explicitly_proved_against_secondary_year_2010":False,
            "pollen_contact_reported_in_primary_article":True,
            "causal_effect_of_specific_animal_taxon_not_isolated_by_diel_bagging":True,
        },
        "limits":{
            "published_original_study_result_not_an_ALAN_discovery":True,
            "original_study_publication_not_independent_replication_of_its_own_fruit_table":True,
            "same_study_area_NOT_equivalent_to_site_year_visitation_event":True,
            "published_bat_hummingbird_identity_not_from_GloBI_species_level_proxy":True,
            "day_vs_night_exclusion_fitness_does_not_singlehandedly_isolate_bat_vs_bird_effect":True,
            "DOI_prefix_as_printed_nonstandard_not_automatically_renewed":True,
            "no_original_botanical_field_coordinates_republished":True,
        },
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_AUTHENTIC_MARQUES_PRIMARY_FIELD_BAT_HUMMINGBIRD_EFFECT_TABLE",
          json.dumps(result,ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--author",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    x=a.parse_args()
    verify(x.author,x.output)
