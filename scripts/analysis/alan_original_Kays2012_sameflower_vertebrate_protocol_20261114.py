#!/usr/bin/env python3
"""Original Kays 2012 balsa same-flower / same-field-season direct guild audit.

Source facts: Kays et al. (2012), Mesoamericana 16(3):56-70, pp. 57,59,65,
official primary PDF https://www.revistas.up.ac.pa/index.php/mesoamericana/article/download/907/777 .
Gamboa original site from p57 and treatment / visited flower outcomes p65.
We do not claim randomly assigned ALAN or exact timestamped visitation for
the entire 301-study source. This explicitly catches the original author's
'complete exclusion' TREATMENT label despite bee-permeable original wire mesh.

Original global source 1652 author treatment rows Git blob and 301
source-prelight SHA are locked. No original site/plant row exported,
only independent paper-level aggregate inference and geometric distances.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from scripts.analysis.alan_original_diel_fourarm_reproductive_function_20261114 import (
    AUTHOR_STUDY_BLOB, exact_source_arms, qualify, descriptive_study_sign
)

SOURCE_301_SHA = "21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af"
STUDY_KEY = "Kays_2012"
REFERENCE = "https://www.revistas.up.ac.pa/index.php/mesoamericana/article/download/907/777"
SOURCE_KAYS_ARMS = {
    "day pollination": (13, 1/13),
    "night pollination": (6, 3/6),
    "open pollination": (6, 4/6),
    "complete exclusion": (13, 0.0),
}
# Published primary p57, Gamboa flower-experiment towers, DMS.
GAMBOA_LAT = 9 + 7/60 + 40.8/3600
GAMBOA_LON = -(79 + 41/60 + 48.8394/3600)
PIPELINE_LAT = 9 + 7/60 + 57/3600
PIPELINE_LON = -(79 + 43/60 + 34.32/3600)
BCI_LAT = 9 + 9/60 + 59.0394/3600
BCI_LON = -(79 + 50/60 + 12.8394/3600)


def require(x, msg):
    if not bool(x):
        raise ValueError("ALAN_KAYS_PRIMARY_SOURCE_PROTOCOL_AUDIT: " + msg)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def gitblob(path):
    a=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(a)).encode()+b"\0"+a).hexdigest()


def haversine_km(la,lo,lb,lob):
    a,b,c,d=np.deg2rad([la,lo,lb,lob])
    x=np.sin((c-a)/2)**2+np.cos(a)*np.cos(c)*np.sin((d-b)/2)**2
    return float(12742*np.arctan2(np.sqrt(x),np.sqrt(1-x)))


def summary_study_sign(a, metric):
    """Paper-medianeffects; any sub-cohorts use only study grain."""
    study=a.groupby("study_key",observed=True)[metric].median()
    nneg=int(study.lt(-1e-12).sum())
    npos=int(study.gt(1e-12).sum())
    return {
        "papers":int(len(study)),
        "plant_taxa":int(a.species.nunique()),
        "original_matched_site_year_metric_outcomes":int(len(a)),
        "study_median":float(study.median()),
        "paper_medians_negative":nneg,
        "paper_medians_positive":npos,
        "paper_medians_tie":int(len(study)-npos-nneg),
        "nominal_two_sided_sign_p":float(binomtest(nneg,nneg+npos,.5).pvalue),
    }


def two_endpoint_holm(a):
    """Whole-publication sensitivity, no new cherry-picked hypothesis test."""
    p={k:summary_study_sign(a,k)["nominal_two_sided_sign_p"]
       for k in ("two_window_nonadditivity","open_relative_to_best")}
    ordered=sorted(p, key=p.get)
    holm={}
    carry=0.
    for i,k in enumerate(ordered):
        carry=max(carry,min(1., (len(ordered)-i)*p[k]))
        holm[k]=float(carry)
    return {
        "nominal":p,
        "two_new_endpoints_Holm":holm,
        "source_study_omission_sensitivity_not_independent_confirmatory_evidence":True
    }


def audit(author, preflight, out):
    out=Path(out)
    require(not out.exists(), "Refuse overwrite of original source verification")
    require(gitblob(author)==AUTHOR_STUDY_BLOB,"Original 1652-row author Git blob changed")
    require(sha256(preflight)==SOURCE_301_SHA,"Original pre-light 301 experiment source changed")
    source=pd.read_csv(preflight,dtype={"study_key":str,"species":str})
    require((len(source),source.study_key.nunique(),source.species.nunique())==
            (301,118,123),"301/118/123 global source no longer original")
    a=pd.read_csv(author,dtype=str,keep_default_na=False)
    require((len(a),a.study_ID.nunique())==(1652,140),
            "Original published experiment-author treatment rows changed")
    k=a.loc[a.study_ID.eq(STUDY_KEY)]
    require(len(k)==4 and k.treatment_effectiveness_metric.eq("fruit set").all(),
            "Original one-field Kays source arms changed")
    require(set(k.treatment_condition)==set(SOURCE_KAYS_ARMS),
            "Kays original published treatment categories changed")
    require(k.year_start.eq("2010").all() and k.month_start.eq("January").all()
            and k.month_end.eq("February").all(),
            "Original meta-analysis field season Jan-Feb 2010 not preserved")
    designs=set(k.study_design_comments)
    require(len(designs)==1 and
            "Bees readily passed through these cages" in next(iter(designs)) and
            "wire mesh cage" in next(iter(designs)),
            "Original author directly reports bee-permeable exclusion")
    for name,(n,v) in SOURCE_KAYS_ARMS.items():
        x=k.loc[k.treatment_condition.eq(name)]
        require(len(x)==1 and int(x.iloc[0].sample_size)==n and
                np.isclose(float(x.iloc[0].effectiveness_value),v,atol=1e-12),
                "Author Kays paper Table 3 numeric treatment mismatch: "+name)
    _,metadata=exact_source_arms(author,source)
    # The same original 153 four-arm published site-year responses as PR #58.
    full,eligibility=qualify(_)
    require((len(full),full.study_key.nunique(),full.species.nunique())==
            (153,79,85),"Original published 153/79/85 4-arm source changed")
    km=full.loc[full.study_key.eq(STUDY_KEY)]
    require(len(km)==1 and km.metric.iloc[0]=="fruit set",
            "Original Kays 2010 source is not unique in four-arm dataset")
    row=km.iloc[0]
    for col,raw in (("day",1/13),("night",3/6),("open",4/6),("closed",0.)):
        require(np.isclose(row[col],raw,atol=1e-12),
                "Source four-arm original Kays contrast reconstructed wrongly "+col)
    d_gam=haversine_km(float(row.lat),float(row.lon),GAMBOA_LAT,GAMBOA_LON)
    d_pipe=haversine_km(float(row.lat),float(row.lon),PIPELINE_LAT,PIPELINE_LON)
    d_bci=haversine_km(float(row.lat),float(row.lon),BCI_LAT,BCI_LON)
    require(d_gam<.05 and d_pipe>3 and d_bci>15 and row.year_exp==2010,
            "Published Gamboa original source cannot be spatially / temporally joined")
    require(row.two_window_nonadditivity>0,
            "Published balsa original 4-arm index should be positive")
    pre={
        metric:summary_study_sign(full,metric)
        for metric in ("two_window_nonadditivity","open_relative_to_best","open_vs_closed")
    }
    require(pre["two_window_nonadditivity"]["paper_medians_negative"]==49 and
            pre["two_window_nonadditivity"]["papers"]==79 and
            np.isclose(pre["two_window_nonadditivity"]["study_median"],
                       -.057470648581547865,atol=1e-10),
            "Frozen original four-arm positive-control baseline does not reproduce")
    clean=full.loc[~full.study_key.eq(STUDY_KEY)].copy()
    post={
        metric:summary_study_sign(clean,metric)
        for metric in ("two_window_nonadditivity","open_relative_to_best","open_vs_closed")
    }
    result={
      "schema":"ALAN_Kays2012_original_sameflower_vertebrate_pollination_protocol_check_v1",
      "status":"SOURCE_VERIFIED_GAMBOA_2010_SAME_FLOWER_VISIT_FRUIT_ORIGINAL_4ARM",
      "original_source":{
        "author_treatment_rows":1652,"source_pre_light_real_pairs":301,
        "source_authentic_four_arm_pairs":153,
        "independent_fourarm_papers":79,
        "published_author_blob":AUTHOR_STUDY_BLOB,
        "source_original_sha256":SOURCE_301_SHA,
        "source_original_FCP_PCNLP_unchanged":True,
      },
      "independently_published_primary_factual_protocol":{
        "published_article_url":REFERENCE,
        "original_author_study_ID":STUDY_KEY,
        "original_meta_analysis_trial_year":2010,
        "original_meta_analysis_months":"January-February",
        "primary_reported_study_flowering_season":"2009-2010",
        "primary_same_study_season_support":"2009-10 spans source Jan-Feb 2010",
        "original_source_Gamboa_distance_km":d_gam,
        "separation_from_other_published_tree_Pipeline_km":d_pipe,
        "separation_from_other_published_tree_BCI_km":d_bci,
        "original_source_treatment_numeric_fruits_out_of_flowers":{
            "day_only":[1,13], "night_only":[3,6], "unrestricted_open":[4,6],
            "wire_mesh_bees_allowed":[0,13]
        },
        "original_insect_passable_wire_mesh_rough_opening_cm":.9,
        "entire_vertebrate_and_large_moth_exclusion_NOT_all_insect_exclusion":True,
        "original_global_label_complete_exclusion_IS_semantically_incomplete":True,
        "original_primary_kinkajou_observed_on_each_of_three_fruiting_night_treatment_flowers":True,
        "primary_flower_visitor_obs_not_full_night_coverage":True,
        "vertebrate_exclusion_not_individual_animal_taxon_manipulation":True,
        "same_experiment_diel_treatment_not_a_moth_specific_experiment":True,
      },
      "fourarm_result":{
        "original_single_balsa_I":float(row.two_window_nonadditivity),
        "original_single_balsa_J":float(row.open_relative_to_best),
        "original_single_balsa_K":float(row.open_vs_closed),
        "full_original_153_pairs":pre,
        "excluding_Kays_incompletely_closed_study":post,
        "full_new_endpoint_Holm":two_endpoint_holm(full),
        "after_exclusion_new_endpoint_Holm":two_endpoint_holm(clean),
        "one_Kays_fourarm_comparison_was_removed":True,
        "no_original_previous_results_rewritten":True,
      },
      "interpretation":{
        "visitor_same_field_site_trial_year_and_individual_fruit_tracking":True,
        "temporal_primary_visitors_span_original_study_season_not_exact_each_full_night":True,
        "only_one_additional_original_independent_paper_not_several_animal_species":True,
        "source_verifies_kinkajou_visitation_not_causal_taxon_specific_pollination":True,
        "original_primary_bee_passable_cage_is_not_true_closed_autogamy_control":True,
        "no_ALAN_manipulation_or_light_effect_from_this_source":True,
        "no_original_study_botanical_coord_year_raw_rows_exported":True,
      }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_KAYS_GAMBOA_2010_TRUE_SAME_FLOWER_DAYNIGHT_GUILD",
          json.dumps(result,ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("author","preflight","out"):
        p.add_argument("--"+n,type=Path,required=True)
    x=p.parse_args()
    audit(x.author,x.preflight,x.out)
