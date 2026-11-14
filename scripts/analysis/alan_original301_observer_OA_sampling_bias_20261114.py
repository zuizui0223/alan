#!/usr/bin/env python3
"""Source-neutral 301 original reproductive studies: visitor-flag/OA selection audit.

Study-level field pollinator observations flag = some observations included
in author paper, NOT known moth guild at study site. Published OA metadata
is accessibility to paper text, NOT observation of a pollinator.

Compare original real night/(night+day) fruit-seed field outcome distributions
for (92 yes-observer-flag papers vs 26 others) and (36 OA paper candidates vs
56 other YES-flag studies). One independent publication contributes one study
median. Source source SHA256 and original 92-study DOI queue are independently
verified. Descriptive 3000 study-bootstrap intervals, NOT causal inference.
Only global aggregate details exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA = "21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af"
AUTHOR_GITBLOB = "39a82d966afcdde670603add1c881ace8791ffe9"


def require(condition,message):
    if not bool(condition):
        raise ValueError("ALAN_301_SOURCE_VISITOR_SELECTIVITY: "+str(message))


def gitblob(path):
    v=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(v)).encode()+b"\0"+v).hexdigest()


def difference(group1,group2,seed=20261114):
    a=np.asarray(group1,dtype=float)
    b=np.asarray(group2,dtype=float)
    require(len(a)>=20 and len(b)>=20 and
            np.all(np.isfinite(a)) and np.all(np.isfinite(b)),
            "Insufficient genuine independent published study groups")
    delta=float(a.mean()-b.mean())
    rng=np.random.default_rng(seed)
    x=np.mean(a[rng.integers(0,len(a),(3000,len(a)))],axis=1)
    y=np.mean(b[rng.integers(0,len(b),(3000,len(b)))],axis=1)
    return {
        "original_independent_studies_groupA":int(len(a)),
        "original_independent_studies_groupB":int(len(b)),
        "groupA_mean_study_median_night_share":float(a.mean()),
        "groupB_mean_study_median_night_share":float(b.mean()),
        "groupA_minus_B_mean_study_median_difference":delta,
        "descriptive_study_bootstrap_95CI":list(map(float,np.quantile(x-y,[.025,.975]))),
        "study_block_bootstrap_draws":3000,
        "published_studies_nonrandom_sample":True,
    }


def run(preflight,author,oa_queue,output):
    output=Path(output)
    require(not output.exists(),"Refuse to overwrite original source audit")
    require(hashlib.sha256(Path(preflight).read_bytes()).hexdigest()==SOURCE_SHA,
            "Original prePCNL field source changed")
    require(gitblob(author)==AUTHOR_GITBLOB,
            "Actual original 1652 author rows changed")
    source=pd.read_csv(preflight,dtype={"study_key":str,"species":str})
    require((len(source),source.study_key.nunique(),source.species.nunique())==
            (301,118,123),"Original global 301 field study cohort changed")
    raw=pd.read_csv(author,dtype=str,keep_default_na=False)
    require(len(raw)==1652 and raw.study_ID.nunique()==140,
            "Original author study rows altered")
    queue=json.loads(Path(oa_queue).read_text())
    require(queue["schema"]=="ALAN_92_real_study_public_bibliography_curation_queue_v1",
            "Original DOI+OA queue schema changed")
    qs=queue["records"]
    require(len(qs)==92 and len({x["original_study_reference"] for x in qs})==92,
            "Original public paper OA queue not 92 disjoint papers")
    oa={x["original_study_reference"]:bool(x.get("openalex_oa_record_claim"))
        for x in qs}
    oa_link={x["original_study_reference"]:bool(
            x.get("openalex_candidate_oa_location_url")) for x in qs}
    require(sum(oa.values())==36 and sum(oa_link.values())==36,
            "Original metadata OA candidate source changed")
    flags={}
    for key, group in raw.groupby("study_ID",observed=True):
        s=set(group.Pollinator_observations_included)
        require(s<=set(["y","n","NA"]), "Original author observation flag changed")
        flags[str(key)]="yes" if "y" in s else (
            "no_flag" if s=={"n"} else "source_unknown")
    source["source_paper_visitor_observation_flag"]=source.study_key.map(flags)
    require(source.source_paper_visitor_observation_flag.notna().all(),
            "Original paper visitor metadata missing")
    paper=source.groupby("study_key",observed=True).agg(
        night_share_median=("night_share","median"),
        year_median=("year_exp","median"),
        abs_lat_median=("lat",lambda x:float(np.median(abs(x)))),
        originals=("source_pair_id","size"),
        source_species=("species","nunique"),
    ).reset_index()
    paper["flag"]=paper.study_key.map(flags)
    paper["OA"]=paper.study_key.map(oa)
    paper["OA_link"]=paper.study_key.map(oa_link)
    obs=paper.loc[paper.flag.eq("yes")].copy()
    no=paper.loc[paper.flag.ne("yes")].copy()
    require((len(obs),len(no))==(92,26),
            "Original 92 observer flag source papers vs 26 others changed")
    require((int(obs.originals.sum()),int(no.originals.sum()))==(247,54),
            "Original 247 observer flag source pairs vs 54 others changed")
    require(obs.OA.notna().all() and no.OA.isna().all(),
            "OA source queue not restricted to visitor flag papers")
    withoa=obs.loc[obs.OA.eq(True)]
    withoutoa=obs.loc[obs.OA.eq(False)]
    require((len(withoa),len(withoutoa))==(36,56),
            "Actual 36 OA candidate papers vs 56 other flagged source changed")
    audit={
        "schema":"alan_source301_observer_reporting_and_OA_metaaccess_bias_audit_v1",
        "status":"REAL_301_SOURCE_OBSERVATION_METADATA_SAMPLING_BIAS_DIAGNOSTIC",
        "provenance":{
            "original_global_301_outcomes_sha256":SOURCE_SHA,
            "original_original_author_1652_gitblob":AUTHOR_GITBLOB,
            "original_reproductive_case_years":[1993,2022],
            "DOI_source_92_observation_publications":92,
            "public_OA_link_candidates":36,
        },
        "representation":{
            "independent_original_studies_total":int(len(paper)),
            "observer_flag_yes_studies":len(obs),
            "observer_flag_yes_original_paired_reproductive_outcomes":int(obs.originals.sum()),
            "observer_flag_no_or_unknown_studies":len(no),
            "observer_flag_no_or_unknown_original_paired_outcomes":int(no.originals.sum()),
            "observed_study_OA_candidate":len(withoa),
            "observed_study_NO_OA_candidate":len(withoutoa),
            "observed_study_OA_candidate_source_outcomes":int(withoa.originals.sum()),
            "observed_study_NO_OA_candidate_source_outcomes":int(withoutoa.originals.sum()),
        },
        "observer_yes_vs_other_papers_study_equal":{
            "night_share":difference(obs.night_share_median,no.night_share_median),
            "original_experiment_year_median":{
                "observer_yes":float(obs.year_median.median()),
                "other":float(no.year_median.median()),
            },
            "original_abs_latitude_median":{
                "observer_yes":float(obs.abs_lat_median.median()),
                "other":float(no.abs_lat_median.median()),
            },
        },
        "OA_link_yes_vs_no_within_92_observer_papers_study_equal":{
            "night_share":difference(withoa.night_share_median,withoutoa.night_share_median),
            "original_experiment_year_median":{
                "OA_candidate":float(withoa.year_median.median()),
                "other_observer_papers":float(withoutoa.year_median.median()),
            },
            "original_abs_latitude_median":{
                "OA_candidate":float(withoa.abs_lat_median.median()),
                "other_observer_papers":float(withoutoa.abs_lat_median.median()),
            },
        },
        "guardrails":{
            "no_flag_in_author_database_not_local_biological_absence":True,
            "OA_link_metadata_not_a_source_verified_moth_or_visitor":True,
            "study_medians_not_301_independent_studies":True,
            "study_bootstrap_descriptive_uncertainty_not_phylogeny":True,
            "study_reporting_status_and_OA_not_random_assignments":True,
            "paper_availability_not_measured_biological_response":True,
            "original_study_species_location_rows_NOT_exported":True,
            "no_ALAN_or_guild_effect_inferred":True,
        },
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(audit,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ACTUAL_ORIGINAL_REPRO_VISITOR_OA_SELECTION",
          json.dumps(audit,ensure_ascii=False),flush=True)
    return audit


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("preflight","author","oa_queue","output"):
        parser.add_argument("--"+k,type=Path,required=True)
    x=parser.parse_args()
    run(x.preflight,x.author,x.oa_queue,x.output)
