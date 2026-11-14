#!/usr/bin/env python3
"""Peer-reviewed 2019 ALAN positive reproduction original-source metadata gate.

This audit distinguishes a *published reproductive finding* from reproduction of
its original data. DOI 10.1002/ecs2.2550 and linked Zenodo 1743917 are queried
live; only the Zenodo metadata/manifest and cautious methods-and-outcomes ledger
are exported. No author raw flower/plant observational data is redistributed.
"""
from __future__ import annotations
import argparse
import datetime
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path

DOI="10.1002/ecs2.2550"
RECORD=1743917
URL=f"https://zenodo.org/api/records/{RECORD}"

PAPERS=[
 {"study_id":"Macgregor_2019_Ecosphere","doi":DOI,"publication_year":2019,
  "study_country":"United Kingdom","organism":"Silene latifolia",
  "research_group":"Macgregor_Pocock_Fox_Evans",
  "light_comparison":"full-night HPS+LED vs unlit; separate part-night vs unlit",
  "genuine_reported_outcome":"female plant flowers successfully pollinated",
  "direction_full_night_vs_unlit":"positive",
  "direction_part_night_vs_unlit":"not resolved",
  "seed_count_full_night_vs_unlit":"not resolved",
  "reproduction_observation_measured":True,
  "direct_moth_visitor_action_measured_in_same_flower":False,
  "night_day_access_experiment":True,
  "separate_experiment_block_count":"4 temporal runs; five plot positions",
  "same_plots_reused_across_light_treatments":True,
  "original_raw_numeric_extract_in_this_QA":False,
  "evidence_tier":"primary_peer_reviewed_published_direction_NOT_original_numeric_refit"},
 {"study_id":"Knop_2017_Nature","doi":"10.1038/nature23288",
  "publication_year":2017,"study_country":"Switzerland",
  "organism":"Cirsium oleraceum",
  "research_group":"Knop_Zoller_Ryser",
  "light_comparison":"streetlight vs unlit",
  "genuine_reported_outcome":"focal plant fruit-set",
  "direction_full_night_vs_unlit":"negative",
  "reported_relative_fruit_difference_percent":-13,
  "reproduction_observation_measured":True,
  "evidence_tier":"primary_peer_reviewed_published_direction_NOT_original_numeric_refit"},
 {"study_id":"Boom_2020_SciReports","doi":"10.1038/s41598-020-75471-1",
  "publication_year":2020,"study_country":"Netherlands",
  "organism":"Silene latifolia","research_group":"Boom_et_al",
  "light_comparison":"white green red vs dark",
  "genuine_reported_outcome":"fruit set, ovule fertilization, larval fruit infestation, estimated seed number",
  "direction_full_night_vs_unlit":"not a prespecified same-spectrum contrast",
  "reported_global_seed_ANOVA_p":0.63,
  "reproduction_observation_measured":True,
  "evidence_tier":"primary_peer_reviewed_aggregates_source_reproduced_in_ALAN_PR67"},
 {"study_id":"Li_2026_bioRxiv","doi":"10.64898/2026.07.10.737511",
  "publication_year":2026,"study_country":"United Kingdom",
  "organism":"three moth-attracting plant taxa, not verified exclusive nocturnal pollination",
  "research_group":"Li_Adeniji_Meah_Clements",
  "light_comparison":"26 lamp linear array vs unlit and 4 lamp point lights",
  "genuine_reported_outcome":"plant any seed production / total seed and capsule counts",
  "direction_full_night_vs_unlit":"compound array positive",
  "reproduction_observation_measured":True,
  "evidence_tier":"preprint_primary_79_record_numeric_refit_in_ALAN_PR74"},
 {"study_id":"Lu_2026_JPE","doi":"10.1093/jpe/rtag160",
  "publication_year":2026,"study_country":"China",
  "organism":"Pedicularis kansuensis / Euphrasia regelii",
  "research_group":"Lu_Wang_Cai_Chen",
  "light_comparison":"ALAN vs natural dark, night/day/open/closed access",
  "genuine_reported_outcome":"conspecific pollen on stigma, seed set",
  "direction_full_night_vs_unlit":"negative for nocturnal-dependent SI, unresolved change in autonomously selfing SC",
  "reported_relative_seed_difference_percent_nocturnal_dependent":-56,
  "reported_nocturnal_pollen_receipt_difference_percent":-53.8,
  "reproduction_observation_measured":True,
  "evidence_tier":"accepted_peer_reviewed_published_abstract_NO_original_experimental_rows"},
 {"study_id":"Lu_2026_EcolEvol","doi":"10.1002/ece3.74368",
  "publication_year":2026,"study_country":"China",
  "organism":"29 night-blooming plants source visiting network",
  "research_group":"Lu_Wang_Cai_Chen",
  "light_comparison":"8000K 4.52lux artificial streetlight vs unlit, five blocks, July-Aug 2024",
  "genuine_reported_outcome":"animal visiting counts NOT stigma pollen or fruit/seed",
  "reported_nocturnal_lepidoptera_visits_change_percent":-57.4,
  "reproduction_observation_measured":False,
  "evidence_tier":"primary_peer_reviewed_five_block_original_source_reproduced_ALAN_PR72"}
]

def require(cond,msg):
    if not cond: raise AssertionError(msg)

def fetch_manifest():
    req=urllib.request.Request(URL,headers={
         "User-Agent":"ALAN-GitHub-peer-reviewed-paper-source-audit/1.0",
         "Accept":"application/json"})
    try:
        with urllib.request.urlopen(req,timeout=45) as q:
            b=q.read(2_000_000)
            status=q.status
    except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
        return {"status":"API_UNAVAILABLE_NO_RAW_DATA_CLAIM",
                "http_status":getattr(e,"code",None),
                "error_type":type(e).__name__,
                "requested_url":URL,"metadata_trusted":False}
    try:
        j=json.loads(b.decode("utf-8"))
    except (UnicodeDecodeError,ValueError):
        return {"status":"NON_JSON_METADATA_UNVERIFIED","requested_url":URL,
                "http_status":status,"metadata_trusted":False}
    declared=str(j.get("id",""))
    if declared!=str(RECORD):
        return {"status":"RECORD_ID_MISMATCH_FAIL_CLOSED",
                "requested_url":URL,"id_returned":declared,"metadata_trusted":False}
    meta=j.get("metadata",{})
    files=j.get("files",{})
    if isinstance(files,dict):
        files=files.get("entries",{})
        files=[{"key":key,**v} for key,v in files.items()] if isinstance(files,dict) else []
    files=files if isinstance(files,list) else []
    manifest=[]
    for f in files:
        manifest.append({"name":f.get("key") or f.get("filename"),
                         "bytes":f.get("size"),
                         "checksum":f.get("checksum"),
                         "is_raw_data_candidate":str(f.get("key") or f.get("filename") or "").lower().endswith((".csv",".xlsx",".rds",".rdata",".tsv",".txt",".zip"))})
    return {"status":"OFFICIAL_ZENODO_METADATA_VERIFIED","requested_url":URL,
            "http_status":status,"metadata_trusted":True,
            "record_id":int(declared),
            "title":meta.get("title"),"doi":j.get("doi") or meta.get("doi"),
            "file_manifest":manifest,"n_files":len(manifest),
            "note":"Metadata alone cannot verify raw data nor reproduce published GLMM"}

def validate(studies):
    require(len(studies)==6,"Study count drift")
    ids=[s["study_id"] for s in studies]
    require(len(ids)==len(set(ids)),"Duplicated paper ID")
    require(len({s["doi"] for s in studies})==len(studies),"Duplicate article DOI")
    require(sum(s.get("reproduction_observation_measured",False) for s in studies)==5,
            "Outcome tier arithmetic changed")
    require(len({s["research_group"] for s in studies})==5,
            "Independence classes changed")
    require(studies[4]["research_group"]==studies[5]["research_group"],
            "Tibet same research-group caution lost")
    require(studies[0]["direction_full_night_vs_unlit"]=="positive","2019 original positive omitted")
    require(studies[1]["direction_full_night_vs_unlit"]=="negative","2017 negative contrast omitted")

def main(out):
    out=Path(out);require(not out.exists(),"Source audit overwrite forbidden")
    validate(PAPERS)
    result={"schema":"ALAN_2017_2026_peer_reviewed_and_preprint_light_reproductive_source_independence_v1",
       "date_utc":datetime.date.today().isoformat(),
       "source_manifest":fetch_manifest(),
       "verified_qualitative_evidence_papers":PAPERS,
       "number_publications":len(PAPERS),
       "source_reproduction_outcome_publications":sum(s["reproduction_observation_measured"] for s in PAPERS),
       "one_visit_only_publication_not_seed_effect":True,
       "original_three_2017_2019_2020_peer_reviewed_studies_independent":True,
       "first_valid_2019_positive_peer_reviewed_reproductive_counterexample":True,
       "independent_study_families":5,
       "critical_source_independence":{
          "Tibetan_JPE_and_EcolEvol_same_four_author_names":True,
          "EcolEvol_five_field_2024_Gannan_site_true":True,
          "JPE_exact_plot_year_moth_visit_overlap_NOT_VERIFIED":True,
          "do_not_pool_JPE_and_EcolEvol_as_independent_effects":True,
          "JPE_measured_pollen_seed_whereas_EcolEvol_measured_network_visits":True},
       "meta_analytic_stopping_rules":{
          "cannot_pool_different_endpoints_as_one_log_response_ratio":True,
          "cannot_use_direction_only_P_values_as_inverse_variance_effects":True,
          "cannot_treat_2019_four_temporal_runs_as_four_independent_research_groups":True,
          "no_Lu_two_paper_source_overlap_certification":True,
          "no_2019_new_effect_size_or_original_raw_reconstruction_if_Zenodo_only_code":True,
          "all_reports_are_associations_or_author_manipulation_not_ecological_mechanism_imputation":True},
       "next_estimate":"Obtain per-study matched plant reproduction numerator denominator, field-block independent means/SE, arm/period/sex and actual light before genuine meta-regression."}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("ALAN_MACGREGOR_2019_REAL_SOURCE_METADATA",json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    main(p.parse_args().out)
