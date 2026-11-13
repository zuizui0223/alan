#!/usr/bin/env python3
"""Original Kendall/Nicholson field-study VISITOR OBSERVATION coverage audit.

Crucial distinction:
  Pollinator_observations_included = y means a published paper contains some
  pollinator observations. It does NOT mean the curator database identifies
  moths, butterflies, pollination efficacy or guild at each fruit/seed trial.

We cross source study IDs only to 301 SHA-audited actual day/night field
reproductive experiment locations and source GloBI v0.11 positive-only guild
plants. A mechanism comment with 'moth'/'bat'/'bee' is a TEXT SCREEN,
NOT an independently verified local visitor observation.

Aggregate counts plus public DOI-only study screening queue; no DOI×taxon/site/year, coordinates or raw comments.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import pandas as pd

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, require
)
from scripts.analysis.alan_genuine_diel_GloBI_v011_flower_guild_overlap_20261113 import (
    original_v011_flags, TIERS
)

AUTHOR_BLOB="39a82d966afcdde670603add1c881ace8791ffe9"


def original_published_observations(path):
    b=Path(path).read_bytes()
    sha=hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()
    require(sha==AUTHOR_BLOB,"Original author study pollinator-observation meta Git blob changed")
    x=pd.read_csv(path,dtype=str,keep_default_na=False)
    require(len(x)==1652 and x.study_ID.nunique()==140,
            "Original 1652 source / 140 papers changed")
    required={"study_ID","plant_species","Pollinator_observations_included",
              "mechanism_comments","Study_comments","DOI"}
    require(required.issubset(x.columns),"Actual original pollinator observation field omitted")
    state={}
    for study,records in x.groupby("study_ID",observed=True):
        vals=set(records.Pollinator_observations_included)
        require(vals<={"y","n","NA"},"Unknown original observation source status")
        if "y" in vals:
            stat="included_some_observations"
        elif vals=={"n"}:
            stat="not_included_in_author_study_source"
        else:
            stat="unknown_or_mixed"
        # Original published DOI is public study-level bibliographic
        # metadata; no botanical study species or coordinates ever exported.
        doi=set(records.DOI)-{"NA",""}
        comment=" ".join(records.mechanism_comments.drop_duplicates().astype(str))
        state[study]={
            "status":stat,
            "n_original_treatment_rows":len(records),
            "source_paper_has_DOI":bool(doi),
            "source_paper_original_DOI_literal":next(iter(doi)) if len(doi)==1 else None,
            "source_paper_DOI_ambiguous":bool(len(doi)>1),
            # SCREENING KEYWORDS ONLY, never coded as actual observer taxa.
            "contains_moth_keyword_in_mechanism_text":bool(re.search(
                r"\b(?:moths?|hawk.?moths?|sphinx)\b",comment,re.I)),
            "contains_bee_keyword_in_mechanism_text":bool(re.search(
                r"\b(?:bees?|bumble.?bees?)\b",comment,re.I)),
            "contains_bat_keyword_in_mechanism_text":bool(re.search(
                r"\bbats?\b",comment,re.I)),
        }
    summary={}
    for name in ("included_some_observations","not_included_in_author_study_source",
                 "unknown_or_mixed"):
        v=[i for i in state.values() if i["status"]==name]
        summary[name]=len(v)
    require(summary=={
        "included_some_observations":111,
        "not_included_in_author_study_source":27,
        "unknown_or_mixed":2
    },"Original 140 published study visitor-observation flags changed")
    return state,summary


def run(joined,receipt,author_study,flags,flag_receipt,out,doi_queue=None):
    out=Path(out)
    if out.exists():raise FileExistsError("Existing independent source audit cannot be overwritten")
    source,_=load_source(joined,receipt)
    require((len(source),source.study_key.nunique(),source.species.nunique())==
            (301,118,123),"Actual field study source altered")
    state,pubsummary=original_published_observations(author_study)
    x=source.copy()
    x["source_obs_state"]=x.study_key.map(lambda s:state[s]["status"] if s in state else "")
    require(x.source_obs_state.ne("").all(),"Field study missing from original author publication")
    names,meta=original_v011_flags(flags,flag_receipt)
    x=x.merge(names[["plant_name","all_moth_any_family",
                    "all_butterfly_any_family",
                    "non_inat_trophich_moth_any_family",
                    "non_inat_trophich_butterfly_any_family"]],
              left_on="species",right_on="plant_name",
              how="left",validate="many_to_one")
    for col in ("all_moth_any_family","all_butterfly_any_family",
                "non_inat_trophich_moth_any_family",
                "non_inat_trophich_butterfly_any_family"):
        x[col]=x[col].fillna(0).astype(int)
    obs=x.source_obs_state.eq("included_some_observations")
    original_studies=set(x.study_key)
    study_source=[state[k] for k in original_studies]
    observed_paper=[a for a in study_source if a["status"]=="included_some_observations"]
    results={
      "schema":"ALAN_original301_DvN_author_study_paper_includes_pollinator_observations_audit_v1",
      "status":"AUTHENTIC_PUBLISHED_STUDY_POLLINATOR_OBSERVATION_AVAILABILITY_SCREEN_ONLY",
      "sources":{
         "original_author_publication_source_Gitblob":AUTHOR_BLOB,
         "original_140_published_studies_observation_status":pubsummary,
         "original_real_field_experiment_pairs":301,
         "original_real_field_experiment_publications":118,
         "original_real_field_species":123,
         "source_GloBI_full_24m_source_sha":meta["sha256"],
         "source_GloBI_field_event_time_and_georef_both_available":
             meta["georeferenced_and_dated_rows"],
      },
      "published_observation_coverage":{
         "field_source_studies_with_some_author_reported_pollinator_observations":
             int(sum(s["status"]=="included_some_observations" for s in study_source)),
         "field_source_studies_without_included_visitor_records_in_author_database":
             int(sum(s["status"]=="not_included_in_author_study_source" for s in study_source)),
         "field_source_studies_observer_status_unknown":
             int(sum(s["status"]=="unknown_or_mixed" for s in study_source)),
         "original_field_daynight_reproductive_pairs_in_studies_reporting_pollinator_observations":
             int(obs.sum()),
         "original_field_species_in_studies_reporting_pollinator_observations":
             int(x.loc[obs,"species"].nunique()),
         "original_study_DOI_available_in_yes_study_source":
             int(sum(s["source_paper_has_DOI"] for s in observed_paper)),
         "original_studies_with_observations_and_moth_keyword_in_MECHANISM_comments_only":
             int(sum(s["contains_moth_keyword_in_mechanism_text"] for s in observed_paper)),
         "original_studies_with_observations_and_bee_keyword_in_MECHANISM_comments_only":
             int(sum(s["contains_bee_keyword_in_mechanism_text"] for s in observed_paper)),
         "original_studies_with_observations_and_bat_keyword_in_MECHANISM_comments_only":
             int(sum(s["contains_bat_keyword_in_mechanism_text"] for s in observed_paper)),
      },
      "cross_source_audit":{
         "field_studies_with_original_observations_FLAG_AND_positive_moth_GloBI_recorded_plants":
             int(x.loc[obs & x.all_moth_any_family.eq(1),"study_key"].nunique()),
         "field_reproductive_pairs_with_original_observations_FLAG_AND_positive_moth_GloBI_recorded_plants":
             int((obs & x.all_moth_any_family.eq(1)).sum()),
         "field_species_with_original_observations_FLAG_AND_positive_moth_GloBI_recorded_plants":
             int(x.loc[obs & x.all_moth_any_family.eq(1),"species"].nunique()),
         "field_studies_with_original_observations_FLAG_AND_positive_butterfly_GloBI_recorded_plants":
             int(x.loc[obs & x.all_butterfly_any_family.eq(1),"study_key"].nunique()),
         "field_studies_with_original_observations_FLAG_AND_strict_moth_GloBI_recorded_plants":
             int(x.loc[obs & x.non_inat_trophich_moth_any_family.eq(1),"study_key"].nunique()),
      },
      "what_is_not_identified":{
         "original_observer_observation_inclusion_does_not_encode_visitor_taxonomic_guild":True,
         "comment_keyword_not_a_source_confirmed_visitor":True,
         "positive_GloBI_at_any_other_site_not_local_original_study_visitor":True,
         "no_moth_nectar_POLLEN_DEPOSITION_service_observed_from_flag":True,
         "source_DOI_screening_does_NOT_republish_original_rare_plant_locality":True,
         "not_native_island_specific_or_ALAN_causal":True,
         "next_curation_gate":"Read original author study DOI evidence / reported pollinator observation table, extract exact observed visitor guild, target plant and same field site/year; keep paper-level observation flag distinct from local guild counts"
      }
    }
    # Open bibliographic DOI queue to procure ACTUAL published visitor field
    # observations at the same study plant/site/time. Study-level metadata
    # ONLY; never publish the original field coordinates or plant identity.
    # This is a curation queue, not a database of verified local pollinators.
    if doi_queue is not None:
        queue=[]
        count_by_study=x.groupby("study_key",observed=True).size().to_dict()
        for study in sorted(original_studies):
            item=state[study]
            if item["status"]!="included_some_observations":continue
            queue.append({
                "original_publication_source_key":str(study),
                "original_source_DOI_literal":item["source_paper_original_DOI_literal"],
                "original_source_DOI_unique":not item["source_paper_DOI_ambiguous"],
                "author_says_paper_includes_pollinator_observations":True,
                "actual_qualified_night_day_reproductive_source_pairs":int(
                    count_by_study[study]),
                "moth_keyword_in_proposed_MECHANISM_notes_UNVERIFIED":
                    item["contains_moth_keyword_in_mechanism_text"],
                "bee_keyword_in_proposed_MECHANISM_notes_UNVERIFIED":
                    item["contains_bee_keyword_in_mechanism_text"],
                "bat_keyword_in_proposed_MECHANISM_notes_UNVERIFIED":
                    item["contains_bat_keyword_in_mechanism_text"],
                "actual_pollinator_guild_verified_in_this_queue":False,
                "same_study_site_year_visitor_record_exists_NOT_YET_CHECKED":True,
            })
        require(len(queue)==92 and sum(
            bool(a["original_source_DOI_literal"]) for a in queue)==91,
            "Original DOI-only 92-study curation opportunities changed")
        destination=Path(doi_queue)
        if destination.exists():raise FileExistsError("Do not replace existing DOI curation queue")
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_text(json.dumps({
            "schema":"ALAN_public_DvN_published_visitor_observation_study_DOI_curation_queue_v1",
            "status":"PUBLIC_BIBLIOGRAPHIC_STUDY_METADATA_NOT_VERIFIED_GUILD",
            "eligible_publications":len(queue),
            "source_papers_with_DOI":91,
            "source_original_301_reproductive_field_outcomes":301,
            "paper_flag_not_certain_site_matched_local_taxon_observation":True,
            "mechanism_comment_keywords_ARE_NOT_TRUE_VISITOR_EVIDENCE":True,
            "needed_evidence_fields":["same_original_field_site","same_year",
                    "same_plant_taxon","original_visitor_species_or_family",
                    "visit_event_time_or_activity_window",
                    "visit_observation_effort","effective_pollen_stigma_deposition"],
            "no_rare_botanical_name_or_source_field_coordinates":True,
            "study_DOI_queue":queue},indent=2,ensure_ascii=False,allow_nan=False)+"\n")
        results["source_metadata_acquisition_queue"]={
            "papers_original_observations_y":92,
            "papers_with_public_original_DOI":91,
            "paper_level_metadata_export_without_species_or_field_coordinates":True,
            "actual_local_animal_guilds_not_yet_curated":True
        }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(results,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ORIGINAL_STUDY_POLINATOR_OBSERVATION_SCREEN",
          json.dumps({key:results[key] for key in ("published_observation_coverage","cross_source_audit")},
                     ensure_ascii=False),flush=True)
    return results


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("joined","receipt","author_study","globi_flags","globi_receipt","out"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--paper_queue_out",type=Path,required=False)
    a=p.parse_args()
    run(a.joined,a.receipt,a.author_study,a.globi_flags,a.globi_receipt,a.out,a.paper_queue_out)
