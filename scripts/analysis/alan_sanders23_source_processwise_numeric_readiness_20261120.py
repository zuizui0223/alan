#!/usr/bin/env python3
"""Independent-source study-grain A/B quantitative readiness, not a pooled meta.

Genuine Sanders 2021 original 126-study / 1304-effect author source (original
CSV SHA256 d11d1e29...; Zenodo author file MD5 pinned in acquisition module).

For all 23 author-KEEP stage-A broad source-study IDs, JOIN original complete
field/lab and numeric mean/SD/n fields to a title-screen-only process
classification. Then count:
  (1) publication/source-study IDs with any records,
  (2) field/lab effect rows,
  (3) effect rows with positive original means and original SD/n,
  (4) original Study IDs containing these numerical candidates.
Original N can count insects/trap observations; no source field plot/lighting
assignment IDs or covariance were verified. These are upper-bound candidate
source papers, NOT admissible k for rma.mv or causal water/land meta-regression.

Do not export raw author Study IDs, individual effect estimates, ecological
observation rows or plot locations. Preserve positive, negative and null data.
"""
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path

from scripts.acquisition.alan_sanders2020_source_126paper_Dryad_preflight_20261119 import zenodo_original_author_mirror
from scripts.analysis.alan_sanders2020_stageAB_source_census_20261119 import original_source_rows
from scripts.analysis.alan_sanders2020_keep_delete_stageAB_primary_screen_20261119 import (
    audit_keep_delete, source_label_stage, original_arm_summaries_positive, doi_norm,
)
from scripts.analysis.alan_sanders23_original_stageA_measurement_process_triage_20261119 import (
    audit as validate_source_process_triage,
)

REGISTRY=Path("config/alan_sanders2020_23paper_stageA_original_endpoint_title_triage_20261119.json")
SOURCE_SHA="d11d1e29c05826f2b4a5c1c2d63b9c195dc2effbedd1ff531afa7b6bdf497693"
SCHEMA="alan_20261120_sanders_original_23paper_stageA_processwise_numeric_candidate_source_census_v1"


def process_quantitative_sufficiency(rows,original_identity,manual_triage,strict=True):
    """Classified source counts, never an experimental/meta-analysis eligibility decision."""
    if original_identity["original_CSV_SHA256"]!=SOURCE_SHA:
        raise ValueError("Original 1304-effect Sanders author CSV SHA changed")
    source=audit_keep_delete(rows,original_identity,strict=strict)
    bibliography=source["candidate_publication_bibliography_for_fulltext_screening"]["A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN"]
    validate_source_process_triage(manual_triage,bibliography,strict=strict)

    codebook={}
    for entry in manual_triage["source_study_entries"]:
        key=tuple(sorted(entry["original_source_DOIs"]))
        if key in codebook:
            raise ValueError("Same original paper DOI set duplicated in source process codebook")
        codebook[key]=entry

    a_keep=[r for r in rows if r["Delete"]=="Keep" and
            source_label_stage(r)["A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN"]]
    studies=defaultdict(list)
    for r in a_keep:
        studies[r["Study"]].append(r)
    if strict and (len(a_keep)!=152 or len(studies)!=23):
        raise ValueError("Original author 152 Keep / 23 source-study Stage A candidate size drift")
    paired=[]
    seen_process=Counter()
    mapped_keys=set()
    for study, source_rows in studies.items():
        dois=tuple(sorted({d for r in source_rows for d in [doi_norm(r.get("DOI"))] if d}))
        if not dois:
            raise ValueError("Original source Study lost all DOI information")
        entry=codebook.get(dois)
        if entry is None:
            raise ValueError("Unmatched author-original DOI cluster in manually titled 23-study review")
        if dois in mapped_keys:
            raise ValueError("Same DOI group assigned to two original source Study IDs")
        mapped_keys.add(dois)
        proc=entry["preliminary_outcome_family_for_screen_only"]
        habitat=entry["preliminary_original_habitat"]
        n_field=sum(r.get("Field.Lab")=="Field" for r in source_rows)
        n_lab=sum(r.get("Field.Lab")=="Laboratory" for r in source_rows)
        n_valid=sum(original_arm_summaries_positive(r) for r in source_rows)
        n_valid_field=sum(r.get("Field.Lab")=="Field" and original_arm_summaries_positive(r) for r in source_rows)
        n_valid_lab=sum(r.get("Field.Lab")=="Laboratory" and original_arm_summaries_positive(r) for r in source_rows)
        if n_field+n_lab!=len(source_rows) or n_valid_field+n_valid_lab!=n_valid:
            raise ValueError("Source field/lab numerical-eligibility accounting differs")
        paired.append({
            "class":proc,
            "habitat_source_hypothesis":habitat,
            "total_effect_rows":len(source_rows),
            "field_rows":n_field,"lab_rows":n_lab,
            "numeric_positive_mean_SD_n_rows":n_valid,
            "numeric_field_rows":n_valid_field,
            "numeric_laboratory_rows":n_valid_lab,
        })
        seen_process[proc]+=1
    if set(codebook)!=mapped_keys:
        raise ValueError("Not all manually classified original sources matched the Keep dataset")
    if strict and sum(p["field_rows"] for p in paired)!=130:
        raise ValueError("Original broad A author Keep field rows no longer 130")
    if strict and sum(p["lab_rows"] for p in paired)!=22:
        raise ValueError("Original broad A author Keep lab rows no longer 22")
    if strict and sum(p["numeric_positive_mean_SD_n_rows"] for p in paired)!=145:
        raise ValueError("Original numerical arm candidate rows no longer 145")

    counts={}
    for family in sorted(seen_process):
        group=[p for p in paired if p["class"]==family]
        counts[family]={
            "source_Study_ID_clusters_from_original_author_keep":len(group),
            "source_Keep_effect_rows":sum(x["total_effect_rows"] for x in group),
            "source_Keep_Field_effect_rows":sum(x["field_rows"] for x in group),
            "source_Keep_Laboratory_effect_rows":sum(x["lab_rows"] for x in group),
            "source_arithmetically_valid_positive_arm_means_SD_n_rows":sum(x["numeric_positive_mean_SD_n_rows"] for x in group),
            "source_Field_arithmetically_valid_positive_arm_means_SD_n_rows":sum(x["numeric_field_rows"] for x in group),
            "source_Laboratory_arithmetically_valid_positive_arm_means_SD_n_rows":sum(x["numeric_laboratory_rows"] for x in group),
            "source_Study_IDs_with_any_Field_candidate":sum(x["field_rows"]>0 for x in group),
            "source_Study_IDs_with_any_numeric_Field_candidate":sum(x["numeric_field_rows"]>0 for x in group),
            "source_habitat_text_hypotheses_from_initial_title_screen_only":dict(sorted(Counter(x["habitat_source_hypothesis"] for x in group).items())),
            "not_an_original_verified_same_endpoint_meta_sample":True,
            "not_a_confirmed_independent_field_illumination_programme_count":True
        }
    if strict and sum(v["source_Study_ID_clusters_from_original_author_keep"] for v in counts.values())!=23:
        raise ValueError("Process categorisation no longer partitions 23 authors' source Study IDs")
    b_result={}
    for label in ("B_flowering_phenology_NOT_reproduction","B_plant_seed_fruit_pollen_reproduction_CANDIDATE","B_invertebrate_possible_flower_visits_not_pollen"):
        keep=[r for r in rows if r["Delete"]=="Keep" and source_label_stage(r).get(label)]
        byid=defaultdict(list)
        for r in keep:
            byid[r["Study"]].append(r)
        b_result[label]={
            "source_original_Keep_effect_rows":len(keep),
            "source_original_Study_ID_clusters":len(byid),
            "source_Field_effect_rows":sum(r["Field.Lab"]=="Field" for r in keep),
            "source_numeric_arms_SD_n_rows":sum(original_arm_summaries_positive(r) for r in keep),
            "source_original_Study_IDs_with_Field_and_numeric_candidate":sum(
                any(r["Field.Lab"]=="Field" and original_arm_summaries_positive(r) for r in rs)
                for rs in byid.values()),
            "not_a_fully_extracted_stigmatic_pollen_or_seed_plant_effect_meta":True
        }
    if strict and (b_result["B_flowering_phenology_NOT_reproduction"]["source_original_Keep_effect_rows"]!=26 or
                   b_result["B_plant_seed_fruit_pollen_reproduction_CANDIDATE"]["source_original_Keep_effect_rows"]!=12):
        raise ValueError("Original B data counts are inconsistent")

    totals={
      "source_original_Keep_StageA_broad_effects":len(a_keep),
      "source_original_Keep_StageA_Study_IDs":len(studies),
      "source_stageA_Field_effect_rows":sum(v["source_Keep_Field_effect_rows"] for v in counts.values()),
      "source_stageA_Field_numerically_arithmetically_extractable_effect_rows":sum(
            v["source_Field_arithmetically_valid_positive_arm_means_SD_n_rows"] for v in counts.values()),
      "source_stageA_positive_original_arm_means_SD_n_rows":sum(
            v["source_arithmetically_valid_positive_arm_means_SD_n_rows"] for v in counts.values()),
    }
    return {
      "schema":SCHEMA,
      "Sanders_original_author_1304row_CSV_SHA256":SOURCE_SHA,
      "source_positive_and_negative_and_author_deletes_retained_in_upstream_screen":True,
      "source_report_not_a_WoS_or_Scopus_PRISMA_hit_count":True,
      "same_original_Study_ID_not_many_experimental_lamps":True,
      "stageA_broad_original_source_sufficiency_totals":totals,
      "stageA_original_23_papers_by_primary_biological_process":counts,
      "stageB_original_authors_candidate_source_numeric_evidence":b_result,
      "confirmed_StageA_aquatic_emergence_and_predator_matched_effect_with_TRUE_independent_lamp_SE_k":None,
      "confirmed_StageB_true_per_plant_seed_lnRR_and_wet_affinity_meta_k":None,
      "confirmed_A_water_x_light_meta_regression_beta":None,
      "confirmed_B_wet_affinity_x_light_meta_regression_beta":None,
      "no_same_study_multiarm_shared_control_covariance_derived":True,
      "original_raw_meta_effect_rows_or_Study_ID_exported":False,
      "study_level_data_are_title_screen_for_process_and_habitat_not_PRISMA_fulltext_inclusion":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--triage",type=Path,default=REGISTRY)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    raw,_=zenodo_original_author_mirror()
    rows,ident=original_source_rows(raw,strict=True)
    m=json.loads(a.triage.read_text())
    d=process_quantitative_sufficiency(rows,ident,m,strict=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
    print("ALAN_20261120_STAGEA_SOURCE_PROCESS_NUMERIC_COUNTER",d["stageA_original_23_papers_by_primary_biological_process"])
    print("ALAN_20261120_STAGEB_SOURCE_OUTCOME_NUMERIC_COUNTER",d["stageB_original_authors_candidate_source_numeric_evidence"])
    print("ALAN_20261120_TRUE_META_MODERATOR_STILL_UNESTIMATED",d["confirmed_A_water_x_light_meta_regression_beta"],d["confirmed_B_wet_affinity_x_light_meta_regression_beta"])


if __name__=="__main__":
    main()
