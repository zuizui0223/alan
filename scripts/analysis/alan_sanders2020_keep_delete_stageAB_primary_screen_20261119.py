#!/usr/bin/env python3
"""Author-original ALAN meta source 1304 rows: Keep/Delete × field/lab × stage.

Sanders et al. Nat Ecol Evol 2021 original 2020 public data,
Zenodo 3998338 independently file-MD5 verified original source:
data CSV SHA d11d1e29c05826f2b4a5c1c2d63b9c195dc2effbedd1ff531afa7b6bdf497693.
Original rows 1304, Sanders' 'Keep' field is a SOURCE exclusion flag, NOT
automatically the inclusion rule for this ALAN biotic interaction review.
This table is a literature source-triage, NOT meta-regression, ecological
variance or proof of 26 independent lighting installations.

Returns source paper bibliographic DOI + Title only for FULLTEXT screening,
never original effect estimate, row identifiers, plant/animal observations,
coordinates, experimental treatment arm mean/SD/n or raw source lines.
"""
import argparse
from collections import Counter,defaultdict
import json
import re
from pathlib import Path

from scripts.analysis.alan_sanders2020_stageAB_source_census_20261119 import (
    original_source_rows,source_panel_census,valnum,
    PLANT_PHYLUM_WORDS,PLANT_CLASS_WORDS,ARTHROPOD_WORDS,
    S1_MARKS,S2_REPRO_MARKS,S2_PHENO_MARKS,S2_VISIT_MARKS
)
from scripts.acquisition.alan_sanders2020_source_126paper_Dryad_preflight_20261119 import (
    zenodo_original_author_mirror,
)

SOURCE_SHA256="d11d1e29c05826f2b4a5c1c2d63b9c195dc2effbedd1ff531afa7b6bdf497693"
SCHEMA="alan_original_sanders_1304_Keep_Delete_stageAB_experimental_unit_gate_v1"

def is_plant(r):
    phy=str(r.get("Phylum") or "").strip().lower()
    cla=str(r.get("Class") or "").strip().lower()
    return any(k in phy for k in PLANT_PHYLUM_WORDS) or any(k in cla for k in PLANT_CLASS_WORDS)

def is_arthropod(r):
    phy=str(r.get("Phylum") or "").strip().lower()
    cla=str(r.get("Class") or "").strip().lower()
    order=str(r.get("Order") or "").strip().lower()
    return (any(k in phy for k in ARTHROPOD_WORDS)
           or any(k in cla for k in ARTHROPOD_WORDS)
           or any(k in order for k in ("araneae","lepidoptera","diptera")))

def source_label_stage(r):
    text=" ".join((str(r.get("Sub.Category") or ""),str(r.get("Measure.description") or "")))
    plant=is_plant(r)
    arthropod=is_arthropod(r)
    return {
      "A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN":arthropod and bool(S1_MARKS.search(text)),
      "B_flowering_phenology_NOT_reproduction":plant and bool(S2_PHENO_MARKS.search(text)),
      "B_plant_seed_fruit_pollen_reproduction_CANDIDATE":plant and bool(S2_REPRO_MARKS.search(text)),
      "B_invertebrate_possible_flower_visits_not_pollen":arthropod and bool(S2_VISIT_MARKS.search(text)),
    }

def original_arm_summaries_positive(r):
    vals=[valnum(r.get(k)) for k in (
        "Control_N","Experimental_N","Control_Mean","Experimental_Mean",
        "Control_SD","Experimental_SD"
    )]
    return (all(isinstance(x,float) for x in vals)
            and vals[0]>=2 and vals[1]>=2 and vals[2]>0 and vals[3]>0
            and vals[4]>=0 and vals[5]>=0)

def doi_norm(raw):
    text=str(raw or "").strip().lower()
    text=re.sub(r"^https?://(?:dx\.)?doi\.org/","",text)
    text=re.sub(r"^doi:\s*","",text)
    text=text.strip().rstrip(" .;,")
    return text if re.fullmatch(r"10\.\d{4,9}/\S+",text) else None

def audit_keep_delete(rows,base,strict=True):
    if base.get("original_CSV_SHA256")!=SOURCE_SHA256:
        raise ValueError("Author original Sanders CSV SHA mismatch")
    if strict and len(rows)!=1304:
        raise ValueError("Published original 1304 source row size changed")
    classes=Counter((r.get("Delete") or "").strip() or "<BLANK>" for r in rows)
    if strict and classes!={"Keep":1109,"Delete":195}:
        raise ValueError("Published Sanders original Keep/Delete provenance drift")
    fields=Counter((r.get("Field.Lab") or "").strip() or "<BLANK>" for r in rows)
    if strict and fields!={"Field":770,"Laboratory":534}:
        raise ValueError("Source Field.Lab 770/534 has changed")
    groups=defaultdict(list)
    for r in rows:
        if not r.get("Study"):
            raise ValueError("Original Study identity not available")
        for label,match in source_label_stage(r).items():
            if match:groups[label].append(r)
    stage_summary={}
    bibliography={}
    for stage,values in sorted(groups.items()):
        kept=[x for x in values if x["Delete"]=="Keep"]
        deleted=[x for x in values if x["Delete"]=="Delete"]
        if len(kept)+len(deleted)!=len(values):
            raise ValueError("Original source status unexpectedly unknown in stage")
        valid=[x for x in kept if original_arm_summaries_positive(x)]
        fieldkept=[x for x in kept if x.get("Field.Lab")=="Field"]
        labkept=[x for x in kept if x.get("Field.Lab")=="Laboratory"]
        identifiers=lambda vals:len(set(x.get("Study") for x in vals))
        stage_summary[stage]={
            "all_source_effect_rows_pre_source_exclusion":len(values),
            "all_source_unique_Study_IDs_pre_source_exclusion":identifiers(values),
            "author_original_Keep_effect_rows":len(kept),
            "author_original_Delete_effect_rows":len(deleted),
            "author_original_Keep_unique_Study_IDs":identifiers(kept),
            "author_original_Delete_unique_Study_IDs":identifiers(deleted),
            "author_original_Keep_Field_effect_rows":len(fieldkept),
            "author_original_Keep_Laboratory_effect_rows":len(labkept),
            "author_original_Keep_Field_Study_IDs":identifiers(fieldkept),
            "author_original_Keep_Laboratory_Study_IDs":identifiers(labkept),
            "author_original_Keep_positive_arm_means_SD_n_candidate_rows":len(valid),
            "Keep_status_is_not_PRISMA_inclusion":True,
            "source_Study_ID_is_not_independent_light_installation":True,
        }
        pubs=defaultdict(list)
        for r in kept:
            pubs[str(r["Study"])].append(r)
        bibliography[stage]=[]
        for paper_rows in pubs.values():
            dois=sorted({p for r in paper_rows for p in [doi_norm(r.get("DOI"))] if p})
            titles=sorted({str(r.get("Title") or "").strip() for r in paper_rows if str(r.get("Title") or "").strip()})
            if len(dois)>1:
                # The paper might cite more than one original research study.
                # Use uncertainty status rather than merging by DOI.
                identity="SOURCE_DOI_CONFLICT"
            else:
                identity="SINGLE_OR_MISSING_SOURCE_DOI"
            rec={
                "source_published_DOIs":dois,
                "source_paper_title":titles[0] if len(titles)==1 else None,
                "title_conflict_or_missing":len(titles)!=1,
                "doi_conflict":len(dois)>1,
                "source_records_in_Keep_status":len(paper_rows),
                "author_2019_source_Decision":"Keep",
                "fulltext_riparian_vs_terrestrial":"UNKNOWN_REQUIRES_PRIMARY_STUDY",
                "independently_randomised_light_installations":"UNKNOWN_REQUIRES_PRIMARY_METHODS",
                "original_effect_standard_error_of_lnRR_from_independent_plot_units":"NOT_VERIFIED",
                "not_an_included_PRISMA_study_yet":True,
            }
            bibliography[stage].append(rec)
        bibliography[stage].sort(key=lambda z:(z["source_published_DOIs"],z["source_paper_title"] or ""))
    result={
        "schema":SCHEMA,
        "original_Sanders_author_CSV_SHA256":SOURCE_SHA256,
        "original_Sanders_review_DOI":"10.1038/s41559-020-01322-x",
        "published_1304_original_rows":len(rows),
        "author_original_keep_delete":dict(classes),
        "author_original_Field_Lab":dict(fields),
        "stage_AB_source_screening_with_original_keep_delete_and_design":stage_summary,
        "candidate_publication_bibliography_for_fulltext_screening":bibliography,
        "source_original_Study_IDs_hidden":True,
        "original_source_effect_rows_not_redistributed":True,
        "original_field_Study_not_equal_physical_light_experiment":True,
        "no_verified_aquatic_vs_terrestrial_meta_moderator":True,
        "no_verified_plant_habitat_moisture_meta_moderator":True,
        "no_random_effect_meta_fitted":True,
        "original_author_Delete_rows_not_assumed_biologically_null":True,
    }
    if strict:
        s1=stage_summary.get("A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN",{})
        if s1.get("all_source_effect_rows_pre_source_exclusion")!=170:
            raise ValueError("Original arthropod broad source screening changed")
        s2=stage_summary.get("B_plant_seed_fruit_pollen_reproduction_CANDIDATE",{})
        if s2.get("all_source_effect_rows_pre_source_exclusion")!=12:
            raise ValueError("Original plant reproductive source screening changed")
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    source,metadata=zenodo_original_author_mirror()
    original,base=original_source_rows(source,strict=True)
    result=audit_keep_delete(original,base,strict=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("SANDERS_ORIGINAL_126PAPER_KEEP_DELETE",result["author_original_keep_delete"])
    print("SANDERS_STAGE_AB_TRUE_KEEP_DELETE_NUMBERS",result["stage_AB_source_screening_with_original_keep_delete_and_design"])
    print("SANDERS_STAGEAB_CITATION_COUNTS_WITH_TRUE_DOI",
          {k:{"paper_keys":len(v),"papers_with_doi":sum(bool(w["source_published_DOIs"]) for w in v)}
           for k,v in result["candidate_publication_bibliography_for_fulltext_screening"].items()})
    print("SANDERS_STAGE_AB_NO_GLOBAL_LIGHT_META",result["no_random_effect_meta_fitted"])

if __name__=="__main__":
    main()
