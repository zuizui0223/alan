#!/usr/bin/env python3
"""Sanders et al. (2020/2021): source-authentic stage A/B candidate audit.

The original published ALAN meta archive contains 1304 Hedges-d effects in
126 source Study IDs. It has original dark/light means, SD and n, but those n
must NOT automatically be treated as independent physical ALAN installations.
Here we merely classify candidate endpoints for full-text screening and count
which original records have arithmetically usable positive arm summary values.

S1 broad: arthropod population/predation responses, not automatically aquatic
emergence or actual predator species. S2 strict: plant reproductive endpoints,
with separate indirect "insect pollination-related" candidates. These are
LITERATURE SCREENING FLAGS, not fully eligible meta effects. Source publication
location codes don't establish true riparian/terrestrial habitats or a
plant's Ellenberg habitat moisture affinity.

Source matched against author original Zenodo 3998338 and original public
MD5 checksums: 10.5061/dryad.wpzgmsbjn. No original author source rows,
rare study IDs, treatment count values or coordinates are exported to git.
"""
import argparse
import csv
from collections import Counter,defaultdict
import io
import json
import math
import re
import zipfile
from pathlib import Path

from scripts.acquisition.alan_sanders2020_source_126paper_Dryad_preflight_20261119 import (
    zenodo_original_author_mirror,parse_zip,PUBLISHED_PAPERS,PUBLISHED_EFFECTS
)

REQUIRED={
    "Study","ID","Phylum","Class","Order","Species","Measure.description",
    "Category ","Sub.Category","Field.Lab","Control_N","Experimental_N",
    "Control_Mean","Control_SD","Experimental_Mean","Experimental_SD",
    "Night_Lux","Location","DOI","Delete",
}
ARTHROPOD_WORDS=("arthropoda","insecta","arachnida","malacostraca")
PLANT_PHYLUM_WORDS=("tracheophyta","magnoliophyta","streptophyta","spermatophyta","plantae","angiosperms","anthophyta")
PLANT_CLASS_WORDS=("magnoliopsida","liliopsida","rosopsida","eudicotyledones","monocots","dicotyledons","eudicots","equisetopsida","monocotyledoneae")
S1_MARKS=re.compile(r"abundan|predat|populat|captur|arthropod|insect|spider|emerg|prey|pitfall",re.I)
S2_REPRO_MARKS=re.compile(r"pollinat|pollen|fruit|seed|ovule|germin|reproduc|stigma",re.I)
S2_VISIT_MARKS=re.compile(r"flower.{0,20}visit|visitation|nectar|foraging.{0,30}flower|pollinat",re.I)
S2_PHENO_MARKS=re.compile(r"flower|phenolog|bloom|anthesis|flowering|bud.?burst",re.I)


def valnum(s):
    try:
        f=float(str(s).strip().replace(",",""))
        return f if math.isfinite(f) else None
    except (ValueError,TypeError):
        return None


def original_source_rows(blob,strict=True):
    base=parse_zip(blob,strict=strict)
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        original=z.read("Sanders_et_al_datafile.csv")
    text=original.decode("utf-8-sig",errors="replace")
    originalreader=csv.DictReader(io.StringIO(text,newline=""))
    fields=originalreader.fieldnames
    if not fields or not REQUIRED.issubset(set(fields)):
        raise ValueError("Sanders original summary-arm and classification columns not source verified")
    rows=list(originalreader)
    if strict and (len(rows)!=PUBLISHED_EFFECTS or
                   len(set(r["Study"] for r in rows if r["Study"]))!=PUBLISHED_PAPERS):
        raise ValueError("Published original 126 study / 1304 source effects mismatch")
    return rows,base


def source_panel_census(rows,source_identity,strict=True):
    s1=[];s2direct=[];s2visits=[]
    category_counts=defaultdict(Counter)
    stage1_class=Counter()
    labelled=Counter()
    candidate_numeric=Counter()
    study_ids=defaultdict(set)
    for row in rows:
        category=(row.get("Category ") or "").strip()
        sub=(row.get("Sub.Category") or "").strip()
        phylum=(row.get("Phylum") or "").strip().lower()
        clazz=(row.get("Class") or "").strip().lower()
        order=(row.get("Order") or "").strip().lower()
        desc=(row.get("Measure.description") or "").strip()
        species=(row.get("Species") or "").strip()
        key=row.get("Study","").strip()
        # Do not turn the original broad "populations" category into a
        # predator abundance signal. Only original endpoint subcategory
        # and original measured variable description can carry that label.
        class_text=" ".join((sub,desc))
        category_counts["original_Category"][category or "<BLANK>"]+=1
        category_counts["original_SubCategory"][sub or "<BLANK>"]+=1
        category_counts["original_Phylum"][phylum or "<BLANK>"]+=1
        category_counts["original_Class"][clazz or "<BLANK>"]+=1
        category_counts["original_Field_Lab"][(row.get("Field.Lab") or "").strip() or "<BLANK>"]+=1
        category_counts["original_Delete"][(row.get("Delete") or "").strip() or "<BLANK>"]+=1
        category_counts["original_Activity"][(row.get("Activity") or "").strip() or "<BLANK>"]+=1
        category_counts["original_Noct"][(row.get("Noct") or "").strip() or "<BLANK>"]+=1
        if not key and strict:
            raise ValueError("Original Sanders Study ID required for dependence")
        plant=(any(x in phylum for x in PLANT_PHYLUM_WORDS) or
               any(x in clazz for x in PLANT_CLASS_WORDS))
        arthropod=(any(x in phylum for x in ARTHROPOD_WORDS) or
                   any(x in clazz for x in ARTHROPOD_WORDS) or
                   "araneae" in order or "lepidoptera" in order or "diptera" in order)
        flags={
            "S1_arthropod_predation_abundance_candidate":arthropod and bool(S1_MARKS.search(class_text)),
            "S2_plant_reproductive_endpoint_candidate":plant and bool(S2_REPRO_MARKS.search(class_text)),
            "S2_plant_flowering_phenology_candidate_NOT_seed_outcome":plant and bool(S2_PHENO_MARKS.search(class_text)),
            "S2_invertebrate_potential_flower_visitation_candidate":arthropod and bool(S2_VISIT_MARKS.search(class_text)),
        }
        nL=valnum(row.get("Experimental_N"))
        nD=valnum(row.get("Control_N"))
        mL=valnum(row.get("Experimental_Mean"))
        mD=valnum(row.get("Control_Mean"))
        sdL=valnum(row.get("Experimental_SD"))
        sdD=valnum(row.get("Control_SD"))
        numeric=(all(isinstance(x,float) for x in (nL,nD,mL,mD,sdL,sdD))
                 and nL>=2 and nD>=2 and mL>0 and mD>0 and sdL>=0 and sdD>=0)
        for name,yes in flags.items():
            if not yes:continue
            labelled[name]+=1
            if numeric:
                candidate_numeric[name]+=1
            study_ids[name].add(key)
            if name.startswith("S1"):
                stage1_class[(phylum or "<BLANK>",clazz or "<BLANK>")]+=1
    broad={
      label:{
        "source_original_effect_records_preliminary_keyword_and_taxonomy":n,
        "source_unique_original_Study_ids_preliminary":len(study_ids[label]),
        "source_positive_arms_SD_n_arithmetically_available_records":candidate_numeric[label],
        "original_physical_illumination_clusters_and_shared_control_covariance_NOT_validated":True,
        "not_a_verified_direct_aquatic_or_plant_pollination_lnRR_meta_input":True,
      }
      for label,n in sorted(labelled.items())
    }
    output={
      "schema":"alan_original_Sanders2020_126Study_1304effects_stageAB_exploratory_census_v1",
      "source_dataset_DOI":"10.5061/dryad.wpzgmsbjn",
      "original_author_CSV_SHA256":source_identity["original_CSV_SHA256"],
      "original_author_R_script_SHA256":source_identity["original_R_code_SHA256"],
      "source_meta_original_hedges_d_not_lnRR":True,
      "source_original_Study_count":len(set(r["Study"] for r in rows if r.get("Study"))),
      "source_original_Study_effect_records":len(rows),
      "literal_original_source_category_counts":{
         field:dict(sorted(v.items(),key=lambda x:(-x[1],x[0])))
         for field,v in category_counts.items()
      },
      "stageAB_broad_record_screening_not_eligible_effects":broad,
      "original_published_Hedges_d_source_effects":1304,
      "original_meta_126papers_not_126_unique_physical_light_programmes_proven":True,
      "source_original_effect_n_values_not_confirmed_independent_lamp_sites":True,
      "requires_publication_fulltext_taxonomy_and_habitat_before_stageA":True,
      "requires_plant_conspecific_pollen_and_true_seed_denominator_before_stageB":True,
      "no_aquatic_vs_terrestrial_moderator_fitted":True,
      "no_wet_plant_moisture_moderator_fitted":True,
      "source_raw_study_IDs_or_rows_exported":False,
      "source_treatment_means_SD_not_exported":True,
    }
    if strict and (output["source_original_Study_count"]!=126 or
                   output["source_original_Study_effect_records"]!=1304):
        raise ValueError("Source original dataset shape changed")
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z,metadata=zenodo_original_author_mirror()
    rows,original=original_source_rows(z,strict=True)
    result=source_panel_census(rows,original,strict=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("SANDERS126_AUTHENTIC_SUBCATEGORY_FREQUENCIES",result["literal_original_source_category_counts"]["original_SubCategory"])
    print("SANDERS126_AUTHENTIC_PHYLUM_COUNTS",result["literal_original_source_category_counts"]["original_Phylum"])
    print("SANDERS126_AUTHENTIC_CLASS_COUNTS",result["literal_original_source_category_counts"]["original_Class"])
    print("SANDERS126_STAGE1_STAGE2_INITIAL_SOURCE_RECORDS",result["stageAB_broad_record_screening_not_eligible_effects"])
    print("SANDERS126_STAGE_META_AUTHOR_ORIGINAL_SOURCE_GOOD_PRISMA_UNSEARCHED",result["source_original_Study_count"])


if __name__=="__main__":
    main()
