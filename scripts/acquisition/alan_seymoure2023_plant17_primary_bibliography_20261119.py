#!/usr/bin/env python3
"""Source-pinned public primary bibliography of the 17 ALAN plant-interaction papers.

This is bibliographic provenance, not a copy of the review's source study rows.
No inference about causal lighting effects or published conclusions is made.
"""
import argparse
from collections import Counter,defaultdict
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

from scripts.analysis.alan_seymoure2023_true_workbook_coding_audit_20261119 import (
    SOURCE_SHA256,SOURCE_URL,SOURCE_FILE_ID,SOURCE_NAME,get,source_rows,
)
from scripts.analysis.alan_seymoure2023_original_direction_study_units_20261119 import (
    _source_scientific_rows,study_key,doi_norm,
)

def source_bibliography(rows,strict=True):
    real,_=_source_scientific_rows(rows,strict=strict)
    matched=[x for x in real if x["interaction_type"] in ("pollination","herbivory")]
    papers=defaultdict(list)
    for r in matched:
        papers[study_key(r)].append(r)
    if strict and not(len(matched)==18 and len(papers)==17):
        raise ValueError("Original source plant-study coverage changed")
    result=[]
    for paper in papers.values():
        first=paper[0]
        doi=doi_norm(first.get("DOI") or first.get("DOI Link"))
        if not doi:
            raise ValueError("Primary source paper lacks DOI; cannot silently infer it")
        # Validate all original source rows of this publication share identity.
        if any(doi_norm(x.get("DOI") or x.get("DOI Link")) not in ("",doi) for x in paper):
            raise ValueError("Original paper DOI identity conflict")
        years={r["Publication Year"] for r in paper}
        if len(years)!=1:raise ValueError("Publication year inconsistency")
        titles={r["Article Title"].strip() for r in paper}
        if len(titles)!=1:raise ValueError("Original article title inconsistency")
        result.append({
            "doi":doi,
            "original_article_title":first["Article Title"].strip(),
            "original_publication_year":first["Publication Year"],
            "original_2023_review_source_types":sorted(set(r["interaction_type"] for r in paper)),
            "original_2023_review_method_string_present":any(bool(r.get("method")) for r in paper),
            "original_2023_review_source_interaction_count":len(paper),
            "verified_original_research_light_installation_units":None,
            "verified_primary_endpoint_effect_size":None,
            "verified_primary_endpoint_SE":None,
            "primary_data_extraction_status":"SOURCE_STOP_NEEDS_PRIMARY_ARTICLE",
        })
    result.sort(key=lambda r:(r["original_publication_year"],r["doi"]))
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    meta=json.loads(get(SOURCE_URL))
    fs=[x for x in meta.get("files",[]) if x.get("id")==SOURCE_FILE_ID and x.get("name")==SOURCE_NAME]
    if len(fs)!=1:raise ValueError("Original author uploaded file not identified")
    b=get(fs[0]["download_url"])
    if hashlib.sha256(b).hexdigest()!=SOURCE_SHA256:
        raise ValueError("Original author workbook SHA mismatch")
    headers,rows=source_rows(b)
    if len(headers)!=20:raise ValueError("Unexpected original source schema")
    bibliography=source_bibliography(rows,strict=True)
    out={
        "schema":"alan_seymoure2023_original_source_17plant_paper_bibliography_v1",
        "original_source_sha256":SOURCE_SHA256,
        "source_original_review_DOI":"10.1098/rstb.2022.0356",
        "coverage":"ALL 6 source pollination + 12 herbivory investigations, 17 distinct papers, not direction-selected",
        "paper_count":len(bibliography),
        "original_guild_membership_distribution":dict(sorted(Counter("+".join(x["original_2023_review_source_types"]) for x in bibliography).items())),
        "original_research_lighting_cluster_SEs_verified":0,
        "source_paper_level_public_bibliography":bibliography,
        "bibliographic_titles_and_DOIs_public_not_original_experiment_rows":True,
        "no_full_original_study_rows_or_experimental_data_exported":True,
    }
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_PLANT_17_PUBLIC_BIBLIOGRAPHY",json.dumps([
        {"doi":r["doi"],"year":r["original_publication_year"],
         "title":r["original_article_title"],
         "source_types":r["original_2023_review_source_types"]}
        for r in bibliography
    ],ensure_ascii=False))
    print("ALAN_PLANT17_ORIGINAL_BIBLIOGRAPHY_SOURCE_GATES_PASS")

if __name__=="__main__":
    main()
