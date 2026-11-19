#!/usr/bin/env python3
"""Original 2023 ALAN biotic review: exact source-coded effect & DOI study audit.

Counts 81 published investigations in 79 papers, *not* random independent
field lighting treatments or comparable effect sizes. Keeps the footer total
separate from scientific records; uses original 0/1 flags without imputation.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import re
from pathlib import Path

from scripts.analysis.alan_seymoure2023_true_workbook_coding_audit_20261119 import (
    SOURCE_SHA256, SOURCE_URL, SOURCE_FILE_ID, SOURCE_NAME, get, source_rows
)

FLAGS = {
    "decrease": "interaction strength (IS) reduction",
    "increase": "IS enhancement",
    "temporal": "temporal mismatch induced",
    "spatial": "spatial mismatch induced",
}
COLUMN_STRENGTH = "strength sum"
COLUMN_TOTAL = "effect sum"


def doi_norm(x):
    x = (x or "").lower().strip()
    x = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", x)
    x = re.sub(r"^doi:\s*", "", x)
    return x.strip().rstrip(" .,;")


def study_key(row):
    doi = doi_norm(row.get("DOI") or row.get("DOI Link"))
    if doi:
        return "doi:"+doi
    authors=row.get("Authors","").casefold().strip()
    title=row.get("Article Title","").casefold().strip()
    year=row.get("Publication Year","").strip()
    if not (authors and title and year):
        raise ValueError("No auditable DOI or stable author-title-year publication key")
    return "fallback:"+repr((authors,title,year))


def _source_scientific_rows(rows, strict=True):
    footers=[row for row in rows if not row.get("interaction_type")]
    if len(footers)!=1:
        raise ValueError(f"Expected one non-biological source footer, got {len(footers)}")
    f=footers[0]
    if not (f.get("DOI","")=="" and f.get("Article Title","")=="" and
            f.get("interaction strength (IS) reduction")=="29" and
            f.get("IS enhancement")=="28" and f.get("strength sum")=="57" and
            f.get("temporal mismatch induced")=="7" and
            f.get("spatial mismatch induced")=="27" and
            f.get("effect sum")==""):
        raise ValueError("The untyped row is not authenticated as source aggregate footer")
    if f["_original_excel_row"]!=max(row["_original_excel_row"] for row in rows):
        raise ValueError("Unspecified source footer not last row")
    scientific=[row for row in rows if row is not f]
    if strict and len(scientific)!=81:
        raise ValueError(f"Expected 81 investigations, got {len(scientific)}")
    return scientific, f


def summarize_source(rows, strict=True):
    scientific, footer=_source_scientific_rows(rows,strict=strict)
    parsed=[]
    for row in scientific:
        if row.get("interaction_type","")=="":
            raise ValueError("Untyped scientific interaction")
        coded={}
        for label,field in FLAGS.items():
            val=row.get(field,"")
            if val not in ("0","1"):
                raise ValueError(f"Original yes/no code missing or changed in {field}")
            coded[label]=int(val)
        if row.get(COLUMN_STRENGTH)!="{}".format(coded["decrease"]+coded["increase"]):
            raise ValueError("Original source strength sum is not equal to original flags")
        num=sum(coded.values())
        if row.get(COLUMN_TOTAL)!="{}".format(num):
            raise ValueError("Original source effect sum is not equal to four original flags")
        coded.update({
            "type":row["interaction_type"],
            "method_original":row.get("method",""),
            "consumer":row.get("Consumer_taxon",""),
            "resource":row.get("Resource_taxon",""),
            "paper":study_key(row),
            "any_change":int(num>0),
            "strength_changed":int(coded["decrease"]+coded["increase"]>0),
            "both_decrease_and_increase":int(coded["decrease"] and coded["increase"]),
            "spatiotemporal_coded":int(coded["spatial"] or coded["temporal"]),
            "year_original":row["Publication Year"],
        })
        parsed.append(coded)
    papers=defaultdict(list)
    for row in parsed:
        papers[row["paper"]].append(row)
    if strict and len(papers)!=79:
        raise ValueError(f"Exact DOI/title publication identity does not produce 79: {len(papers)}")
    type_counts=Counter(row["type"] for row in parsed)
    types={}
    for typ,n in sorted(type_counts.items()):
        rr=[row for row in parsed if row["type"]==typ]
        types[typ]={
            "source_interactions":n,
            "source_publications":len(set(r["paper"] for r in rr)),
            "source_coded_reduction_interactions":sum(r["decrease"] for r in rr),
            "source_coded_enhancement_interactions":sum(r["increase"] for r in rr),
            "source_coded_both_directions_interactions":sum(r["both_decrease_and_increase"] for r in rr),
            "source_coded_any_strength_interactions":sum(r["strength_changed"] for r in rr),
            "source_coded_temporal_mismatch_interactions":sum(r["temporal"] for r in rr),
            "source_coded_spatial_mismatch_interactions":sum(r["spatial"] for r in rr),
            "source_coded_any_change_interactions":sum(r["any_change"] for r in rr),
            "original_method_missing_interactions":sum(not r["method_original"] for r in rr),
        }
    interaction_counts={
        "investigations":len(parsed),
        "coded_reduction":sum(r["decrease"] for r in parsed),
        "coded_enhancement":sum(r["increase"] for r in parsed),
        "both_reduction_and_enhancement":sum(r["both_decrease_and_increase"] for r in parsed),
        "reduction_only":sum(bool(r["decrease"] and not r["increase"]) for r in parsed),
        "enhancement_only":sum(bool(r["increase"] and not r["decrease"]) for r in parsed),
        "no_reduction_or_enhancement_coded":sum(not r["strength_changed"] for r in parsed),
        "any_strength_coded":sum(r["strength_changed"] for r in parsed),
        "coded_temporal":sum(r["temporal"] for r in parsed),
        "coded_spatial":sum(r["spatial"] for r in parsed),
        "coded_temporal_and_spatial":sum(bool(r["temporal"] and r["spatial"]) for r in parsed),
        "coded_time_or_space":sum(r["spatiotemporal_coded"] for r in parsed),
        "coded_strength_and_time_or_space":sum(bool(r["strength_changed"] and r["spatiotemporal_coded"]) for r in parsed),
        "coded_time_or_space_without_strength":sum(bool(r["spatiotemporal_coded"] and not r["strength_changed"]) for r in parsed),
        "coded_any_of_four":sum(r["any_change"] for r in parsed),
        "no_four_flags_recorded":sum(not r["any_change"] for r in parsed),
        "total_flagged_source_effect_entries_not_independent":sum(sum(r[f] for f in FLAGS) for r in parsed),
        "original_method_empty":sum(not r["method_original"] for r in parsed),
    }
    study_counts={
        "publications":len(papers),
        "papers_with_multiple_interaction_rows":sum(len(rr)>1 for rr in papers.values()),
        "source_rows_by_article_multiplicity":dict(sorted(Counter(len(rr) for rr in papers.values()).items())),
        "coded_reduction":sum(any(r["decrease"] for r in rr) for rr in papers.values()),
        "coded_enhancement":sum(any(r["increase"] for r in rr) for rr in papers.values()),
        "coded_both_directions_in_same_article":sum(
            any(r["decrease"] for r in rr) and any(r["increase"] for r in rr)
            for rr in papers.values()),
        "coded_any_strength":sum(any(r["strength_changed"] for r in rr) for rr in papers.values()),
        "coded_any_temporal":sum(any(r["temporal"] for r in rr) for rr in papers.values()),
        "coded_any_spatial":sum(any(r["spatial"] for r in rr) for rr in papers.values()),
        "coded_any_time_or_space":sum(any(r["spatiotemporal_coded"] for r in rr) for rr in papers.values()),
        "coded_any_of_four":sum(any(r["any_change"] for r in rr) for rr in papers.values()),
        "source_studies_all_methods_blank":sum(all(not r["method_original"] for r in rr) for rr in papers.values()),
    }
    out={
        "schema":"alan_seymoure2023_original_81interactions_79publications_codebook_stage_audit_v1",
        "original_source_SHA256":SOURCE_SHA256,
        "verified_footer_original_row":footer["_original_excel_row"],
        "verified_original_footer_29minus_28plus_7temporal_27spatial":True,
        "source_type_direction_and_mismatch_counts":types,
        "interaction_level_original_coded_counts":interaction_counts,
        "publication_level_DOI_aggregated_counts":study_counts,
        "source_method_original_value_counts":dict(sorted(Counter(r["method_original"] or "<BLANK>" for r in parsed).items())),
        "source_consumer_taxon_original_value_counts":dict(sorted(Counter(r["consumer"] or "<BLANK>" for r in parsed).items())),
        "source_resource_taxon_original_value_counts":dict(sorted(Counter(r["resource"] or "<BLANK>" for r in parsed).items())),
        "source_author_text_comparison":{
            "paper_text_predator_prey":56,
            "source_original_predator_prey":type_counts.get("predator-prey",0),
            "paper_text_herbivory":11,
            "source_original_herbivory":type_counts.get("herbivory",0),
            "paper_text_29_enhanced_28_reduced":"Original paper prose states 29 enhancement, 28 reduction; original source has 29 reduction and 28 enhancement",
            "published_75_of_79_altered_interactions":"Cannot derive this original abstract count from only four source-coded effect flags; different outcome definition or discrepancy",
            "paper_text_54_strength_studies":"Compare to DOI-aggregated source coded_any_strength, not 81 investigation rows",
        },
        "source_limits":{
            "source_review_was_not_full_systematic_literature_search":True,
            "paper_is_not_random_light_experiment":True,
            "original_flag_0_means_not_coded_positive_not_verified_ecological_absence":True,
            "source_has_no_effect_size_SE_for_random_effect_meta":True,
            "source_method_missingness_limits_design_moderator":True,
            "direction_and_mismatch_not_independent_effects":True,
            "ALAN_pollen_seed_causal_mediation_not_estimated":True,
            "no_original_study_rows_exported":True,
        },
    }
    # Checks from the source author footer (not journal prose).
    if strict:
        ic=out["interaction_level_original_coded_counts"]
        for name,expected in [("coded_reduction",29),("coded_enhancement",28),
                               ("coded_temporal",7),("coded_spatial",27)]:
            if ic[name]!=expected:
                raise ValueError(f"Source footer count drift: {name}: {ic[name]} != {expected}")
        if ic["total_flagged_source_effect_entries_not_independent"]!=91:
            raise ValueError("Source total flags inconsistent with original footer")
    return out


def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--out",type=Path,required=True)
    args=arg.parse_args()
    meta=json.loads(get(SOURCE_URL))
    files=[f for f in meta.get("files",[]) if f.get("id")==SOURCE_FILE_ID and f.get("name")==SOURCE_NAME]
    if len(files)!=1:
        raise ValueError("Source Figshare file metadata changed")
    original=get(files[0]["download_url"])
    if hashlib.sha256(original).hexdigest()!=SOURCE_SHA256:
        raise ValueError("Original Figshare source SHA drift")
    hdr,rows=source_rows(original)
    if len(hdr)!=20 or len(rows)!=82:
        raise ValueError("Expected original 20 columns and 82 nonempty rows with footer")
    report=summarize_source(rows,strict=True)
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("SEYMOURE_TRUE_79PAPERS",report["publication_level_DOI_aggregated_counts"])
    print("SEYMOURE_TRUE_81INTERACTIONS",report["interaction_level_original_coded_counts"])
    print("SEYMOURE_TRUE_DIRECTION_BY_INTERACTION_TYPE",report["source_type_direction_and_mismatch_counts"])
    print("SEYMOURE_TRUE_METHODS",report["source_method_original_value_counts"])


if __name__=="__main__":
    main()
