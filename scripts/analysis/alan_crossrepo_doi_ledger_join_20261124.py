#!/usr/bin/env python3
"""Join frozen 62-repo DOI-string inventory to ALAN research admission ledgers.

Snapshot only; DOI strings are NOT validated literature or independent studies.
Read ONLY source DOI identifiers, repo paths, and stored ALAN evidence tables.
Do not turn unexplained identifiers into meta effects or PRISMA inclusions.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict, Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPECTED_BASELINE = {
    "repositories_reported": 62,
    "repositories_with_full_tree": 62,
    "priority_files_actually_read": 539,
    "distinct_literal_DOIs_in_read_sources": 676,
}
CONTEXT = re.compile(
    r"alan|night.?light|light.?pollut|nocturn|artificial.?light|"
    r"pollinat|arthropod|insect|predat|herbivor|"
    r"biotic|trophic|seed|fruit|flower|reproduct|"
    r"experimental.?light|ecolog", re.I
)
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$", re.I)

def csv_rows(path: Path) -> list[dict[str,str]]:
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if rows and any(None in row for row in rows):
        raise ValueError(f"STOP: CSV row has more columns than header: {path}")
    return rows

def normalize_doi(x: str) -> str:
    """DOI-string join normalization, never scientific identity validation."""
    z = str(x).strip().lower()
    for prefix in ("https://dx.doi.org/", "https://doi.org/", "http://doi.org/", "doi:"):
        if z.startswith(prefix):
            z = z[len(prefix):]
    return z.strip().rstrip(".,;: \t")

def register(known: dict[str,set[str]], doi: str, label: str) -> None:
    key=normalize_doi(doi)
    if DOI_RE.fullmatch(key):
        known[key].add(label)

def load_existing(root: Path=ROOT) -> tuple[dict[str,set[str]], dict[str,dict]]:
    data = root / "database/alan_experiments"
    known: dict[str,set[str]]=defaultdict(set)
    metadata: dict[str,dict]={}

    def ingest_csv(path, field, label, skip_empty=True):
        rows=csv_rows(data/path)
        for row in rows:
            raw=row.get(field, "")
            if not raw:
                if skip_empty:continue
                raise ValueError(f"STOP: empty DOI {path}")
            register(known,raw,label)
        return len(rows)

    ingest_csv(Path("studies.csv"),"doi","ALAN_MAIN_STUDY_REGISTER")
    ingest_csv(Path("arthropod_expansion_20261122/records.csv"),"doi","ALAN_EXTRA_341_SEARCH")
    ingest_csv(Path("reproductive_meta_admission_candidates_20261124.csv"),"doi","ALAN_REPRODUCTIVE_ADMISSION")
    ingest_csv(Path("cross_repo_fulltext_design_queue_20261124.csv"),"doi","ALAN_PENDING_55")
    ingest_csv(Path("cross_repo_iwe_alan_source_crosswalk_20261124.csv"),"doi","IWE_9_CROSSWALK")

    reviews=csv_rows(data/"search/unified_review_inventory.csv")
    if len(reviews)!=553:
        raise ValueError(f"STOP: unified review inventory changed: {len(reviews)}")
    for row in reviews:
        key = row["report_key"]
        family=row["publication_family_key"]
        for token in (key,family):
            register(known,token,"ALAN_BASE_553_REVIEW_KEY")
        id_doi=normalize_doi(key)
        if DOI_RE.fullmatch(id_doi):
            x=metadata.setdefault(id_doi,{"primary_read_documented":False,
                                         "pending_documented_fulltext_followup":False,
                                         "coverage_status":set()})
            x["primary_read_documented"]|=row["primary_read_documented"].lower()=="true"
            x["pending_documented_fulltext_followup"]|=row["pending_documented_fulltext_followup"].lower()=="true"
            x["coverage_status"].add(row["coverage_status"])

    return known,metadata

def make_crosswalk(inventory: Path, root: Path=ROOT):
    summary=json.loads((inventory/"summary.json").read_text(encoding="utf-8"))
    for key,expected in EXPECTED_BASELINE.items():
        if summary.get(key)!=expected:
            raise ValueError(f"STOP: crossrepo frozen inventory identity drift {key}: {summary.get(key)} vs {expected}")
    if summary.get("all_references_full_text_adjudicated") is not False:
        raise ValueError("STOP: inventory cannot claim completed full-text review")
    if summary.get("all_default_branch_text_files_scanned") is not False:
        raise ValueError("STOP: expected sample-limited DOI inventory")
    occurrences=csv_rows(inventory/"doi_occurrence_sources.csv")
    repo_coverage=csv_rows(inventory/"repository_coverage.csv")
    if len(repo_coverage)!=62:
        raise ValueError("STOP: not 62 repository coverage rows")
    if len(set(z["repository"].lower() for z in repo_coverage))!=62:
        raise ValueError("STOP: duplicate repository coverage")
    if len(occurrences)!=summary["DOI_source_occurrences_distinct_by_file"]:
        raise ValueError("STOP: source DOI occurrence row count changed")
    grouped=defaultdict(list)
    for row in occurrences:
        key=normalize_doi(row["doi"])
        if not DOI_RE.fullmatch(key):
            raise ValueError(f"STOP: syntactically invalid observed DOI: {key}")
        if not row.get("repo") or not row.get("path") or not row.get("git_blob_sha") or not row.get("line"):
            raise ValueError("STOP: occurrence lacks traceable Git source")
        grouped[key].append(row)
    if len(grouped)!=676:
        raise ValueError(f"STOP: expected 676 observed DOI string identifiers, got {len(grouped)}")

    known,metadata=load_existing(root)
    cross=[]
    for doi,items in sorted(grouped.items()):
        repos=sorted({x["repo"] for x in items})
        filepaths=sorted({x["repo"]+"@"+x["git_blob_sha"]+":"+x["path"]+":"+x["line"] for x in items})
        source_flags=sorted(known.get(doi,set()))
        contextual=[
            row for row in items if CONTEXT.search(row["path"]) or
            row["repo"].lower() in ("zuizui0223/alan","zuizui0223/iwe")
        ]
        has_alan="zuizui0223/alan" in (r.lower() for r in repos)
        other_count=sum(r.lower() not in ("zuizui0223/alan","zuizui0223/iwe") for r in repos)
        status=("EXACT_ID_ALREADY_IN_ALAN_OR_IWE_LEDGER" if source_flags
                else "UNMATCHED_DOI_STRING_SOURCE_LEAD_NOT_VERIFIED")
        lane=("ALAN_OR_IWE_CONTEXT" if contextual
              else "ONLY_OTHER_REPOSITORY_CONTEXT_NOT_ALAN_VALIDATED")
        priority=("P1_EXACT_UNMATCHED_WITH_ALAN_IWE_SOURCE" if not source_flags and has_alan
                  else "P2_EXACT_UNMATCHED_IN_RELEVANT_DOC_PATH" if not source_flags and contextual
                  else "P3_EXACT_UNMATCHED_OTHER_CONTEXT" if not source_flags
                  else "KNOWN_REFERENCE_NOT_NEW_STUDY")
        review=metadata.get(doi,{})
        cross.append({
            "doi_string":doi,
            "exact_ledger_status":status,
            "matched_source_ledgers":";".join(source_flags),
            "triage_priority":priority,
            "source_context_gate":lane,
            "repo_count":len(repos),
            "source_files":len(filepaths),
            "source_repositories":";".join(repos),
            "source_file_locations":";".join(filepaths[:5]),
            "remaining_source_locations_not_displayed":max(0,len(filepaths)-5),
            "already_primary_read_in_ALAN_review":review.get("primary_read_documented",False),
            "pending_primary_fulltext_in_ALAN_review":review.get("pending_documented_fulltext_followup",False),
            "novel_primary_article_confirmed":False,
            "independent_experiment_programme_confirmed":False,
            "meta_eligible_effect_confirmed":False,
        })

    counts=Counter(row["triage_priority"] for row in cross)
    source_status=Counter(row["exact_ledger_status"] for row in cross)
    detail={
        "snapshot":"2026-11-24",
        "source_run_id":45328217854,
        "source_repository_sampled_not_full_content":True,
        "original_source_inventory_distinct_doi_strings":676,
        "doi_join_rows":len(cross),
        "exact_match_rows":source_status.get("EXACT_ID_ALREADY_IN_ALAN_OR_IWE_LEDGER",0),
        "not_exact_matched_doi_strings":source_status.get("UNMATCHED_DOI_STRING_SOURCE_LEAD_NOT_VERIFIED",0),
        "triage_priority_counts":dict(counts),
        "source_repositories_visited":62,
        "source_files_read":539,
        "all_PDF_and_supplementary_fulltexts_read":False,
        "abstract_title_screen_and_author_data_validate_pending":True,
        "new_meta_programmes_added":0,
        "accepted_effects_added":0,
        "strict_warning":"Exact DOI strings absent from ledgers are NOT automatically relevant or new papers. Joining selected path-only citations cannot adjudicate title, source integrity, physical independence or study endpoints. Full-source screening, missing non-DOI reports, unmerged PRs and WoS/Scopus remain.",
        "source_warning":"The 553 review inventory report_keys may have no DOI; exact matched tally underestimates report overlap. Do not call absence a new experiment.",
    }
    return cross,detail

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--inventory",required=True,type=Path)
    ap.add_argument("--output-dir",default=Path("results/alan_crossrepo_doi_ledger_join_20261124"),type=Path)
    args=ap.parse_args()
    rows,summary=make_crosswalk(args.inventory)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/"doi_crossrepo_ALAN_ledger_join.csv").open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),quoting=csv.QUOTE_MINIMAL)
        w.writeheader();w.writerows(rows)
    for priority,short_name in (
        ("P1_EXACT_UNMATCHED_WITH_ALAN_IWE_SOURCE","p1_alan_iwe_unmatched_doi_source_leads.csv"),
        ("P2_EXACT_UNMATCHED_IN_RELEVANT_DOC_PATH","p2_relevant_path_unmatched_doi_leads.csv"),
    ):
        subset=[r for r in rows if r["triage_priority"]==priority]
        with (args.output_dir/short_name).open("w",newline="",encoding="utf-8") as file:
            w=csv.DictWriter(file,fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(subset)
        if priority.startswith("P1"):
            print("ALAN_CROSSREPO_P1_SOURCE_LEADS",
                  json.dumps([{"doi":r["doi_string"],"repo":r["source_repositories"],
                               "source":r["source_file_locations"].split(";")[0]}
                              for r in subset],ensure_ascii=False),flush=True)
    (args.output_dir/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("ALAN_CROSSREPO_DOI_JOIN",json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
