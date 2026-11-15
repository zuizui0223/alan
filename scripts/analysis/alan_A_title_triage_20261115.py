#!/usr/bin/env python3
"""Title-stage TRIAGE, never a full-text inclusion decision.

Every title is retained. The triage queues only order human full-text reading;
unknown titles are not forced into an exclusion bin.
"""
from __future__ import annotations

import csv
import json
import pathlib
import re
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE = "openalex_paged_metadata.csv"
OUTPUT = "openalex_title_triage.csv"
FIELDS = (
    "openalex_id", "doi", "year", "title", "query_ids",
    "screen_queue", "document_type_hint", "registered_id",
    "registered_design", "registered_A_stage", "source_work_url",
    "fulltext_reviewed", "effect_variance_verified",
    "independent_dependence_verified", "eligibility_confirmed",
)

REPRO = re.compile(
    r"\b(seed set|seed production|seed number|seed yield|seedling recruitment|"
    r"fruit set|fruit production|fruit infestation|fruiting|fecundity|"
    r"pollination success|pollination and fitness|reproductive output|"
    r"plant fitness|reproductive capacity|plant reproductive)\b", re.I
)
PLANT_POLLIN = re.compile(
    r"\b(plant|plants|vegetation|floral|flower|flowers|flowering|pollinat|"
    r"pollen|moth|hawkmoth|nectar|orchid|silene|yucca|pollinator|"
    r"gecko|lizard|nectarivore|angiosperm|urban forest|phenolog)\w*", re.I
)
DOCUMENT = re.compile(r"^(raw data|supplementary material|dataset)|\bdeposited data\b", re.I)
REVIEW = re.compile(r"\b(review|meta.analysis|bibliometric|perspective|editorial|thesis)\b", re.I)


def queue_for(title: str, registry: dict | None) -> tuple[str, str]:
    title = title or ""
    kind = "unresolved_publication_type"
    if DOCUMENT.search(title):
        kind = "likely_secondary_data_record"
    elif REVIEW.search(title):
        kind = "likely_synthesis_or_thesis"
    else:
        kind = "original_or_unchecked"
    if registry is not None:
        stage = registry.get("a_primary_stage", "")
        if stage.startswith("exclude_A_no_reproduction"):
            return "REGISTERED_MECHANISM_ONLY", kind
        if stage.startswith("exclude_A_not_manipulated"):
            return "REGISTERED_OBSERVATIONAL_B_ONLY", kind
        return "REGISTERED_REPRODUCTIVE_FULLTEXT_PRIORITY", kind
    if DOCUMENT.search(title) or REVIEW.search(title):
        return "SECONDARY_RECORD_CHECK_PARENT", kind
    if REPRO.search(title):
        return "UNREGISTERED_REPRO_MAY_BE_RELEVANT", kind
    if PLANT_POLLIN.search(title):
        return "PLANT_POLLINATOR_MECHANISM_OR_UNSPECIFIED", kind
    return "UNCERTAIN_REQUIRES_TITLE_ABSTRACT_REVIEW", kind


def read_registry(path: pathlib.Path) -> dict:
    registry = json.loads(path.read_text(encoding="utf-8"))
    return {str(p["doi"]).lower(): p for p in registry["entries"]}


def screen_rows(rows: list[dict], known: dict) -> list[dict]:
    out = []
    for row in rows:
        doi = str(row["doi"]).lower()
        entry = known.get(doi)
        queue, kind = queue_for(row["title"], entry)
        out.append({
            "openalex_id": row["openalex_id"], "doi": doi, "year": row["year"],
            "title": row["title"], "query_ids": row["query_ids"],
            "screen_queue": queue, "document_type_hint": kind,
            "registered_id": "" if entry is None else entry["id"],
            "registered_design": "" if entry is None else entry["light_design"],
            "registered_A_stage": "" if entry is None else entry["a_primary_stage"],
            "source_work_url": row["openalex_id"],
            "fulltext_reviewed": "no", "effect_variance_verified": "unknown",
            "independent_dependence_verified": "unknown",
            "eligibility_confirmed": "NOT_ASSESSED",
        })
    return out


def triage(input_dir: pathlib.Path, registry_path: pathlib.Path) -> dict:
    with (input_dir / SOURCE).open(newline="", encoding="utf-8") as f:
        source = list(csv.DictReader(f))
    known = read_registry(registry_path)
    rows = screen_rows(source, known)
    if len(source) != len(rows):
        raise AssertionError("Screen must retain every raw discovery record")
    with (input_dir / OUTPUT).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(r["screen_queue"] for r in rows)
    matched = {r["doi"]: r["registered_id"] for r in rows if r["registered_id"]}
    missing = sorted(set(known) - set(matched))
    summary = {
        "schema": "ALAN_A_title_triage_v1",
        "source": SOURCE,
        "n_source_unique": len(source), "n_retained": len(rows),
        "queue_counts": dict(counts),
        "known_registry_size": len(known),
        "known_doi_recovered": len(matched),
        "known_doi_missing": missing,
        "unregistered_repro_candidates": [
            {"doi": r["doi"], "title": r["title"], "openalex_id": r["openalex_id"]}
            for r in rows if r["screen_queue"] == "UNREGISTERED_REPRO_MAY_BE_RELEVANT"
        ],
        "non_auto_exclusion": True,
        "fulltext_screened": 0,
        "effect_variance_extracted": 0,
        "meta_regression_beta_V": None,
        "caveat": "Queues use titles only, never establish experimental design or suitability. Every record remains available for full-text assessment."
    }
    (input_dir / "openalex_title_triage_receipt.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    r = triage(pathlib.Path("openalex_a_discovery"), ROOT / "config/alan_A_reproductive_screen_20261115.json")
    print(json.dumps({
        "n_source_unique": r["n_source_unique"], "queue_counts": r["queue_counts"],
        "known_doi_missing": r["known_doi_missing"],
        "fulltext_screened": r["fulltext_screened"],
    }, ensure_ascii=False))
