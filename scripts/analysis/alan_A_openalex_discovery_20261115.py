#!/usr/bin/env python3
"""ALAN A: source-visible OpenAlex title/abstract discovery, not automatic eligibility.

Each query retrieves only the first 100 hits without an API key. If OpenAlex
reports more, mark as truncated; do not claim a systematic census or PRISMA.
The known-study registry is independent of discovery, and DOI/title matches do
not automatically qualify a study for meta-regression.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import pathlib
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

CUTOFF = "2026-11-15"
QUERIES = (
    "\"artificial light at night\" \"fruit set\"",
    "\"artificial light at night\" \"seed set\"",
    "\"artificial light at night\" pollination reproduction",
    "\"artificial light at night\" plant reproduction",
    "\"light pollution\" \"fruit set\" pollination",
    "\"light pollution\" \"seed set\" plant",
    "\"street lighting\" pollination fruit",
    "streetlights moth pollination seed",
    "\"nocturnal lighting\" plant fecundity",
    "\"artificial light at night\" reproductive success",
    "\"artificial light at night\" pollination",
    "\"artificial light at night\" seed",
    "\"artificial light at night\" fruit",
    "\"artificial light at night\" flowering",
    "\"artificial light at night\" ecosystem functioning",
    "\"artificial light at night\" reproductive capacity",
    "\"artificial light at night\" plant-pollinator",
    "\"light pollution\" hawkmoth",
    "\"artificial night lighting\" plant reproductive",
    "\"light pollution\" nocturnal pollination",
)
FIELDS = ("openalex_id", "doi", "title", "year", "publication_date", "query_count", "query_ids")


def normalize_doi(s: str | None) -> str:
    return str(s or "").strip().lower().replace("https://doi.org/", "").replace("http://dx.doi.org/", "").removeprefix("doi:")


def normalize_record(work: dict, query_num: int) -> dict | None:
    if not isinstance(work, dict):
        return None
    pub_date = str(work.get("publication_date") or "")
    if pub_date and pub_date > CUTOFF:
        return None
    doi = normalize_doi(work.get("doi"))
    work_id = str(work.get("id") or "")
    if not doi and not work_id:
        return None
    return {
        "openalex_id": work_id,
        "doi": doi,
        "title": str(work.get("title") or work.get("display_name") or ""),
        "year": work.get("publication_year") or "",
        "publication_date": pub_date,
        "query_count": 1,
        "query_ids": str(query_num),
    }


def merge_rows(query_payloads: list[tuple[int, dict]]) -> list[dict]:
    seen: dict[str, dict] = {}
    for query_num, payload in query_payloads:
        for work in payload.get("results", []):
            row = normalize_record(work, query_num)
            if row is None:
                continue
            key = row["doi"] or row["openalex_id"].lower()
            if key not in seen:
                seen[key] = row
            else:
                old = seen[key]
                old["query_count"] += 1
                old["query_ids"] += "," + str(query_num)
    return sorted(seen.values(), key=lambda x: (str(x["year"]), x["doi"], x["openalex_id"]))


def run_query(text: str) -> dict:
    params = {
        "search.title_and_abstract": text,
        "filter": "to_publication_date:" + CUTOFF,
        "per_page": 100,
        "page": 1,
    }
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
    headers = {"User-Agent": "ALAN-evidence-screen/2026-11-15 (+https://github.com/zuizui0223/alan)"}
    if os.environ.get("OPENALEX_API_KEY"):
        params["api_key"] = os.environ["OPENALEX_API_KEY"]
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
    last_error = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=35) as response:
                data = json.load(response)
            if not isinstance(data.get("results"), list) or not isinstance(data.get("meta"), dict):
                raise ValueError("Malformed OpenAlex response")
            return data
        except (urllib.error.URLError, ValueError, TimeoutError) as exc:
            last_error = exc
            time.sleep(2 + 3 * attempt)
    raise RuntimeError("OpenAlex query failed: " + repr(last_error))


def discover(out_dir: pathlib.Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    payloads, per_query, errors = [], [], []
    for i, query in enumerate(QUERIES, 1):
        try:
            payload = run_query(query)
            meta = payload["meta"]
            n = len(payload["results"])
            total = int(meta.get("count", 0))
            payloads.append((i, payload))
            per_query.append({"query_id": i, "query": query, "total_reported": total,
                              "returned_page_one": n, "truncated_after_first_page": total > n})
        except RuntimeError as exc:
            errors.append({"query_id": i, "query": query, "error": str(exc)})
    rows = merge_rows(payloads)
    with (out_dir / "openalex_firstpage_metadata.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    registry_path = pathlib.Path("config/alan_A_reproductive_screen_20261115.json")
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    known_dois = {normalize_doi(x["doi"]): x["id"] for x in registry["entries"]}
    matched = {x["doi"]: known_dois[x["doi"]] for x in rows if x["doi"] in known_dois}
    result = {
        "schema": "ALAN_A_OpenAlex_first_page_discovery_v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cutoff": CUTOFF,
        "record_status": "SEARCH_DISCOVERY_ONLY_NOT_ELIGIBILITY_OR_PRISMA",
        "query_count": len(QUERIES),
        "queries_succeeded": len(payloads),
        "query_details": per_query,
        "errors": errors,
        "records_unique_within_pages": len(rows),
        "truncated_any": any(x["truncated_after_first_page"] for x in per_query),
        "known_registry_doi_matches": matched,
        "eligible_independent_ALAN_experiments": None,
        "eligible_reproductive_effects_with_se": None,
        "eligible_external_night_dependence": None,
        "full_eligibility_inference_allowed": False,
        "limitation": "Title and abstract retrieval only, first page per query. No systematic full-text, WoS/Scopus, or citation chasing; no automatic causal eligibility."
    }
    (out_dir / "openalex_firstpage_receipt.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="openalex_a_discovery")
    args = parser.parse_args()
    report = discover(pathlib.Path(args.output))
    print(json.dumps({k: report[k] for k in (
        "query_count", "queries_succeeded", "records_unique_within_pages", "truncated_any",
        "eligible_independent_ALAN_experiments", "errors")}, ensure_ascii=False))
    if not report["queries_succeeded"]:
        raise SystemExit("All OpenAlex queries failed; no source scan performed.")


if __name__ == "__main__":
    main()
