"""Supplementary arthropod search; no eligibility or effect estimation.

Usage: python search_snapshot.py /absolute/path/to/private_raw_cache
Raw abstracts remain outside the repository. Re-execution must use a new
snapshot directory; it never silently overwrites an existing manifest.
"""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sys
import urllib.parse
import urllib.request

OUT = Path(__file__).resolve().parent
BASE = OUT.parent / "search"
ENDPOINT = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
QUERY = (
    '(TITLE_ABS:"artificial light at night" OR TITLE_ABS:"light pollution" '
    'OR TITLE_ABS:"street lighting" OR TITLE_ABS:streetlight* '
    'OR TITLE_ABS:"nocturnal illumination" OR TITLE_ABS:"night-time light") '
    'AND (TITLE_ABS:insect* OR TITLE_ABS:arthropod* OR TITLE_ABS:invertebrate* '
    'OR TITLE_ABS:spider* OR TITLE_ABS:moth* OR TITLE_ABS:beetle* '
    'OR TITLE_ABS:Diptera OR TITLE_ABS:Chironomid* OR TITLE_ABS:mayfl* '
    'OR TITLE_ABS:caddisfl*) AND FIRST_PDATE:[1900-01-01 TO 2026-11-22]'
)


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(name, rows, fields):
    with (OUT / name).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    cache = Path(sys.argv[1]).resolve()
    repo_root = OUT.parents[2]
    if repo_root == cache or repo_root in cache.parents:
        raise ValueError("Raw responses must remain outside the public repository")
    if (OUT / "search_manifest.json").exists():
        raise FileExistsError("Snapshot already exists; preserve it and use a new directory")
    cache.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    pages, raw_records, seen = [], [], set()
    cursor, total = "*", None
    for index in range(100):
        if cursor in seen:
            raise ValueError("Repeated cursor")
        seen.add(cursor)
        params = {"query": QUERY, "format": "json", "pageSize": 1000,
                  "resultType": "core", "cursorMark": cursor}
        url = ENDPOINT + "?" + urllib.parse.urlencode(params)
        request = urllib.request.Request(url, headers={"User-Agent": "ALAN-research-evidence-search/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            raw = response.read()
        data = json.loads(raw)
        if total is None:
            total = data["hitCount"]
        if data["hitCount"] != total:
            raise ValueError("Index changed during retrieval")
        path = cache / f"europepmc_{index:03}.json"
        path.write_bytes(raw)
        batch = data["resultList"]["result"]
        pages.append({"page": index, "url": url, "filename": path.name,
                      "sha256": sha(path), "records": len(batch)})
        raw_records.extend(batch)
        print(json.dumps({"page": index, "received": len(batch),
                          "accumulated": len(raw_records), "hit_count": total}), flush=True)
        if len(raw_records) >= total:
            break
        if not batch:
            raise ValueError("Premature empty page")
        cursor = data["nextCursorMark"]
    if len(raw_records) != total:
        raise ValueError("Incomplete retrieval")
    ids = [r["source"] + ":" + r["id"] for r in raw_records]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate index IDs")

    inventory = {r["report_key"] for r in read_csv(BASE / "unified_review_inventory.csv")}
    appearances = read_csv(BASE / "combined_search_appearances.csv")
    by_id = {r["record_id"]: r["report_key"] for r in appearances}
    by_pmc = {}
    for r in appearances:
        if r["pmcid"]:
            by_pmc.setdefault(r["pmcid"], set()).add(r["report_key"])
    records = []
    for r in raw_records:
        ident = r["source"] + ":" + r["id"]
        doi = r.get("doi", "").strip().lower().removeprefix("https://doi.org/")
        key = "doi:" + doi if doi else ident
        matches = set()
        if key in inventory:
            matches.add(key)
        if ident in by_id:
            matches.add(by_id[ident])
        matches.update(by_pmc.get(r.get("pmcid", ""), set()))
        status = "identity_conflict_review" if len(matches) > 1 else ("existing_report" if matches else "unmatched_candidate_report")
        matched = next(iter(matches)) if len(matches) == 1 else ""
        canonical = matched or key
        records.append({"record_id": ident, "doi": doi, "pmcid": r.get("pmcid", ""),
            "title": r.get("title", ""), "year": r.get("pubYear", ""),
            "first_publication": r.get("firstPublicationDate", ""),
            "journal": r.get("journalInfo", {}).get("journal", {}).get("title", ""),
            "source_url": "https://europepmc.org/article/" + r["source"] + "/" + r["id"],
            "has_abstract": bool(r.get("abstractText")), "match_status": status,
            "matched_existing_report_key": matched, "candidate_report_key": canonical,
            "identity_candidates": ";".join(sorted(matches)),
            "screening_status": "prior_history_preserved" if matched else "pending"})
    fields = list(records[0]) if records else ["record_id", "candidate_report_key", "screening_status"]
    write_csv("records.csv", records, fields)
    write_csv("unmatched_records.csv", [r for r in records if r["match_status"] != "existing_report"], fields)
    write_csv("search_pages.csv", pages, ["page", "url", "filename", "sha256", "records"])
    new_keys = {r["candidate_report_key"] for r in records if r["match_status"] == "unmatched_candidate_report"}
    manifest = {"search_id": "arthropod_europepmc_20261122", "database": "Europe PMC",
        "role": "supplementary_index_search", "query": QUERY, "endpoint": ENDPOINT,
        "started_utc": started, "completed_utc": datetime.now(timezone.utc).isoformat(),
        "cutoff": "2026-11-22", "raw_hits": total, "retrieved_records": len(records),
        "existing_reports_in_baseline_inventory": len(inventory),
        "records_matching_baseline": sum(r["match_status"] == "existing_report" for r in records),
        "unmatched_records": sum(r["match_status"] == "unmatched_candidate_report" for r in records),
        "unmatched_candidate_report_keys": len(new_keys),
        "identity_conflict_records": sum(r["match_status"] == "identity_conflict_review" for r in records),
        "records_without_abstract": sum(not r["has_abstract"] for r in records),
        "baseline_input_sha256": {n: sha(BASE/n) for n in ["unified_review_inventory.csv", "combined_search_appearances.csv"]},
        "raw_response_sha256": {p["filename"]: p["sha256"] for p in pages},
        "final_prisma": False, "final_eligibility": False, "independent_experiments_estimated": False,
        "limit": "Additional candidate reports require screening and identity/experiment-family review. Europe PMC is not WoS or Scopus. No abstracts redistributed in metadata CSVs."}
    (OUT / "search_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
