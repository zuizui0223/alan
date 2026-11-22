"""Rebuild the supplementary citation ledger from the three archived search manifests."""
import csv
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parent / "search"
SOURCES = {
    "meta": ("meta_source_search.json", "candidates"),
    "aquatic": ("aquatic_experiment_search.json", "records"),
    "terrestrial": ("terrestrial_experiment_search.json", "records"),
}
# These AGR entries omit DOIs. The November 22 Europe PMC core metadata (authors,
# venue, year, title and abstract) was individually compared to each primary source.
# This is a local, auditable ID crosswalk, not an automated title-based merge.
INDEX_VERIFICATIONS = {
    "AGR:IND609454928": (
        "10.1016/j.biocon.2026.111711",
        "Five authors in order, journal/volume 315, title and four-site camera-trap experiment agree.",
        "https://pure.knaw.nl/portal/en/datasets/data-from-timing-of-attraction-to-light-of-nocturnal-insects-is-s/",
    ),
    "AGR:IND609426870": (
        "10.1016/j.biocon.2026.111722",
        "Fourteen authors in order, journal, March 2026, title and four-treatment Chihuahuan Desert trial agree.",
        "https://www.sciencedirect.com/science/article/pii/S0006320726000303?dgcid=rss_sd_all",
    ),
    "AGR:IND605878255": (
        "10.1111/gcb.13615",
        "Six authors in order, journal/volume 23/issue 7, July 2017, title and three-year experiment agree.",
        "https://pubmed.ncbi.nlm.nih.gov/28139040/",
    ),
}


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file))


def flatten(rows, fields):
    values = []
    for row in rows:
        for field in fields:
            value = row.get(field)
            if isinstance(value, list):
                values.extend(str(item) for item in value if item)
            elif value:
                values.append(str(value))
    return " | ".join(dict.fromkeys(values))


def collect_sources():
    entries = defaultdict(list)
    manifests = {}
    for name, (filename, section) in SOURCES.items():
        doc = json.loads((HERE / "citation_sources" / filename).read_text(encoding="utf-8"))
        paths = [(f"/{section}/{i}", record) for i, record in enumerate(doc[section])]
        if name == "aquatic":
            paths.append(("/additional_primary_followup", doc["additional_primary_followup"]))
        manifests[name] = {
            "path": f"citation_sources/{filename}",
            "query_list_pointer": "/queries",
            "query_count": len(doc["queries"]),
            "contributing_report_records": len(paths),
        }
        for pointer, raw in paths:
            doi = raw.get("doi", "").strip().lower()
            key = f"doi:{doi}" if doi else f"{name}:{pointer}"
            entries[key].append({
                **raw,
                "source_location": f"citation_sources/{filename}#{pointer}",
                "query_manifest": f"citation_sources/{filename}#/queries",
            })
    return entries, manifests


def main():
    baseline = {row["report_key"].lower(): row for row in read_csv(SEARCH / "unified_review_inventory.csv")}
    records = read_csv(HERE / "records.csv")
    search_by_key = defaultdict(list)
    for record in records:
        record_id = record["record_id"]
        evidence = defaultdict(list)
        if record["doi"]:
            evidence["doi:" + record["doi"].lower()].append("index DOI")
        if record["matched_existing_report_key"]:
            key = record["matched_existing_report_key"].lower()
            if key not in baseline:
                raise ValueError(f"Unrecognized baseline key for {record_id}: {key}")
            evidence[key].append("existing verified record-to-report crosswalk")
        if record_id in INDEX_VERIFICATIONS:
            doi, note, url = INDEX_VERIFICATIONS[record_id]
            key = f"doi:{doi}"
            if evidence and key not in evidence:
                raise ValueError(f"Conflicting identity evidence for {record_id}")
            evidence[key].append(f"verified index ID ({note} {url})")
        for key, basis in evidence.items():
            search_by_key[key].append((record, " / ".join(basis)))

    entries, manifests = collect_sources()
    output = []
    for key, sources in sorted(entries.items()):
        matches = search_by_key.get(key, [])
        in_baseline = key in baseline
        title = flatten(sources, ["title"])
        year = flatten(sources, ["year"])
        if not title and matches:
            title, year = matches[0][0]["title"], matches[0][0]["year"]
        output.append({
            "report_key": key,
            "doi": key.removeprefix("doi:") if key.startswith("doi:") else "",
            "title": title,
            "year": year,
            "source_entries": flatten(sources, ["source_location"]),
            "source_labels": flatten(sources, ["key", "id"]),
            "query_manifests": flatten(sources, ["query_manifest"]),
            "baseline_match": baseline[key]["report_key"] if in_baseline else "",
            "expansion_record_ids": ";".join(r["record_id"] for r, _ in matches),
            "expansion_match_evidence": " | ".join(f"{r['record_id']}: {basis}" for r, basis in matches),
            "incremental_status": ("in_prior_inventory" if in_baseline else
                                   "in_expansion_search_only" if matches else "additional_citation_candidate"),
            "design_and_units": flatten(sources, ["design", "unit", "design_and_units"]),
            "endpoint": flatten(sources, ["endpoint"]),
            "screening_evidence": flatten(sources, ["verification", "screen_status", "status", "data_status"]),
            "location_and_coordinates": flatten(sources, ["site_coordinates", "location", "coordinate_status"]),
            "programme_and_overlap": flatten(sources, ["programme", "identity_warning", "context",
                                                     "family_or_overlap", "possible_programme_overlap"]),
            "primary_urls": flatten(sources, ["primary_url", "primary_sources", "source", "sources"]),
            "data_or_document": flatten(sources, ["data_or_document", "data_or_supplement", "public_data_or_access"]),
            "chat_search_refs": flatten(sources, ["source_ref", "source_refs", "web_refs"]),
        })
    with (HERE / "citation_ledger.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    def count_status(status):
        return sum(row["incremental_status"] == status for row in output)

    summary = {
        "date_jst": "2026-11-22",
        "sources": manifests,
        "appearances_including_followup": sum(len(s) for s in entries.values()),
        "unique_citation_reports": len(output),
        "unique_citation_reports_with_doi": sum(bool(row["doi"]) for row in output),
        "baseline_reports_compared": len(baseline),
        "expansion_europe_pmc_records_compared": len(records),
        "exact_baseline_report_matches": count_status("in_prior_inventory"),
        "matches_to_expansion_including_baseline_reports": sum(bool(row["expansion_record_ids"]) for row in output),
        "matches_to_expansion_outside_baseline": count_status("in_expansion_search_only"),
        "additional_citation_report_candidates_outside_both": count_status("additional_citation_candidate"),
        "additional_doi_candidates_outside_both": sum(
            bool(r["doi"]) and r["incremental_status"] == "additional_citation_candidate" for r in output),
        "checks": {
            "kuhne_icab149": next(r["expansion_record_ids"] for r in output if r["doi"] == "10.1093/icb/icab149"),
            "explicitly_verified_index_ids": INDEX_VERIFICATIONS,
        },
        "interpretation": [
            "Reports are deduplicated by exact case-insensitive article DOI. An undated or non-DOI record would retain its source identity. Distinct papers may share an experiment.",
            "Overlap with Europe PMC uses an exact DOI or verified index-record mapping. Source-level title suggestions are preserved in the archived originals and do not automatically merge reports.",
            "Each source entry points into the unaltered archival manifest with a JSON pointer. The manifests list search queries; they do not attribute each query to a particular record. No per-query hit assignment was reconstructed.",
            "Primary and data URLs provide portable provenance. Chat search references are auxiliary session labels. Publicly linked material may be unavailable for download.",
            "These are supplementary search leads, including ineligible designs, companion papers, and a search-follow-up for Kühne. Counts do not establish extraction readiness, sample size or independent experimental families.",
        ],
    }
    (HERE / "citation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
