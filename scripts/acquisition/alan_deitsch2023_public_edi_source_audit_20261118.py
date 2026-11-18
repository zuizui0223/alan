#!/usr/bin/env python3
"""Read original EDI knb-lter-hbr.368.2, audit provenance and source schema.

Only aggregate metadata, hashes, column labels and number of independent
reported plot IDs enter version control/CI artifacts. The source's actual
site-level and caterpillar-level records remain ephemeral in Actions.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

PACKAGE = "knb-lter-hbr.368.2"
BASE = "https://pasta.lternet.edu/package"
TIMEOUT = 40
MAX_READ = 12_000_000


def get(url):
    req = urllib.request.Request(url, headers={
        "User-Agent": "ALAN-original-source-audit/2026-11-18 (academic source verification)",
        "Accept": "*/*"
    })
    with urllib.request.urlopen(req, timeout=TIMEOUT) as response:
        data = response.read(MAX_READ + 1)
        if len(data) > MAX_READ:
            raise ValueError("Original data entity exceeds bounded-source audit limit")
        return data


def source_table_summary(content):
    """Never return original source rows, study site IDs or person identifiers."""
    s = content.decode("utf-8-sig", "replace")
    first = s[:2000]
    if "," not in first and "\t" not in first and ";" not in first:
        return {"format": "not_structured_csv"}
    head = s.partition("\n")[0]
    if head.count("\t") > head.count(","):
        sep = "\t"
    elif head.count(";") > head.count(","):
        sep = ";"
    else:
        sep = ","
    try:
        reader = csv.reader(io.StringIO(s), delimiter=sep)
        columns = [c.strip() for c in next(reader)]
        if len(columns) < 2 or len(columns) > 100:
            return {"format": "not_structured_csv"}
        # Keep schema, raw row count and *number* of distinct physical IDs;
        # do not export original source rows, site IDs or raw treatments.
        names = [v.lower().replace("_", "").replace(" ", "") for v in columns]
        plot_cols = [i for i, k in enumerate(names) if k in (
            "plot", "plotid", "site", "siteid", "treatmentplot", "plotnumber")]
        treat_cols = [i for i, k in enumerate(names) if k in (
            "treatment", "treat", "light", "illumination", "experimental")]
        summary = {"format": "delimited_text", "header_columns": columns,
                   "record_rows": 0, "possible_plot_ID_column_count": len(plot_cols),
                   "possible_light_treatment_column_count": len(treat_cols)}
        groups = [set() for i in plot_cols]
        treatment_levels = [set() for i in treat_cols]
        for row in reader:
            if len(row) != len(columns):
                continue
            summary["record_rows"] += 1
            for target, ix in zip(groups, plot_cols):
                if row[ix].strip():
                    target.add(row[ix].strip())
            for target, ix in zip(treatment_levels, treat_cols):
                if row[ix].strip():
                    target.add(row[ix].strip())
        summary["unique_physical_plot_ID_counts_only"] = [len(q) for q in groups]
        summary["treatment_levels_counts_only"] = [len(q) for q in treatment_levels]
        return summary
    except (csv.Error, StopIteration, IndexError):
        return {"format": "not_structured_csv"}



def figshare_supplement_inspection(report):
    """Cross-check official Figshare supplementary collection; not EDI raw proof."""
    base = "https://api.figshare.com/v2"
    request_url = base + "/collections/6440235/articles"
    try:
        items = json.loads(get(request_url))
        if not isinstance(items, list):
            raise ValueError("Figshare collection was not a public article array")
        report["figshare_source"] = {
            "status": "PUBLIC_COLLECTION_METADATA_RETRIEVED",
            "collection": "10.6084/m9.figshare.c.6440235",
            "articles_listed": len(items), "article_details": [],
        }
        for item in items[:8]:
            article_id = item.get("id")
            if not isinstance(article_id, int):
                continue
            d = json.loads(get(f"{base}/articles/{article_id}"))
            title = d.get("title", "")
            match_source = "caterpillar" in title.lower() or "light at night" in title.lower()
            files = d.get("files", [])
            entry = {
                "article_id": article_id, "title": title, "paper_subject_match": match_source,
                "source_files": [
                    {"name": v.get("name"), "size": v.get("size"), "reported_md5": v.get("supplied_md5")}
                    for v in files
                ],
            }
            report["figshare_source"]["article_details"].append(entry)
            if not match_source:
                continue
            for v in files[:15]:
                name = v.get("name", "")
                size = v.get("size", -1)
                url = v.get("download_url")
                if not isinstance(size, int) or not 0 < size < 4_000_000 or not url:
                    continue
                if not name.lower().endswith((".csv", ".tsv", ".txt")):
                    continue
                try:
                    content = get(url)
                    digest = {"entity_source_url": url, "original_source_file_name": name,
                              "raw_SHA256": hashlib.sha256(content).hexdigest(),
                              "byte_count": len(content), "source_type": "FIGSHARE_SUPPLEMENT_NOT_EDI"}
                    digest.update(source_table_summary(content))
                    report["figshare_source"].setdefault("small_file_schemas", []).append(digest)
                except (ValueError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
                    report["figshare_source"].setdefault("small_file_SOURCE_STOPs", []).append(
                        {"file_name": name, "error_class": type(e).__name__})
    except (ValueError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as e:
        report["figshare_source"] = {
            "status": "SUPPLEMENT_SOURCE_STOP", "url": request_url, "error_class": type(e).__name__}

def inspect_source():
    report = {"schema": "alan_edi_deitsch2023_authentic_source_preflight_20261118_v1",
              "source": {"package_ID": PACKAGE, "dataset_DOI": "10.6073/pasta/71525c9c724db0a21249584afca9ddbd",
                         "original_article_DOI": "10.1098/rspb.2023.0153",
                         "field_design_from_primary_paper": "4 paired ALAN and dark plots (8 total); repeated arthropods and clay caterpillars are nested observations"},
              "resource_requests": {}, "entity_summaries": [],
              "raw_plot_rows_committed": False, "raw_coordinates_committed": False,
              "one_study_not_many_effects": True, "new_model_of_predation_mediation_fitted": False}
    paths = {"resource_map": f"{BASE}/eml/knb-lter-hbr/368/2",
             "entity_names": f"{BASE}/name/eml/knb-lter-hbr/368/2",
             "metadata_eml": f"{BASE}/metadata/eml/knb-lter-hbr/368/2"}
    downloaded = {}
    for key, url in paths.items():
        try:
            data = get(url)
            downloaded[key] = data
            report["resource_requests"][key] = {
                "status": "SOURCE_BYTES_RETRIEVED", "url": url, "sha256": hashlib.sha256(data).hexdigest(),
                "byte_count": len(data)}
        except (ValueError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            report["resource_requests"][key] = {
                "status": "SOURCE_STOP", "url": url, "error_class": type(e).__name__}
    table_info = []
    if "metadata_eml" in downloaded:
        try:
            root = ET.fromstring(downloaded["metadata_eml"])
            for table in root.findall(".//{*}dataTable"):
                name = table.findtext(".//{*}entityName", default="unlabelled")
                attrs = [t.text for t in table.findall(".//{*}attributeName") if t.text]
                table_info.append({"table_name": name, "published_attribute_names": attrs})
        except ET.ParseError:
            report["metadata_XML_source_parse"] = "STOP_MALFORMED_SOURCE"
    report["published_EML_data_tables"] = table_info
    urls = set()
    for key in ("resource_map", "entity_names"):
        raw = downloaded.get(key, b"").decode("utf-8", "replace")
        for s in re.findall(r"https?://[^\s<>\"']+", raw):
            s = s.rstrip(",;)")
            if re.search(r"/package/data/eml/knb-lter-hbr/368/2/[a-zA-Z0-9_-]+$", s):
                urls.add(s)
    if not urls and "resource_map" in downloaded:
        report["entity_listing_no_valid_source_urls"] = True
    for url in sorted(urls)[:20]:
        try:
            b = get(url)
            d = {"entity_source_url": url,
                 "raw_SHA256": hashlib.sha256(b).hexdigest(), "byte_count": len(b)}
            if len(b) <= 4_000_000:
                d.update(source_table_summary(b))
            else:
                d["format"] = "large_entity_no_full_schema_decode"
            report["entity_summaries"].append(d)
        except (ValueError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            report["entity_summaries"].append({"entity_source_url": url, "status": "SOURCE_STOP",
                                               "error_class": type(e).__name__})
    report["source_status"] = (
        "RAW_SOURCE_AND_SCHEMA_VERIFIED" if any(v.get("format") == "delimited_text"
                                             for v in report["entity_summaries"])
        else "SOURCE_METADATA_ONLY_OR_STOP")
    figshare_supplement_inspection(report)
    report["independent_raw_plot_effect_fit"] = False
    report["no_invented_se_or_light_pair_ids"] = True
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    d = inspect_source()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(d, ensure_ascii=False, indent=2) + "\n")
    print("DEITSCH2023_EDI_SOURCE", d["source_status"])
    print("DEITSCH2023_EDI_ENDPOINT_STATUSES",
          [(k, v["status"]) for k, v in d["resource_requests"].items()])
    print("DEITSCH2023_EDI_PUBLIC_ENTITY_SCHEMAS", len(d["entity_summaries"]))


if __name__ == "__main__":
    main()
