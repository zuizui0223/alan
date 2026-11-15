#!/usr/bin/env python3
"""Original Figshare registered source probe, not an ALAN manipulation.

Blareau & Requier 2025 online / 2026 journal print:
DOI 10.1111/jen.70019 ; original data DOI 10.6084/m9.figshare.29930195.v1 .
Only source provenance, file kinds and safe schema are collected. No individual
plant/site observations are committed. Public download statuses are factual.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import pathlib
import urllib.error
import urllib.request
from datetime import datetime, timezone

FIGSHARE_ARTICLE_ID = 29930195
API = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE_ID}"
DOI_PAPER = "10.1111/jen.70019"
DOI_DATA = "10.6084/m9.figshare.29930195.v1"
MAX_FILE_BYTES = 20_000_000

def get_public(url: str, max_bytes: int = MAX_FILE_BYTES) -> tuple[dict, bytes]:
    meta = {"url": url, "status": None, "bytes": None, "sha256": None, "error": None}
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ALAN-SOURCE-AUDIT-2026-11-15 (https://github.com/zuizui0223/alan)"})
        with urllib.request.urlopen(req, timeout=35) as resp:
            data = resp.read(max_bytes + 1)
            meta["status"] = resp.status
            if len(data) > max_bytes:
                meta["error"] = "OVER_MAX_FILE_BYTES"
                return meta, b""
            meta["bytes"] = len(data)
            meta["sha256"] = hashlib.sha256(data).hexdigest()
            return meta, data
    except urllib.error.HTTPError as exc:
        meta["status"] = exc.code
        meta["error"] = "HTTP_" + str(exc.code)
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        meta["error"] = type(exc).__name__ + ": " + str(exc)[:200]
    return meta, b""

def csv_schema(data: bytes) -> dict:
    for encoding in ("utf-8-sig", "latin-1"):
        try:
            text = data.decode(encoding)
            reader = csv.reader(io.StringIO(text))
            cols = next(reader, [])
            n = sum(1 for row in reader if any(cell.strip() for cell in row))
            return {"file_encoding": encoding, "columns": cols, "row_count_nonempty": n}
        except UnicodeDecodeError:
            continue
    return {"error":"UNKNOWN_ENCODING"}

def discover() -> dict:
    metadata_receipt, metadata_data = get_public(API, max_bytes=5_000_000)
    report = {
        "schema":"ALAN_P1_B1_original_fruit_crop_2025_figshare_source_manifest_v1",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "parent_article_doi":DOI_PAPER,
        "data_doi":DOI_DATA,
        "source_type":"NIGHT_DAY_EXCLOSURE_ONLY_NOT_ALAN_MANIPULATION",
        "study_calendar_years":[2020,2024],
        "known_published_species":["Rubus idaeus","Ribes rubrum","Ribes nigrum"],
        "published_total_flowers":850,
        "published_total_plants":36,
        "paper_source":"https://onlinelibrary.wiley.com/doi/full/10.1111/jen.70019",
        "dataset_metadata_receipt":metadata_receipt,
        "files":[],
        "moderator_V_source_usable":None,
        "ALAN_effect_estimated":False,
        "claim":"Source discovery only. No source-data analysis is asserted without compatible fields and independent sampling units."
    }
    if metadata_data:
        try:
            metadata=json.loads(metadata_data)
            report["figshare_metadata_title"]=metadata.get("title")
            report["figshare_metadata_doi"]=metadata.get("doi")
            report["license"]=metadata.get("license")
            report["number_of_public_files"]=len(metadata.get("files") or [])
            for f in metadata.get("files") or []:
                entry={"id":f.get("id"),"name":f.get("name"),"size_recorded":f.get("size"),"md5_recorded":f.get("supplied_md5")}
                url=f.get("download_url")
                if url and int(f.get("size") or 0)<=MAX_FILE_BYTES:
                    status,data=get_public(url)
                    entry["receipt"]=status
                    if data and f.get("name","").lower().endswith(".csv"):
                        entry["csv_data_structure"]=csv_schema(data)
                    else:
                        entry["original_data_analyzed"]=False
                else:
                    entry["download_not_attempted"]="NO_URL_OR_LARGE_FILE"
                report["files"].append(entry)
        except ValueError as exc:
            report["metadata_parse_error"]=str(exc)
    target=pathlib.Path("results/alan_P1_B1_Bareau_Requier_original_figshare_manifest_20261115.json")
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return report

if __name__=="__main__":
    x=discover()
    print(json.dumps({"api_status":x["dataset_metadata_receipt"]["status"],
                      "files":[{"name":f["name"],"size":f.get("size_recorded"),
                                "download_status":f.get("receipt",{}).get("status"),
                                "csv":f.get("csv_data_structure")}
                               for f in x["files"]],
                      "ALAN_effect_estimated":x["ALAN_effect_estimated"]},ensure_ascii=False))
