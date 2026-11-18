#!/usr/bin/env python3
"""Eckhartt & Ruxton 2022 original Dryad 7-lamppost prey-bait field data gate.

Only DOI metadata, SHA256, sheet titles, column headings, row counts, and
counts of independent lamppost units may leave a runner. Individual bait
locations, plot IDs, outcome rows, coordinates, and raw workbook bytes remain
ephemeral. Do NOT treat individual baits as experimental light replicates.
"""
import argparse
import hashlib
import io
import json
import re
import urllib.error
import urllib.request
from pathlib import Path

DOI = "10.5061/dryad.7m0cfxpxv"
ARTICLE = "10.1098/rsbl.2022.0281"
API = "https://datadryad.org/api/v2"
FIXED_SOURCE_FILE = "https://datadryad.org/downloads/file_stream/1919546"
MAX_BYTES = 2_000_000


def get(url):
    req = urllib.request.Request(url, headers={
        "User-Agent": "ALAN-source-verification/1.0 (scientific replication)",
        "Accept": "application/json,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*",
    })
    with urllib.request.urlopen(req, timeout=25) as f:
        body = f.read(MAX_BYTES + 1)
        if len(body) > MAX_BYTES:
            raise ValueError("Original source exceeds bounded preflight limit")
        return body


def workbook_schema(original):
    from openpyxl import load_workbook
    if not original.startswith(b"PK\x03\x04"):
        raise ValueError("Source bytes are not an XLSX ZIP")
    wb = load_workbook(io.BytesIO(original), read_only=True, data_only=True)
    sheets = []
    for ws in wb.worksheets:
        rows = ws.iter_rows(values_only=True)
        candidates = []
        for i, row in enumerate(rows):
            if i >= 12:
                break
            cells = [str(v).strip() if v is not None else "" for v in row]
            if sum(bool(v) for v in cells) >= 3:
                candidates.append((i+1,cells))
        if candidates:
            header_row, headers = candidates[0]
        else:
            header_row, headers = (None,[])
        names=[x.casefold().replace(" ","").replace("_","").replace(".","") for x in headers]
        possible_plot=[i for i,n in enumerate(names) if re.search(r"(lamp|plot|site|post|pole|station|loc)",n)]
        # The source has seven lampposts but the original column identity has
        # to be read, not guessed. Never print the actual observed IDs.
        distinct = {i:set() for i in possible_plot}
        n=0
        for ridx,row in enumerate(ws.iter_rows(values_only=True),start=1):
            if header_row is None or ridx<=header_row:continue
            if not any(x is not None for x in row):continue
            n+=1
            for i in possible_plot:
                if len(row)>i and row[i] is not None:
                    distinct[i].add(str(row[i]))
        sheets.append({
            "title":ws.title,
            "header_row":header_row,
            "header":headers,
            "source_rows_after_header":n,
            "possible_lamppost_ID_column_indices":[i for i in possible_plot],
            "possible_lamppost_ID_unique_counts_only":{str(i):len(s) for i,s in distinct.items()},
        })
    wb.close()
    return sheets


def audit():
    result={
        "schema":"alan_20261118_eckhartt2022_original_dryad_source_schema_v1",
        "published_article_doi":ARTICLE,
        "frozen_dataset_doi":DOI,
        "published_source_design":{
            "lampposts":7,"study_days":13,"distance_levels":5,
            "published_bait_observations":830,"independent_light_installations_not_830":True,
            "pre_post_dusk_not_randomized_light_on_off":True,
        },
        "fetch_attempts":[],
        "source_status":"SOURCE_STOP",
        "raw_bait_or_lamppost_IDs_exported":False,
        "individual_data_rows_exported":False,
        "causal_alternation_of_light_identified":False,
    }
    doi="doi%3A10.5061%2Fdryad.7m0cfxpxv"
    for label,url in [
        ("dryad_dataset_api",f"{API}/datasets/{doi}"),
        ("dryad_dataset_versions",f"{API}/datasets/{doi}/versions"),
    ]:
        try:
            b=get(url)
            j=json.loads(b)
            result["fetch_attempts"].append({"label":label,"status":"SUCCESS",
                "bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()})
            if label=="dryad_dataset_api":
                result["api_dataset_title"]=str(j.get("title",""))[:220]
                result["api_embedded_version_status"]=str(j.get("versionStatus",""))[:100]
            if label=="dryad_dataset_versions":
                result["api_version_total"]=j.get("total")
        except (ValueError,KeyError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
            result["fetch_attempts"].append({"label":label,"status":"SOURCE_STOP",
                "error_type":type(exc).__name__,
                **({"http_status":exc.code} if isinstance(exc,urllib.error.HTTPError) else {})})
    try:
        original=get(FIXED_SOURCE_FILE)
        result["fetch_attempts"].append({"label":"original_workbook_file_stream_1919546",
            "status":"SUCCESS", "bytes":len(original),"sha256":hashlib.sha256(original).hexdigest()})
        result["source_xlsx_sha256"]=hashlib.sha256(original).hexdigest()
        result["source_original_workbook_bytes"]=len(original)
        result["source_original_workbook_sheets"]=workbook_schema(original)
        result["source_status"]="RAW_XLSX_SCHEMA_VERIFIED"
    except (ValueError,ImportError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as exc:
        result["fetch_attempts"].append({"label":"original_workbook_file_stream_1919546",
            "status":"SOURCE_STOP","error_type":type(exc).__name__,
            **({"http_status":exc.code} if isinstance(exc,urllib.error.HTTPError) else {})})
        if any(x["status"]=="SUCCESS" for x in result["fetch_attempts"]):
            result["source_status"]="API_METADATA_ONLY"
    result["raw_7_lamppost_spatiotemporal_model_fitted"]=False
    result["no_unobserved_precise_effect_estimated"]=True
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    x=audit()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2)+"\n")
    print("ECKHARTT_2022_ORIGINAL_STATUS",x["source_status"])
    print("ECKHARTT_2022_ENDPOINTS",x["fetch_attempts"])
    print("ECKHARTT_2022_SCHEMA",x.get("source_original_workbook_sheets"))


if __name__=="__main__":
    main()
