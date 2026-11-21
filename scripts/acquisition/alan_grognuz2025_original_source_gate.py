#!/usr/bin/env python3
"""Original published Grognuz/Gisler/Knop flowering-phenology source availability gate.

Data: DOI 10.17632/xjmn5jhxyn.2, original illuminated vs dark Swiss
wildflower strips (14 paired locations). This script reads public Mendeley
metadata and optionally a small number of public original files. It DOES NOT
fit plant-level or site-level effects until source period, treatment assignment,
species name, response, and true installation identity have been independently
verified. Never use number of frames/photos as independent replicate count.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

DATASET="xjmn5jhxyn"
VERSION=2
ARTICLE_DOI="10.1016/j.biocon.2025.111650"
SOURCE_DOI="10.17632/xjmn5jhxyn.2"
API="https://api.data.mendeley.com"
MAX_EXAMINED_FILE=25*1024*1024
TIMEOUT=50

def get(url):
    req=urllib.request.Request(url,headers={
        "Accept":"application/json, application/vnd.mendeley-public-dataset.1+json",
        "User-Agent":"ALAN original ecology source audit (GitHub research workflow)"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=TIMEOUT) as r:
                return r.read(),r.geturl(),r.headers.get("Content-Type","")
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            if isinstance(e,urllib.error.HTTPError) and e.code in (401,403,404,451):
                raise RuntimeError(f"source_http_{e.code}") from e
            if attempt==2:raise
            time.sleep(2*(attempt+1))

def json_request(url):
    blob,_,content_type=get(url)
    assert len(blob)<5*1024*1024, "Unexpectedly large source metadata"
    return json.loads(blob)

def flatten_files(payload):
    if isinstance(payload,list):return payload
    if isinstance(payload,dict):
        for key in ("files","results","data","items"):
            if isinstance(payload.get(key),list):return payload[key]
    raise ValueError("Unknown Mendeley public files API schema")

def inspect_file(blob,name):
    lower=name.lower()
    if lower.endswith(".zip") and zipfile.is_zipfile(io.BytesIO(blob)):
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            membernames=z.namelist()
            return {"type":"zip","members":len(membernames),
                    "member_descriptions":[{"filename":q.filename,"uncompressed_size":q.file_size}
                                           for q in z.infolist()[:100]]}
    if lower.endswith((".csv",".tsv",".txt")):
        head=blob[:400000].decode("utf-8-sig",errors="replace")
        delimiter="\t" if lower.endswith(".tsv") else ","
        rows=list(csv.reader(io.StringIO(head),delimiter=delimiter))
        return {"type":"delimited","headers":rows[0] if rows else [],
                "first_nonempty_field_counts":[len(x) for x in rows[1:6]],
                "line_estimate":blob.count(b"\n")}
    if lower.endswith(".gz"):
        try:
            head=gzip.decompress(blob[:])
            return {"type":"gzip","uncompressed_bytes":len(head),
                    "header_bytes_preview":head[:90].decode("utf-8",errors="replace")}
        except Exception:
            return {"type":"gzip","status":"unreadable_compressed_file"}
    return {"type":"other","extension":Path(name).suffix.lower()}

def run(try_data=False):
    metadata_url=f"{API}/datasets/{DATASET}?version={VERSION}"
    files_url=f"{API}/datasets/publics/{DATASET}/files?version={VERSION}"
    out={"schema":"alan_grognuz2025_original_source_gate_v1",
         "source_article_doi":ARTICLE_DOI,"source_dataset_doi":SOURCE_DOI,
         "question":"Does verified experimental ALAN flowering shift differ by plant species MOIST?",
         "treatment_installations_from_article":"14 lit strips + 14 paired dark controls",
         "independent_light_installations_are_not_images_or_species":True,
         "metadata_endpoints":[metadata_url,files_url],"downloaded_raw_data_not_committed":True,
         "stage":"source_only"}
    try:
        meta=json_request(metadata_url)
        fs=flatten_files(json_request(files_url))
    except Exception as e:
        out.update(status="source_metadata_stop",reason=str(e))
        return out
    doi=str((meta.get("doi") or {}).get("id") or meta.get("doi") or "")
    out.update(metadata_name=meta.get("name"),metadata_id=meta.get("id"),
               metadata_version=meta.get("version"),metadata_doi=doi,
               metadata_files_count=len(fs))
    if not fs:
        out.update(status="source_files_empty_stop",reason="No published files returned")
        return out
    if doi and SOURCE_DOI not in doi and DATASET not in doi:
        out.update(status="source_doi_mismatch_stop")
        return out
    rows=[]
    for f in fs:
        det=f.get("content_details") or {}
        fid=str(f.get("id") or "")
        name=str(f.get("filename") or f.get("name") or "")
        size=det.get("size",f.get("size"))
        row={"id":fid,"filename":name,"source_size_bytes":size,
             "source_sha256":det.get("sha256_hash"),"status":"published_file_metadata"}
        if try_data and fid and size is not None and int(size)<=MAX_EXAMINED_FILE:
            url=f"{API}/datasets/{DATASET}/files/{fid}/file_downloaded?version={VERSION}"
            try:
                buf,_,_=get(url)
                calc=hashlib.sha256(buf).hexdigest()
                row["fetched_sha256"]=calc
                if det.get("sha256_hash") and det["sha256_hash"].lower()!=calc.lower():
                    row["status"]="hash_mismatch_stop"
                else:
                    row["status"]="source_bytes_verified"
                    row["layout"]=inspect_file(buf,name)
            except Exception as e:
                row.update(status="download_stop",reason=str(e))
        elif try_data:
            row["status"]="size_over_25mb_or_unreported_skip"
        rows.append(row)
    out["source_files"]=rows
    out["status"]="source_metadata_verified"
    out["next_inference_gate"]={
        "site_pair_id":"required original treatment installation identifier",
        "per_plant_taxon":"true species name and precise taxonomic join",
        "response":"flower start and peak DOY / floral resources by species",
        "before_effect_inspection":"validate exact lit/dark pair structure and full 2022/2023 source years",
        "wet_trait":"match independently sourced species EIV MOIST",
        "minimum":"at least 4 wet and 4 drier species with lit and dark paired observations",
        "if_gate_fails":"STOP; retain original study as evidence without moderator coefficient"}
    return out

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--inspect-small-public-files",action="store_true")
    a.add_argument("--output",required=True)
    z=a.parse_args()
    result=run(z.inspect_small_public_files)
    p=Path(z.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="source_files"},ensure_ascii=False,sort_keys=True))
    print("SOURCE_FILES",json.dumps(result.get("source_files",[]),ensure_ascii=False,sort_keys=True))

if __name__=="__main__":main()
