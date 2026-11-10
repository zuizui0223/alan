#!/usr/bin/env python3
"""Acquire audited public field data for opposing local vs broad ALAN effects.

DO NOT model observations without verifying original units/replicate keys.
Source files are transient; only source checksums and schema counts public.
"""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path
import pandas as pd
import requests

RESOURCES=(
 ("macgregor_2016_eidc","https://data-package.ceh.ac.uk/data/31cc5cec-d33b-4dd6-a932-061ff947e708",None),
 ("grenis_2023_zenodo","https://zenodo.org/records/8349596/files/Grenis_Moth%20Species%20Compiled.xlsx?download=1","896366aabd01efa61fbea68bd20e1bbb"),
)
SOURCE_DOCS={
 "macgregor_2016_eidc":"https://doi.org/10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708",
 "grenis_2023_zenodo":"https://doi.org/10.5281/zenodo.8349596"
}

def sha256(b):return hashlib.sha256(b).hexdigest()
def md5(b):return hashlib.md5(b).hexdigest()
def pull(session,url):
    response=session.get(url,timeout=(35,110),allow_redirects=True)
    response.raise_for_status()
    return response

def audit_table(raw, name):
    if name.lower().endswith(".csv"):
        try: d=pd.read_csv(io.BytesIO(raw),nrows=10000,low_memory=False)
        except Exception: return {"name":name,"status":"not_parseable_csv","file_bytes":len(raw)}
        return {"name":name,"status":"readable_csv_first_10000","columns":list(d.columns)[:90],
                "preview_read_rows":len(d),"file_bytes":len(raw)}
    if name.lower().endswith((".xlsx",".xls")):
        out={"name":name,"file_bytes":len(raw),"sheets":[]}
        with pd.ExcelFile(io.BytesIO(raw)) as x:
            for sh in x.sheet_names:
                d=pd.read_excel(x,sheet_name=sh)
                out["sheets"].append({"name":sh,"rows":len(d),"columns":list(map(str,d.columns))[:100],
                                      "min_missing_percent":float(d.isna().mean().min()*100) if len(d.columns) else None,
                                      "max_missing_percent":float(d.isna().mean().max()*100) if len(d.columns) else None})
        return out
    return {"name":name,"status":"other_file","file_bytes":len(raw)}

def run(outdir):
    outdir=Path(outdir)
    outdir.mkdir(parents=True,exist_ok=True)
    session=requests.Session()
    session.headers.update({"User-Agent":"ALAN-open-ecology-original-source-2026/1.0",
                            "Accept":"application/zip,application/octet-stream,*/*"})
    result={"schema":"ALAN_independent_moth_field_light_bait_and_pollen_raw_data_audit_20261110_v1",
            "no_empirical_biological_model_fitted":True,"source_results":{}}
    for name,url,expected_md5 in RESOURCES:
        d={"publication":SOURCE_DOCS[name],"download_url":url}
        try:
            response=pull(session,url)
            b=response.content
            d["status"]="original_source_downloaded"
            d["bytes"]=len(b);d["sha256"]=sha256(b);d["md5"]=md5(b)
            d["http_status"]=response.status_code
            d["content_type"]=response.headers.get("Content-Type")
            d["final_url"]=response.url
            if expected_md5:
                d["expected_published_md5"]=expected_md5
                if md5(b)!=expected_md5:
                    d["status"]="published_original_MD5_mismatch"
                    result["source_results"][name]=d
                    continue
            if len(b)<150 or b[:100].lstrip().startswith(b"<!DOCTYPE"):
                d["status"]="unexpected_html_or_empty"
            elif zipfile.is_zipfile(io.BytesIO(b)):
                d["archive_members"]=[]
                with zipfile.ZipFile(io.BytesIO(b)) as z:
                    for member in z.infolist():
                        if member.is_dir():continue
                        v={"name":member.filename,"file_bytes":member.file_size}
                        if member.file_size<=10_000_000 and member.filename.lower().endswith((".csv",".xlsx",".xls")):
                            v["table"]=audit_table(z.read(member),member.filename)
                        d["archive_members"].append(v)
                d["status"]="verified_original_archive_and_schema"
            elif b[:2]==b"PK" or url.endswith("xlsx"):
                d["direct_workbook"]=audit_table(b,name+".xlsx")
                d["status"]="verified_original_workbook_and_schema"
            elif b[:1]==b"{":
                d["status"]="unexpected_JSON_response"
            else:
                d["status"]="unrecognized_download"
        except (requests.RequestException,zipfile.BadZipFile,ValueError,OSError,KeyError) as e:
            d["status"]="source_unavailable_or_unverified"
            d["exception"]=str(e)[:450]
        result["source_results"][name]=d
        print("ALAN_ORIGINAL_INDEPENDENT_MOTH_SOURCE",json.dumps({name:d},sort_keys=True),flush=True)
    (outdir/"real_moth_public_original_schema_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,default=Path("source-audit"))
    a=parser.parse_args()
    run(a.out)
