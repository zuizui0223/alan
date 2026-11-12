#!/usr/bin/env python3
"""Read a *public original* Dryad moth-trap dataset to audit the light-attraction hypothesis.

This is a SOURCE/POSITIVITY preflight only, not a new biological effect.
Datasets from individual moth traps cannot be treated as evidence of flowers'
effective pollination. Keep the original workbooks outside Git.
"""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path
import urllib.error
import urllib.request
import pandas as pd

DOI="10.5061/dryad.6djh9w1d2"
BASE="https://datadryad.org/api/v2"
DATASET="doi%3A10.5061%2Fdryad.6djh9w1d2"
API=BASE+"/datasets/"+DATASET
ARCHIVE=API+"/download"
EXPECTED_SHEETS={"Sample","Trap","Site","LedBait","Cond","FMSample","FMSample01","Trait"}
DATA_HINTS=("NearestSL","SL400","Skyglow","Trap","BA400","Forest400","Day")


def retrieve(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-Source-Audit/1.0 (open reproducible moth ecology)",
        "Accept":"application/zip,application/json,*/*"})
    with urllib.request.urlopen(req,timeout=100) as r:
        assert r.status==200
        return r.read(),r.headers.get("Content-Type","")


def source_report(path, report_path, optional_api=True):
    receipt={"schema":"ALAN_true_moth_light_vs_bait_trap_original_Dryad_source_audit_v1",
             "source_doi":DOI,"license":"As listed on authoritative Dryad source",
             "source_download_url":ARCHIVE,
             "local_source_sha256":None,"genuine_moth_specimen_rows_claimed":False,
             "not_plant_pollination_nor_global_moth_abundance":True,
             "status":"not_checked"}
    archive=Path(path)
    if archive.is_file() and archive.stat().st_size:
        content=archive.read_bytes()
        receipt["status"]="source_archive_provided"
    else:
        try:
            content, type_header=retrieve(ARCHIVE)
            receipt["source_content_type"]=type_header
            receipt["status"]="public_dryad_archive_downloaded"
        except (OSError, urllib.error.URLError, ValueError) as exc:
            receipt["status"]="source_access_blocked"
            receipt["access_exception"]=str(exc)[:350]
            content=None
            if optional_api:
                try:
                    meta,ctype=retrieve(API)
                    d=json.loads(meta)
                    receipt["metadata_access"]="ok"
                    receipt["metadata_title"]=d.get("title")
                    receipt["metadata_version"]={
                        k:v for k,v in d.items() if k in ("versionNumber","versionStatus",
                                                "publicationDate","identifier","id")}
                    receipt["metadata_keys"]=list(d)[:30]
                except (OSError,urllib.error.URLError,ValueError) as e:
                    receipt["metadata_access"]="blocked: "+str(e)[:200]
    if content is not None:
        receipt["local_source_sha256"]=hashlib.sha256(content).hexdigest()
        receipt["downloaded_bytes"]=len(content)
        try:
            with zipfile.ZipFile(io.BytesIO(content)) as z:
                receipt["download_zip_members"]=[
                    {"name":x.filename,"size":x.file_size}
                    for x in z.infolist() if not x.is_dir()]
                books=[]
                for item in z.infolist():
                    if item.filename.lower().endswith((".xlsx",".xls")):
                        books.append((item.filename,z.read(item)))
                    elif item.filename.lower().endswith(".zip"):
                        with zipfile.ZipFile(io.BytesIO(z.read(item))) as sub:
                            for subitem in sub.infolist():
                                if subitem.filename.lower().endswith((".xlsx",".xls")):
                                    books.append((item.filename+"/"+subitem.filename,
                                                  sub.read(subitem)))
                receipt["found_original_workbooks"]=[]
                for name,data in books:
                    xl=pd.ExcelFile(io.BytesIO(data))
                    book={"filename":name,"sha256":hashlib.sha256(data).hexdigest(),
                          "sheets":[]}
                    for sheet in xl.sheet_names:
                        df=pd.read_excel(xl,sheet_name=sheet)
                        fields=list(map(str,df.columns))
                        book["sheets"].append({
                            "name":sheet,"n_original_rows":len(df),"n_columns":len(fields),
                            "columns":fields[:70],
                            "first_20_dtypes":{str(k):str(v) for k,v in
                                               df.dtypes.head(20).items()}})
                    receipt["found_original_workbooks"].append(book)
                receipt["expected_author_sheets_found"]=[
                    sh for b in receipt["found_original_workbooks"]
                    for sh in EXPECTED_SHEETS if sh in [x["name"] for x in b["sheets"]]]
                receipt["status"]="original_workbook_schema_audited" if books else \
                    "zip_downloaded_no_actual_workbook"
        except (OSError,zipfile.BadZipFile,ValueError) as exc:
            receipt["status"]="nonzip_or_invalid_original_download"
            receipt["archive_parse_error"]=str(exc)[:350]
    rpath=Path(report_path)
    rpath.parent.mkdir(parents=True,exist_ok=True)
    rpath.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("ALAN_MOTH_TRAP_TRUE_SOURCE_PREFLIGHT",json.dumps(receipt,sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--archive",default="source/EVdS_et_al-EnvPol_full_source.zip")
    p.add_argument("--receipt",default="results/alan_true_moth_lightbait_source_preflight.json")
    a=p.parse_args()
    source_report(a.archive,a.receipt)
