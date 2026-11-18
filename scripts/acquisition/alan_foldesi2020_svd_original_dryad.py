#!/usr/bin/env python3
"""Recover published 2020 Dryad single-visit pollen deposition meta-analysis source.

Földesi et al. J Appl Ecol 2021 DOI 10.1111/1365-2664.13798.
Data DOI 10.5061/dryad.05qfttf1n, 28 source studies, 127 pollen
deposition observations for bees (94) and Diptera (33) on 30 plant species.

The original meta includes Apis mellifera controls, NOT Lepidoptera
moths. Therefore the source cannot calibrate q_fly/q_moth directly.
This script fetches a published original file, audits its exact header,
study and species support, but never creates moth coefficients from bees.
No undownloadable-source effect is invented. Public original data stay
in local CI /tmp only; aggregate source schema is committed.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path

DOI="10.5061/dryad.05qfttf1n"
ARTICLE_DOI="10.1111/1365-2664.13798"
SCHEMA="alan_foldesi2020_original_127_SVD_FlyBee_Dryad_source_gate_v1"
SOURCE_URLS=[
 "https://datadryad.org/downloads/file_stream/478779",
 "https://datadryad.org/api/v2/files/478779/download",
 "https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.05qfttf1n/download"
]
API="https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.05qfttf1n"

def fetch(url):
    req=urllib.request.Request(url,headers={
     "User-Agent":"Mozilla/5.0 (compatible; open-science original-source audit)",
     "Accept":"text/csv,application/zip,application/json,*/*"})
    with urllib.request.urlopen(req,timeout=45) as resp:
        payload=resp.read(12<<20)
        return payload,resp.geturl(),resp.headers.get("Content-Type","")

def decode_csv(blob):
    for charset in ("utf-8-sig","utf-8","latin-1"):
        try:
            text=blob.decode(charset)
            if "\x00" in text:continue
            separator=","
            try:separator=csv.Sniffer().sniff(text[:4000],delimiters=",;\t").delimiter
            except csv.Error:pass
            frame=list(csv.DictReader(io.StringIO(text),delimiter=separator))
            header=list(frame[0]) if frame else []
            if len(header)<3:continue
            return header,frame,separator
        except (UnicodeError,ValueError):
            continue
    raise ValueError("Downloaded source not a valid CSV with at least 3 named columns")

def read_candidates(payload):
    candidates=[]
    if zipfile.is_zipfile(io.BytesIO(payload)):
        with zipfile.ZipFile(io.BytesIO(payload)) as z:
            if z.testzip() is not None:raise ValueError("Source ZIP CRC corrupt")
            for info in z.infolist():
                if not info.is_dir() and info.filename.lower().endswith(".csv") and info.file_size<3<<20:
                    candidates.append((info.filename,z.read(info)))
    else:candidates=[("Metadata.csv",payload)]
    return candidates

def summary(csv_name,blob):
    header,rows,delim=decode_csv(blob)
    numeric=[]
    categoric={}
    for k in header:
        values=[str(r.get(k,"")).strip() for r in rows]
        n_nonmissing=sum(v!="" and v.lower() not in ("na","nan","null") for v in values)
        parsed=[]
        for v in values:
            try: parsed.append(float(v))
            except ValueError:pass
        if len(parsed)>=.75*max(n_nonmissing,1):
            numeric.append({"column":k,"numeric_n":len(parsed),
                "min":min(parsed) if parsed else None,"max":max(parsed) if parsed else None})
        else:
            x=Counter(v for v in values if v and v.lower() not in ("na","nan","null"))
            categoric[k]={"nonmissing":n_nonmissing,"distinct_values":len(x),
                "top_source_coded_values":dict(x.most_common(12))
                if len(x)<=30 and max((len(z) for z in x),default=0)<50 else None}
    return {
       "source_filename":csv_name,
       "original_csv_sha256":hashlib.sha256(blob).hexdigest(),
       "original_csv_bytes":len(blob),
       "original_data_rows":len(rows),
       "delimiter":delim,
       "source_columns":header,
       "source_numeric_columns":numeric,
       "source_categorical_column_support":categoric,
       "source_no_moth_coefficients_may_be_invented":True
    }

def audit():
    out={"schema":SCHEMA,"dataset_doi":DOI,
       "article_doi":ARTICLE_DOI,
       "published_meta_studies":28,"published_meta_obs":127,
       "published_fly_observations":33,"published_bee_observations":94,
       "published_plant_species":30,
       "relative_pollen_effectiveness_moth_vs_fly_source_support":False,
       "two_guild_plant_specific_deposition_q_moth_vs_q_fly_estimated":False,
       "no_study_level_source_pseudoreplication_allowed":True,
       "source_requests":[]}
    try:
        raw,u,mime=fetch(API)
        src=json.loads(raw)
        out["source_api"]={
         "http_api_url":API,
         "dataset_title":src.get("title"),
         "identifier":src.get("identifier"),
         "version":src.get("versionNumber"),
         "source_doi_in_payload":src.get("identifier"),
         "metadata_request_bytes":len(raw),
         "available_file_links":str(src.get("_links",{}))[:900]}
    except Exception as e:
        out["source_api"]={"status":"metadata_access_STOP","reason_type":type(e).__name__,
                           "http_status":e.code if isinstance(e,urllib.error.HTTPError) else None}
    for url in SOURCE_URLS:
        test={"url":url}
        try:
            blob,actual,mime=fetch(url)
            test.update(downloaded_bytes=len(blob),content_type=mime,
                        redirected_to=actual if actual!=url else None)
            options=read_candidates(blob)
            summaries=[]
            for name,data in options:
                summaries.append(summary(name,data))
            if not summaries:
                test["status"]="no_csv_in_original_source_STOP"
            else:
                test["status"]="original_csv_source_recovered"
                out["original_source_csvs"]=summaries
                out["status"]="original_SVD_source_schema_recovered"
                out["source_requests"].append(test)
                break
        except urllib.error.HTTPError as e:
            test.update(status="source_http_stop",status_code=e.code)
        except (ValueError,urllib.error.URLError,TimeoutError,zipfile.BadZipFile) as e:
            test.update(status="source_download_or_schema_stop",
                        reason_type=type(e).__name__)
        out["source_requests"].append(test)
    else:
        out["status"]="original_SVD_source_unavailable_STOP"
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    d=audit()
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(d,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("FOLDESI_PUBLIC_SVD_DRYAD_SOURCE_GATE",json.dumps(
      {k:v for k,v in d.items() if k!="source_requests"},ensure_ascii=False,sort_keys=True),
      flush=True)
    print("FOLDESI_DRYAD_REQUESTS",json.dumps(d["source_requests"],sort_keys=True),flush=True)

if __name__=="__main__":
    main()
