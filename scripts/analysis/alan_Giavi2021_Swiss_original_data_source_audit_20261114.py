#!/usr/bin/env python3
"""Independent Swiss 2016 ALAN / daytime plant-insect source acquisition audit.

Published Giavi, Fontaine & Knop (2021), Nature Communications 12:1690,
DOI 10.1038/s41467-021-22011-8, Zenodo version 4540407.
Original 12 fields (6 lit/6 dark), 37 paired sample events, 21 flowering
plants, 3 pollinating insect orders. Verify true Zenodo-signed MD5 and
read Rdata object schema and original README. No raw source study rows
or synthetic comparisons made here.

Rdata Zenodo source is freely reused with proper CC BY 4.0 attribution;
public artifact remains aggregate-only metadata.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import hashlib
import json
import urllib.request
import urllib.error
import time

DATASET_DOI="10.5281/zenodo.4540407"
ARTICLE_DOI="10.1038/s41467-021-22011-8"
SRC={
 "data.Rdata":"f60b21d82c163cbb53520fda0bc61985",
 "readme":"3fccd78016a9cd41f122907e1561fa19",
 "script.R":"1ba8c516766a879b326d0b14bba69f1e",
}
ENDPOINTS=[
 "https://zenodo.org/records/4540407/files/{name}?download=1",
 "https://zenodo.org/api/records/4540407/files/{name}/content",
 "https://zenodo.org/api/files/404e570d-ef87-4976-bbce-7b3d9c684252/{name}"
]

def load_file(name):
    expected=SRC[name]
    errors=[]
    for template in ENDPOINTS[:2]:
        url=template.format(name=name)
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"alan-source-audit/0.1 (github.com/zuizui0223/alan)",
                "Accept":"*/*",
            })
            with urllib.request.urlopen(req,timeout=25) as resp:
                b=resp.read(3_000_000)
            m=hashlib.md5(b).hexdigest()
            if m!=expected:
                errors.append({"endpoint":url,"error":"md5_mismatch",
                               "returned_md5":m,"source_size":len(b)})
                continue
            return b,{"source_url":url,"original_md5":expected,
                      "retrieved_bytes":len(b),
                      "sha256":hashlib.sha256(b).hexdigest(),
                      "source_authenticated":True}
        except (OSError,ValueError,urllib.error.HTTPError,TimeoutError) as e:
            errors.append({"endpoint":url,"error":type(e).__name__,
                           "http_code":getattr(e,"code",None)})
    return None,{"source_authenticated":False,
                 "original_expected_md5":expected,
                 "attempts":errors}

def dataframe_info(d):
    import pandas as pd
    o={"rows":int(len(d)),
       "columns":[str(x) for x in d.columns],
       "dtype_names":{str(k):str(t) for k,t in d.dtypes.items()},
       "source_missing_fraction_per_col":{
         str(k):round(float(d[k].isna().mean()),4)
         for k in d.columns},
       "candidate_design_fields":{},
       "individual_source_observation_rows_exported":False}
    for col in d.columns:
        name=str(col)
        u=d[col].dropna().astype(str).unique()
        if len(u)<=60 and any(t in name.lower() for t in (
            "treat","plant","spec","guild","group","order","insect",
            "pair","site","plot","sample","survey","day","date",
            "year","light","abund","obs","transect")):
            o["candidate_design_fields"][name]={
               "distinct_count":int(len(u)),
               "levels":sorted(set(u))[:60]}
    return o

def run(out):
    out=Path(out)
    if out.exists():raise FileExistsError("Do not overwrite previous source QA")
    files={}
    original={}
    for name in SRC:
        b,receipt=load_file(name)
        files[name]=receipt
        if b is not None:original[name]=b
    extracted={}
    if "data.Rdata" in original:
        import pyreadr
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"data.Rdata"
            p.write_bytes(original["data.Rdata"])
            try:
                r=pyreadr.read_r(str(p))
                extracted["r_object_schemas"]={
                    str(key):dataframe_info(value) if hasattr(value,"columns") else {
                      "python_type":type(value).__name__}
                    for key,value in r.items()
                }
                extracted["objects_read_success"]=True
            except Exception as e:
                extracted["objects_read_success"]=False
                extracted["error_type"]=type(e).__name__
                extracted["error_string"]=str(e)[:180]
    if "readme" in original:
        # The original tiny README is a public paper methods note. Retain
        # exact source byte hash; extract only the sections about data design.
        r=original["readme"].decode("utf-8",errors="replace")
        extracted["readme_preview"]=r[:6000]
        extracted["readme_original_chars"]=len(r)
    if "script.R" in original:
        s=original["script.R"].decode("utf-8",errors="replace")
        extracted["analysis_script_object_assignments"]=[
            x[:180] for x in s.splitlines()
            if any(k in x for k in ("load(","lmer(","glmer(","read","model","plant","species"))][:65]
    out.parent.mkdir(parents=True,exist_ok=True)
    result={"schema":"ALAN_Giavi2021_Swiss_manipulated_daytime_network_original_source_metadata_v1",
            "study":{"article_DOI":ARTICLE_DOI,"zenodo_DOI":DATASET_DOI,
                     "original_experimental_year":2016,"paired_field_sites":6,
                     "number_field_sites":12,"sampling_paired_events":37,
                     "original_plant_taxa_in_published_analysis":21,
                     "orders":["Diptera","Hymenoptera","Coleoptera"],
                     "original_interaction_events_reported":2384},
            "source_receipts":files,"original_data_structure":extracted,
            "status":"ORIGINAL_ZENODO_DATA_SOURCE_AUTHENTICATED" if
                files.get("data.Rdata",{}).get("source_authenticated") else
                "NO_AUTHENTICATED_ORIGINAL_RDATA_RETRIEVED",
            "guardrails":{"no_unrecorded_plant_inferred_absent":True,
                  "site_number_not_number_independent_plant_studies":True,
                  "not_the_same_region_or_protocol_as_Lu2026_Tibet":True,
                  "no_predation_or_seed_from_Giavi2021_reported":True,
                  "no_ALAN_x_island_shift_phenology_built_from_two_papers":True,
                  "no_original_plant_x_plot_raw_rows_uploaded":True}}
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_GIAVI2021_ORIGINAL_SWISS_RAW_MANIFEST",json.dumps({
        "status":result["status"],"source_receipts":files,
        "original_data_structure":extracted},ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.out)
