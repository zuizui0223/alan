#!/usr/bin/env python3
"""EIDC original Macgregor 2016 relational schema, keys and sampling-support audit.

Original OGL data: DOI 10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708
No missing collection or non-carrier inferred from incomplete key.
Archive ZIP stays ephemeral; only schema, SHA, categories and supports emitted.
"""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

URL="https://data-package.ceh.ac.uk/data/31cc5cec-d33b-4dd6-a932-061ff947e708.zip"
FILES={key:f"data/F{idx}_{name}.csv" for idx,key,name in [
 (1,"sites","sites"),(2,"moths","moths"),(3,"flight","overhead"),(4,"pollen","pollen")]}
MIN_FIELDS={
"sites":{"Site_Number","Site_Name","Number_of_lights_within_100m"},
"moths":{"Date","Site","Method","Treatment","Reference","BinomialName"},
"flight":{"Date","Site","Treatment","FirstCount","SecondCount","ThirdCount"},
"pollen":{"MothReference","BinomialName"}}
EXPECT={"sites":22,"moths":1599,"flight":90,"pollen":361}
SCHEMA="ALAN_EIDC_Macgregor2016_unmodified_relational_sampling_support_v1"

def dig(b):return hashlib.sha256(b).hexdigest()
def kinds(s):
    # Do not encode missing as an environmental category.
    z=s.dropna().astype(str).str.strip()
    freq=z.value_counts(dropna=False)
    return {"n_nonmissing":int(len(z)),"n_unique":int(z.nunique()),
            "frequencies_top_35":{str(k):int(v) for k,v in freq.head(35).items()}}

def normalized(x):
    return x.dropna().astype(str).str.strip()

def run(out):
    out=Path(out)
    if out.exists():raise FileExistsError(out)
    resp=requests.get(URL,timeout=(35,110),headers={"User-Agent":"ALAN-original-moth-field-data-audit/1.0"})
    resp.raise_for_status()
    raw=resp.content
    if not zipfile.is_zipfile(io.BytesIO(raw)): raise ValueError("EIDC actual original source is not a ZIP")
    tables={}
    hashes={}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip(): raise ValueError("Original source CRC invalid")
        for key,name in FILES.items():
            b=z.read(name)
            hashes[key]={"source_filename":name,"file_sha256":dig(b),"file_bytes":len(b)}
            d=pd.read_csv(io.BytesIO(b))
            if not MIN_FIELDS[key].issubset(d.columns) or len(d)!=EXPECT[key]:
                raise ValueError(f"Original EIDC file drift: {name}, n={len(d)}, columns={d.columns.tolist()}")
            tables[key]=d
    sites,moths,flight,pollen=[tables[x] for x in ("sites","moths","flight","pollen")]
    summary={"schema":SCHEMA,"original_data_url":URL,
             "source_DOI":"10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708",
             "source_outer_zip_SHA256":dig(raw),"source_zip_bytes":len(raw),
             "original_inner_filenames_and_SHA":hashes,
             "exact_raw_rows":{key:len(tables[key]) for key in FILES},
             "real_source_relational_tables":{}}
    for key,frame in tables.items():
        d={"source_columns":list(map(str,frame.columns)),"n_original_rows":len(frame),
            "n_all_duplicate_rows":int(frame.duplicated().sum())}
        for field in ("Site","Site_Number","Method","Treatment","Date","Group","MothReference"):
            if field in frame.columns:
                d[field]=kinds(frame[field])
        if key=="moths":
            d["reference_unique"]=int(frame.Reference.nunique())
            d["reference_missing"]=int(frame.Reference.isna().sum())
            d["distinct_site_date_method_treatment"]=int(frame[["Site","Date","Method","Treatment"]].drop_duplicates().shape[0])
            d["sample_count_by_method_treatment"]=[
               {"Method":str(a),"Treatment":str(b),"individual_moths":int(c)}
               for (a,b),c in frame.groupby(["Method","Treatment"]).size().items()]
            d["site_treatment_cross_tab"]={
                str((a,b)):int(v)
                for (a,b),v in frame.groupby(["Site","Treatment"]).size().items()}
            d["site_date_with_both_treatments_count"]=int(
                frame.groupby(["Site","Date"]).Treatment.nunique().ge(2).sum())
        if key=="flight":
            cols=["FirstCount","SecondCount","ThirdCount"]
            vals=frame[cols].apply(pd.to_numeric,errors="coerce")
            d["three_count_missing"]={str(c):int(vals[c].isna().sum()) for c in cols}
            d["total_overhead_moth_passes_recorded"]=int(vals.sum().sum())
            d["cross_site_dates_with_both_treatments"]=int(
                frame.groupby(["Site","Date"]).Treatment.nunique().ge(2).sum())
            d["original_overhead_counts_by_treatment"]={
                str(t):{"n_rows":len(g),
                  "sum_three_reading_counts":int(g[cols].apply(pd.to_numeric,errors="coerce").sum().sum())}
                for t,g in frame.groupby("Treatment")}
        if key=="pollen":
            d["n_pollen_source_references"]=int(frame.MothReference.nunique())
            d["n_pollen_rows_without_matching_F2_reference"]=int(
                (~frame.MothReference.astype(str).isin(moths.Reference.astype(str))).sum())
            counts=frame.drop(columns=["MothReference","BinomialName"]).apply(pd.to_numeric,errors="coerce").fillna(0)
            d["n_pollen_types_with_any_count"]=int(counts.gt(0).any().sum())
            d["sum_original_identified_pollen_grains"]=float(counts.sum().sum())
            d["number_of_nonzero_pollen_grain_rows"]=int(counts.gt(0).any(axis=1).sum())
        summary["real_source_relational_tables"][key]=d
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("ALAN_EIDC_MACGREGOR_ORIGINAL_SOURCE_KEY_AUDIT",json.dumps(summary,sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",default="source-audit/macgregor_2016_original_keys.json")
    args=parser.parse_args()
    run(args.out)
