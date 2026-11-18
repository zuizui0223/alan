#!/usr/bin/env python3
"""Identify original 2021 Page et al. 1564-visit SVE Figshare research project data.

Original paper DOI 10.1002/ajb2.1764, Page/Nicholson et al. AJB 2021.
Publication DATA AVAILABILITY explicitly names Figshare project 120921.
Original source reported 168 studies, 1564 effects, 240 plant taxa;
frequency-matched observations from 69 studies. Need to inspect file
columns and explicit insect order FIRST to determine if moth vs fly SVD
on the same plant is available. Never assert these are direct contrasts
until authentic source and metadata are inspected.

Source-only: public metadata/file name/size/sha256 and exploratory structure
of separately downloaded original CSV where below 6 MB. No raw study
observation rows committed into GitHub. Study IDs and organism taxonomy
are public ecological metadata, source file remains in CI /tmp only.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import urllib.error
import urllib.request
from pathlib import Path

PROJECT=120921
DOI="10.1002/ajb2.1764"
API="https://api.figshare.com/v2"
MAX_SOURCE=6<<20
SCHEMA="alan_page2021_single_visit_pollen_1564_figshare_source_gate_v1"

def get(url):
    req=urllib.request.Request(url,headers={
     "User-Agent":"Original pollinator efficacy open-data source audit 2026",
     "Accept":"application/json,text/csv,*/*"})
    with urllib.request.urlopen(req,timeout=60) as r:
        return r.read(MAX_SOURCE+1),r.geturl(),r.headers.get("Content-Type","")

def jget(url):
    blob,redir,mime=get(url)
    if len(blob)>MAX_SOURCE:raise ValueError("Source metadata too large")
    return json.loads(blob)

def inspect_text_csv(buf):
    for enc in ("utf-8-sig","latin-1"):
        try:
            string=buf.decode(enc)
            if string.strip().lower().startswith(("<!doctype","<html")):
                return {"status":"original_url_downloaded_HTML_not_source"}
            dialect=","
            try:dialect=csv.Sniffer().sniff(string[:5000],delimiters=",;\t").delimiter
            except csv.Error:pass
            rows=list(csv.DictReader(io.StringIO(string),delimiter=dialect))
            if not rows or len(rows[0])<5:continue
            keys=list(rows[0])
            if not any(any(t in k.lower() for t in (
                  "study","plant","poll","species","single","sve","visit",
                  "order","family","effect")) for k in keys):continue
            counts={}
            for key in keys:
                if not any(z in key.lower() for z in ("order","guild","group","insect","type","taxon")):
                    continue
                vals={r.get(key,"").strip() for r in rows}
                counts[key]={"n_levels":len(vals),"levels":sorted(vals)[:50]}
            return {"status":"original_ecological_csv_verified",
                "source_rows":len(rows),"source_columns":keys,
                "source_taxonomic_group_fields":counts,
                "encoding":enc,"delimiter":dialect}
        except (UnicodeError,ValueError):
            continue
    return {"status":"unknown_original_non_csv_format"}

def audit():
    out={"schema":SCHEMA,"paper_DOI":DOI,
        "project":PROJECT,
        "original_abstract_reported_studies":168,
        "original_abstract_reported_single_visit_effects":1564,
        "original_abstract_reported_plants":240,
        "original_abstract_reported_frequency_matched_studies":69,
        "moth_vs_fly_same_plant_original_contrasts_estimated":False,
        "published_source_is_an_original_effect_meta_not_an_ALAN_treatment":True,
        "source_requests":[],"status":"not_yet_accessed"}
    try:
        project=jget(f"{API}/projects/{PROJECT}")
        out["original_figshare_project_meta"]={
          "id":project.get("id"),"title":project.get("title"),
          "doi":project.get("doi"),"published_date":project.get("published_date")}
    except Exception as e:
        out.update(status="project_metadata_source_STOP",error_type=type(e).__name__)
        return out
    try:
        articles=jget(f"{API}/projects/{PROJECT}/articles")
        if isinstance(articles,dict):articles=articles.get("items",[])
        if not isinstance(articles,list):raise RuntimeError("Unknown official Figshare articles type")
        out["source_article_count"]=len(articles)
        out["source_article_metadata"]=[]
        for item in articles[:100]:
            aid=item["id"]
            try:detail=jget(f"{API}/articles/{aid}")
            except Exception as e:
                out["source_requests"].append({"article_id":aid,
                        "status":"article_detail_STOP","reason":type(e).__name__})
                continue
            entry={"id":aid,"title":detail.get("title"),
                   "doi":detail.get("doi"),"url":detail.get("url_public_html"),
                   "files":[]}
            for f in detail.get("files",[]):
                x={"id":f.get("id"),"name":f.get("name"),
                   "bytes":f.get("size"),"supplied_MD5":f.get("supplied_md5")}
                url=f.get("download_url")
                if (isinstance(url,str) and
                    (f.get("size") or 0)>0 and
                    int(f.get("size"))<=MAX_SOURCE and
                    x["name"].lower().endswith((".csv",".tsv",".txt"))):
                    try:
                        source,_,_=get(url)
                        if len(source)<=MAX_SOURCE:
                            x["source_sha256"]=hashlib.sha256(source).hexdigest()
                            x["source_schema"]=inspect_text_csv(source)
                        else:
                            x["source_schema"]={"status":"too_large_SKIP"}
                    except Exception as e:
                        x["source_schema"]={"status":"source_file_access_STOP",
                                            "reason":type(e).__name__}
                entry["files"].append(x)
            out["source_article_metadata"].append(entry)
        out["status"]="official_figshare_project_and_file_schema_audited"
    except Exception as e:
        out["status"]="figshare_project_articles_source_STOP"
        out["error_type"]=type(e).__name__
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    d=audit()
    path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("PAGE2021_1564_POLLINATOR_SOURCE_FIGSHARE_PROJECT_AUDIT",
      json.dumps(d,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
