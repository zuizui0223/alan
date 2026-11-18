#!/usr/bin/env python3
"""Source-only access to 2015 Dryad actual heathland single-visit pollen network.

Ballantyne et al. (2015) Proc R Soc B DOI 10.1098/rspb.2015.1130
Dryad official dataset DOI 10.5061/dryad.17pp3

Directly matches actual same-community observed visitation and
single-visit stigmatic conspecific pollen (SVD) for Erica tetralix
(wet heath), Erica cinerea, Calluna vulgaris, Ulex minor/europaeus.
The field observation is one Dorset heathland community across 2013-14,
NOT an ALAN manipulation, five independent ecosystems, or replicated
artificial-light plots. Absence of SVD for a visitor is NOT zero PE.
Do not classify Erica tetralix as aquatic/macrophyte or assign current
soil-water content from general species habitat affinity.

Outputs only source provenance, columns, counts and codebook.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import urllib.request
import zipfile
from pathlib import Path

DATASET="doi%3A10.5061%2Fdryad.17pp3"
API="https://datadryad.org/api/v2/datasets/"+DATASET
SCHEMA="alan_ballantyne2015_wet_heath_single_visit_pollen_dryad_original_schema_v1"

def fetch(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"Open-science heathland pollination source audit",
        "Accept":"application/json,application/zip,text/csv,*/*"})
    with urllib.request.urlopen(req,timeout=60) as r:
        return r.read(8<<20),r.geturl(),r.headers.get("Content-Type","")

def read_data(raw):
    # All fields are ecological, not volunteer personal data.
    for enc in ("utf-8-sig","latin-1"):
        try:
            text=raw.decode(enc)
            if text.lower().lstrip().startswith(("<html","<!doctype")):
                raise ValueError("HTML anti-bot masquerading as botanical CSV")
            sep=","
            try:sep=csv.Sniffer().sniff(text[:3000],delimiters=",;\t").delimiter
            except csv.Error:pass
            frame=list(csv.DictReader(io.StringIO(text),delimiter=sep))
            if not frame or len(frame[0])<3:continue
            columns=list(frame[0])
            if not any(any(k in col.lower() for k in (
                "visit","plant","pollen","insect","pollin","flower","taxon","species"))
                       for col in columns):
                continue
            return enc,sep,frame
        except UnicodeError:continue
    raise RuntimeError("No ecological source CSV content found")

def audit():
    m,_,_=fetch(API)
    d=json.loads(m)
    links=d.get("_links",{})
    href=links.get("stash:version",{}).get("href")
    if not href or not href.startswith("/api/v2/versions/"):
        raise RuntimeError("Dryad latest version identifier absent")
    version=int(href.rsplit("/",1)[-1])
    url="https://datadryad.org/api/v2/versions/"+str(version)+"/download"
    payload,actual,mime=fetch(url)
    if not zipfile.is_zipfile(io.BytesIO(payload)):
        raise RuntimeError("Original Dryad full dataset must be verifiable ZIP")
    with zipfile.ZipFile(io.BytesIO(payload)) as z:
        if z.testzip() is not None:raise RuntimeError("Original Dryad ZIP CRC failed")
        names=[v.filename for v in z.infolist() if not v.is_dir()]
        csvs=[k for k in names if k.lower().endswith(".csv")]
        guides=[k for k in names if k.lower().endswith(".txt")]
        if len(csvs)!=1 or not guides:
            raise RuntimeError("Expected exactly one original pollen-and-visit CSV and README")
        name=csvs[0]
        raw=z.read(name)
        enc,sep,frame=read_data(raw)
        guide=z.read(guides[0]).decode("utf-8-sig",errors="replace")
    cols=list(frame[0])
    notes=[s[:150] for s in guide.splitlines() if s.strip()][:24]
    # Metadata only; original individual interactions stay off GitHub.
    summary={
        "schema":SCHEMA,
        "source_doi":"10.5061/dryad.17pp3",
        "article_doi":"10.1098/rspb.2015.1130",
        "official_public_api_version_id":version,
        "source_dataset_title":d.get("title"),
        "source_csv_filename":name,
        "source_csv_sha256":hashlib.sha256(raw).hexdigest(),
        "source_csv_byte_length":len(raw),
        "source_original_rows":len(frame),
        "source_original_header":cols,
        "source_csv_delimiter":sep,
        "source_csv_encoding":enc,
        "source_guidance_filename":guides[0],
        "source_guidance_short_lines":notes,
        "source_valid_zip_CRC":True,
        "source_independent_field_community_count":1,
        "source_ALAN_treatment_manipulated":False,
        "plant_pollen_effectiveness_and_visits_data_present":True,
        "source_no_visitor_absence_inferred_from_missing_single_visit_deposition":True,
        "source_no_independent_light_effect":True,
        "original_sample_rows_not_redistributed":True,
    }
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    data=audit()
    target=Path(args.out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("BALLANTYNE_HEATHLAND_TRUE_POLLEN_VISIT_SOURCE",json.dumps(data,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
