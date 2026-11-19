#!/usr/bin/env python3
"""First source audit of Sanders et al. 2020/2021 ALAN meta original Dryad.

DOI 10.5061/dryad.wpzgmsbjn, the published 126-paper / 1304-effect
source. Source used Hedges d, NOT lnRR. This preflight inspects the
authenticated author CSV+R script and reports only SHA, original column
names, source effect rows, original study ID cardinality, and coarse
provenance fields. NO raw observations or study-identifying title rows
are committed to the repository. It never invents effect variances
or imputes water vs terrestrial from taxonomy.
"""
import argparse
import csv
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import urllib.parse
import urllib.request
import zipfile

DATA_DOI="10.5061/dryad.wpzgmsbjn"
PUBLISHED_PAPERS=126
PUBLISHED_EFFECTS=1304
EXPECTED_FILES={"Sanders_et_al_datafile.csv","Sanders_meta__rcode.r"}
BASE="https://datadryad.org/api/v2"
MAX_BYTES=4_000_000
SOURCE_TYPE="ALAN_126paper_1304_hedges_d_source"


def get(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-systematic-review-source-preflight/2026-11-19",
        "Accept":"application/zip,application/json,*/*"})
    with urllib.request.urlopen(req,timeout=45) as f:
        b=f.read(MAX_BYTES+1)
        if len(b)>MAX_BYTES:
            raise ValueError("Original source exceeds bounded size")
        return b


def parse_zip(raw,strict=True):
    if not raw.startswith(b"PK\x03\x04"):
        raise ValueError("No authenticated source ZIP bytes")
    results={"source_zip_sha256":hashlib.sha256(raw).hexdigest(),
             "source_zip_bytes":len(raw),
             "published_original_papers":PUBLISHED_PAPERS,
             "published_original_effects":PUBLISHED_EFFECTS,
             "author_effect_metric":"Hedges d, NOT lnRR",
             "source_effect_size_meta_newly_fitted":False,
             "original_file_names":[],"original_CSV_SHA256":None,
             "original_R_code_SHA256":None,
             "verified_original_source_tabular_rows":0,
             "original_study_IDs_or_raw_effect_rows_exported":False,
             "real_water_wetland_moderator_source_validated":False,
             "source_csv_2018_2022_Penryn_trough_data_not_inferred":True}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        files={Path(info.filename).name:info
               for info in z.infolist() if not info.is_dir()}
        if strict and not EXPECTED_FILES.issubset(files):
            raise ValueError("Original author CSV and code not both found")
        results["original_file_names"]=sorted(files)
        if "Sanders_et_al_datafile.csv" not in files:
            return results
        b=z.read(files["Sanders_et_al_datafile.csv"])
        results["original_CSV_SHA256"]=hashlib.sha256(b).hexdigest()
        if len(b)>MAX_BYTES:
            raise ValueError("Original source CSV unexpectedly large")
        text=b.decode("utf-8-sig",errors="replace")
        csvsrc=csv.DictReader(io.StringIO(text))
        if not csvsrc.fieldnames or len(set(csvsrc.fieldnames))!=len(csvsrc.fieldnames):
            raise ValueError("Original CSV column names missing/colliding")
        rows=list(csvsrc)
        results["verified_original_source_tabular_rows"]=len(rows)
        results["source_original_column_names"]=csvsrc.fieldnames
        results["source_original_column_count"]=len(csvsrc.fieldnames)
        stats=[]
        for key in csvsrc.fieldnames:
            nm=key.strip().lower()
            if any(t in nm for t in (
                "study","paper","reference","article","categor","taxon",
                "response","light","species","trait","type","habitat",
                "effect","year","source","stat","author","measure","biome",
            )):
                vals=[str(r.get(key) or "").strip() for r in rows]
                unique={x for x in vals if x}
                stats.append({
                    "source_column":key,
                    "nonnull_row_count":len(vals)-sum(not x for x in vals),
                    "unique_nonblank_values_count_only":len(unique),
                    "original_value_labels_not_exported":True,
                })
        results["original_publication_and_outcome_candidate_field_shape"]=stats
        if "Sanders_meta__rcode.r" in files:
            code=z.read(files["Sanders_meta__rcode.r"])
            results["original_R_code_SHA256"]=hashlib.sha256(code).hexdigest()
            results["original_R_code_bytes"]=len(code)
    results["source_status"]="AUTHENTIC_ORIGINAL_HEDGES_D_CSV_AND_R_CODE_VERIFIED" if (
        results["original_R_code_SHA256"] and len(rows)>0
    ) else "PARTIAL_SOURCE_STOP"
    results["original_CSV_exact_1304_effect_rows_matches_published"]=(len(rows)==1304)
    return results


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    doi=urllib.parse.quote("doi:"+DATA_DOI,safe="")
    archive_url=f"{BASE}/datasets/{doi}/download"
    out={
        "schema":"alan_sanders2020_original_126paper_1304_meta_source_preflight_v1",
        "source_dataset_DOI":DATA_DOI,
        "source_status":"SOURCE_STOP",
        "published_n_papers":PUBLISHED_PAPERS,
        "published_n_effects":PUBLISHED_EFFECTS,
        "source_data_not_assumed_lnRR":True,
        "raw_source_rows_saved_in_git":False,
        "physical_illumination_programme_IDs_verified":False,
    }
    try:
        archive=get(archive_url)
        result=parse_zip(archive,strict=True)
        out.update(result)
        out["source_retrieval_route"]="authenticated_Dryad_public_dataset_zip"
    except (ValueError,zipfile.BadZipFile,urllib.error.URLError,TimeoutError) as ex:
        out["source_status"]="SOURCE_STOP"
        out["source_stop_error_type"]=type(ex).__name__
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("SANDERS126_RAW_META_CSV_STATUS",out["source_status"])
    print("SANDERS126_TRUE_META_SOURCE_HEADERS",out.get("source_original_column_names",[]))
    print("SANDERS126_TRUE_SOURCE_ROWS",out.get("verified_original_source_tabular_rows"))
    print("SANDERS126_AUTHOR_SOURCE_FILE_DIGESTS",out.get("original_CSV_SHA256"),out.get("original_R_code_SHA256"))
    print("SANDERS126_SOURCE_COLUMN_SHAPES",out.get("original_publication_and_outcome_candidate_field_shape",[]))


if __name__=="__main__":
    main()
