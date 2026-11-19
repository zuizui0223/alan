#!/usr/bin/env python3
"""Bennie et al. 2015 primary original grassland aphid-light Excel access gate.

The study has ALAN spectral treatments, pea aphid, floral resource, predator
exclusion. Only original workbook SHA, sheet/column names, record counts,
and *number* of physical light plots are safe output. Do not infer from
aphid or flower-head repeated rows that they are independently lit plots.
"""
import argparse
import hashlib
import io
import json
import urllib.error
import urllib.request
from pathlib import Path

ARTICLE_DOI="10.1098/rstb.2014.0131"
PMCID="PMC4375371"
EXPECTED_FILE="rstb20140131supp1.xlsx"
SOURCE_LOCATORS=[
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC4375371/bin/rstb20140131supp1.xlsx",
  "https://pmc.ncbi.nlm.nih.gov/articles/instance/4375371/bin/rstb20140131supp1.xlsx",
  "https://royalsocietypublishing.org/doi/suppl/10.1098/rstb.2014.0131/suppl_file/rstb20140131supp1.xlsx"
]
MAX_BYTES=1_500_000

def get(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-original-grassland-2015-source-gate/1.0 (academic reproducibility)",
        "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"})
    with urllib.request.urlopen(req,timeout=28) as f:
        b=f.read(MAX_BYTES+1)
        if len(b)>MAX_BYTES:raise ValueError("Original source exceeds bounded preflight")
        return b

def workbook_schema(source):
    from openpyxl import load_workbook
    if not source.startswith(b"PK\x03\x04"):
        raise ValueError("Not valid XLSX file signature")
    wb=load_workbook(io.BytesIO(source),read_only=True,data_only=True)
    tabs=[]
    for ws in wb.worksheets:
        first_nonempty=[]
        n_nonempty=0
        for idx,row in enumerate(ws.iter_rows(values_only=True),start=1):
            z=[str(x).strip() if x is not None else "" for x in row]
            if not any(z):continue
            n_nonempty+=1
            if len(first_nonempty)<6:
                first_nonempty.append({"original_row":idx,
                    "n_nonempty_cells":sum(bool(x) for x in z),
                    # Header candidates only; source values are intentionally
                    # NOT included if row number greater than first six.
                    "column_count":len(z)})
        # Never publish source individual plot-treatment records or biology.
        tabs.append({"sheet":ws.title,"data_nonempty_rows":n_nonempty,
                     "first_rows_shape_only":first_nonempty,
                     "max_cols":ws.max_column,
                     "original_source_plot_IDs_not_exported":True})
    wb.close()
    return tabs

def audit():
    out={
        "schema":"alan_bennie2015_original_flora_aphid_supplement_source_gate_v1",
        "original_primary_publication_DOI":ARTICLE_DOI,
        "original_PMC_ID":PMCID,
        "original_publication_suppl_xlsx":EXPECTED_FILE,
        "source_original_design_hint":"Grassland spectral lighting with aphid and predation factors, original lighting installation IDs NOT inferred",
        "source_status":"SOURCE_STOP",
        "attempts":[],
        "published_one_physical_lamp_assumed_per_plant_or_aphid":False,
        "original_figshare_review_DOI_mismatch_if_any_is_not_a_true_study":True,
        "no_plant_reproductive_fitness_measured":True,
        "actual_lighting_effect_fit_from_original_samples":False,
        "no_original_study_rows_or_site_IDs_exported":True
    }
    for url in SOURCE_LOCATORS:
        try:
            b=get(url)
            original=workbook_schema(b)
            if not original:raise ValueError("No original sheets")
            out["source_status"]="RAW_SOURCE_XLSX_SCHEMA_VERIFIED"
            out["source_SHA256"]=hashlib.sha256(b).hexdigest()
            out["source_bytes"]=len(b)
            out["source_original_sheet_structure"]=original
            out["source_success_URL"]=url
            out["attempts"].append({"URL":url,"status":"SOURCE_RETRIEVED","size":len(b)})
            break
        except (ValueError,ImportError,urllib.error.URLError,urllib.error.HTTPError,TimeoutError) as e:
            out["attempts"].append({"URL":url,"status":"SOURCE_STOP",
                                    "error_type":type(e).__name__,
                                    **({"http_status":e.code} if isinstance(e,urllib.error.HTTPError) else {})})
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    x=audit()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_BENNIE2015_ORIGINAL_PRIMARY_SOURCE",x["source_status"])
    print("ALAN_BENNIE2015_PRIMARY_ORIGINAL_SOURCE_ATTEMPTS",x["attempts"])
    print("ALAN_BENNIE2015_PRIMARY_WORKBOOK_SCHEMAS",x.get("source_original_sheet_structure",[]))

if __name__=="__main__":
    main()
