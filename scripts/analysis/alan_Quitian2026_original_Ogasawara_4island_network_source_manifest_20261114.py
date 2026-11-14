#!/usr/bin/env python3
"""Recover legitimate original four Ogasawara island plant–pollinator network workbook.

Quitian et al. 2026, Conservation Biology, DOI 10.1111/cobi.70304.
Original SOURCE Zenodo DOI 10.5281/zenodo.19221853 published 2026-03-25:
Data_Quitianetal2026_B.xlsx MD5=68f295b91605ea82bce0c05bffc7f884.

The published paper reports 4,174 interactions, 55 plant spp and 115 source
pollinator species over four islands and 2021 May/Jul/Sep. Legitimate visits
involve touching floral reproductive parts but pollen efficacy/fruit set is
not assessed for every guild. The author studies invasive anoles + honey bees,
but NOT an ALAN light experiment.

Source-manifest ONLY step: inspect 129.5KB authentic raw workbook schema and
categorical levels. No protected original interaction rows or locations
exported. No fabricated native guild label or pollutant effect.
"""
from __future__ import annotations
import argparse,hashlib,io,json,urllib.request,zipfile
from pathlib import Path
from openpyxl import load_workbook
URL="https://zenodo.org/records/19221853/files/Data_Quitianetal2026_B.xlsx?download=1"
SOURCE_MD5="68f295b91605ea82bce0c05bffc7f884"
DOI="10.5281/zenodo.19221853"
PAPER="10.1111/cobi.70304"
def source():
    req=urllib.request.Request(URL,headers={
       "User-Agent":"ALAN-independent-source-island-floral-network-audit/0.1",
       "Accept":"application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=40) as f:
        b=f.read(1_000_000)
    if hashlib.md5(b).hexdigest()!=SOURCE_MD5:
        raise ValueError("Zenodo original author data MD5 source changed")
    if not zipfile.is_zipfile(io.BytesIO(b)):
        raise ValueError("Original author file not original XLSX zip")
    return b

def summary(out):
    out=Path(out)
    if out.exists():raise FileExistsError("Refuse overwrite")
    b=source()
    w=load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    sheets=[]
    for s in w.worksheets:
        actual=[]
        for row in s.iter_rows(min_row=1,max_row=min(3,s.max_row),values_only=True):
            actual.append([str(v).strip()[:150] if v is not None else ""
                           for v in row[:80]])
        hdr=actual[0] if actual else []
        categorical={}
        names={k.lower().strip():i for i,k in enumerate(hdr) if k}
        for col in ("island","season","anole","invasional context","forest status","pollinator order","pollinator class","pollinator family","habitat"):
            if col in names:
                idx=names[col]
                values=set()
                for row in s.iter_rows(min_row=2,max_row=min(5000,s.max_row),
                                       min_col=idx+1,max_col=idx+1,values_only=True):
                    x=row[0]
                    if x is not None:values.add(str(x).strip())
                categorical[col]={"levels":sorted(values)[:55],
                                   "number_unique":len(values)}
        sheets.append({"sheet":s.title,"nrows":s.max_row,"ncols":s.max_column,
             "column_names":hdr,"categorical_source_levels":categorical,
             "individual_interaction_rows_not_exported":True})
    w.close()
    status={
        "schema":"ALAN_Ogasawara2021_four_islands_original_legitimate_pollination_source_manifest_v1",
        "status":"ZENODO_AUTHENTIC_4ISLAND_RAW_XLSX_SCHEMA_VERIFIED_NO_NEW_FIT",
        "paper_DOI":PAPER,"original_dataset_DOI":DOI,"file":URL,
        "publisher_author_MD5_verified":SOURCE_MD5,
        "source_file_SHA256":hashlib.sha256(b).hexdigest(),
        "source_file_bytes":len(b),
        "sheets":sheets,
        "paper_reported_not_refitted":{
          "archipelago":"Ogasawara",
          "source_islands":["Chichijima","Hahajima","Anijima","Ototojima"],
          "original_survey_year":2021,
          "months":["May","July","September"],
          "reported_legitimate_plant_visitor_interactions":4174,
          "reported_pollinator_taxa":115,
          "reported_flowering_plant_species":55,
          "anole_predator_nuisance_reported":True,
          "honeybee_invasive_floral_competitor_reported":True,
          "source_field_flower_abundance_estimated":True
        },
        "no_fabricated_claims":{
          "no_ALAN_lux_or_light_treatment_in_this_island_study":True,
          "anole_presence_not_numerical_predation_rate":True,
          "reproductive_structure_touch_not_measured_pollen_deposition":True,
          "plant_visitor_phenology_month_only_not_first_flowering_date":True,
          "no_natural_pollinator_shift_baseline_inferred_from_one_year":True,
          "no_raw_original_species_location_site_records_exported":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(status,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_OGASAWARA_2021_REAL_AUTHOR_ISLAND_NETWORK_SCHEMA",
          json.dumps(status,ensure_ascii=False),flush=True)
    return status

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();summary(a.out)
