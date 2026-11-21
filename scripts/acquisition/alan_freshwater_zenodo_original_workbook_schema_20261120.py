#!/usr/bin/env python3
"""Inspect source authentic Zenodo Campos et al 2026 worksheet grids.

Frozen original XLSX SHA256 b5c36d9b4c184a3f563388d56a9072fcfe9fa8c988600e9979054402083a7d6d
and original journal source doi 10.3389/fevo.2026.1830790.
Shows sheet headers and numeric support; no ecological effect estimated.
Uses stdlib XML reader from source-verified original Lu workbook workflow,
not Excel's formatted display or manually invented trial replicates.
"""
import argparse, collections, hashlib, io, json, urllib.request
from pathlib import Path
from scripts.acquisition.alan_lu2026_original_workbook_schema import sheet_data

SOURCE="https://zenodo.org/api/records/19022243/files/data_Campos%20et%20al%202026.xlsx/content"
SHA="b5c36d9b4c184a3f563388d56a9072fcfe9fa8c988600e9979054402083a7d6d"
SCHEMA="alan_Campos2026_Zenodo19022243_6sheet_original_structure_v1"
def run():
    with urllib.request.urlopen(urllib.request.Request(
            SOURCE,headers={"User-Agent":"ALAN-source-original-sheet-audit"}),timeout=60) as r:
        original=r.read(2<<20)
    if hashlib.sha256(original).hexdigest()!=SHA:
        raise RuntimeError("ORIGINAL_ZENODO_XLSX_SHA_CHANGED_STOP")
    sheets=sheet_data(original)
    if len(sheets)!=6:
        raise RuntimeError("ORIGINAL_ZENODO_EXPECTED_SIX_SHEETS_STOP")
    result={"schema":SCHEMA,"source_sha256":SHA,
            "source_spreadsheet_bytes":len(original),
            "original_worksheet_count":len(sheets),
            "source_article_doi":"10.3389/fevo.2026.1830790",
            "original_data_doi":"10.5281/zenodo.19022243",
            "status":"SHEET_SCHEMA_SOURCE_VERIFIED",
            "reproduced_ALAN_effects":False,
            "wet_plants_measurements":False,
            "worksheets":{}}
    for label,rows in sheets.items():
        sizes=collections.Counter(len(x) for x in rows)
        data=[r for r in rows if any(str(c).strip() for c in r)]
        first=data[:7]
        result["worksheets"][label]={
          "raw_nonempty_sheet_rows":len(data),
          "field_count_distribution":dict(sizes),
          "header_and_first_source_rows":[
             [str(c)[:90] for c in row[:22]] for row in first],
          "names_not_imputed_or_assumed_units":True,
          "no_collected_animal_treated_as_independent_light_unit":True,
        }
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    r=run()
    q=Path(a.out);q.parent.mkdir(parents=True,exist_ok=True)
    q.write_text(json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("CAMPO2026_6SHEET_OFFICIAL_SOURCE_STRUCTURE "+json.dumps(r,sort_keys=True),flush=True)
if __name__=="__main__":main()
