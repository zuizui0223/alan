#!/usr/bin/env python3
"""Independent public original-workbook replay for 2026 Xue 6 ALAN light zones.

The main research repository earlier extracted source 77 cylinder records
into 24 zone×species biomass means. Here verify the ORIGINAL ScienceDB
Monoculture_Growth.xlsx in an independently executed new CI workflow:
original file bytes digest, row counts, treatment zoning, surviving plant
IDs, 24 original treatment-zone-species cell means and counts; only
the audit receipt (no original individual biomass records) is emitted.

If the Science DataBank origin cannot be re-downloaded, source status is
explicitly STOP and the existing Git-pinned 24-cell analysis remains
reproducible on its own. Do NOT infer raw validation from successful tests.
"""
from __future__ import annotations
import argparse, hashlib, io, json, urllib.request, urllib.error
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_xue2026_original6zone_native_invader_divergence import (
    read_source, validate_zones, ORIGINAL_SHA
)

SCHEMA="alan_Xue2026_original_source_second_independent_ScienceDB_workbook_replay_v1"
URL="https://www.scidb.cn/api/sdb-download-service/downloadFileMole?fileId=721cbfeb236f2d2fc05bdaa21f9b23d6"
MAX=3<<20

def execute():
    report={
      "schema":SCHEMA,
      "official_original_dataset":"10.57760/sciencedb.35715",
      "expected_original_XLSX_SHA256":ORIGINAL_SHA,
      "original_source_url":URL,
      "original_extracted_24_cells_git_verified":True,
      "author_original_workbook_SHA_verified_this_run":False,
      "original_surviving_cylinder_records_validated_this_run":False,
      "source_numeric_zone_means_identical_to_original_this_run":False,
      "source_independent_light_zones":6,
      "no_new_plant_effect_fitted_from_unverified_bytes":True,
      "status":"ORIGINAL_WORKBOOK_SECOND_FETCH_STOP",
    }
    src=read_source()
    fixed=src["zone"].copy()
    validate_zones(fixed)
    try:
        req=urllib.request.Request(URL,headers={"User-Agent":"ALAN-Xue-independent-original-data-verify/2026",
                        "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"})
        with urllib.request.urlopen(req,timeout=90) as response:
            original=response.read(MAX+1)
        if len(original)>MAX:
            raise RuntimeError("SOURCE_XLSX_EXCEEDS_3MB_STOP")
        actual=hashlib.sha256(original).hexdigest()
        report["retrieved_original_file_bytes"]=len(original)
        report["observed_SHA256"]=actual
        if actual!=ORIGINAL_SHA:
            report["status"]="ORIGINAL_XLSX_SOURCE_SHA_DRIFT_STOP"
            return report
        report["author_original_workbook_SHA_verified_this_run"]=True
        d=pd.read_excel(io.BytesIO(original),engine="openpyxl")
        reqcols={"Light","Zone","Species","ID","Rand","Total_biomass"}
        if not reqcols.issubset(d):
            raise ValueError("ORIGINAL_XUE_REQUIRED_COLUMNS_CHANGED")
        if len(d)!=77 or d.ID.nunique()!=77 or d.Rand.nunique()!=77:
            raise ValueError("ORIGINAL_XUE_77_UNIQUE_SURVIVORS_CHANGED")
        if d.Zone.nunique()!=6 or d.groupby("Zone").Light.nunique().ne(1).any():
            raise ValueError("ORIGINAL_XUE_LIGHT_INSTALLATION_DRIFT")
        d["Species"]=d.Species.str.strip()
        grouped=d.groupby(["Light","Zone","Species"],as_index=False).agg(
          biomass_g=("Total_biomass","mean"),
          n_surviving_cylinders=("ID","size"))
        cols=["Light","Zone","Species","biomass_g","n_surviving_cylinders"]
        given=fixed[cols].sort_values(["Light","Zone","Species"]).reset_index(drop=True)
        now=grouped[cols].sort_values(["Light","Zone","Species"]).reset_index(drop=True)
        if len(now)!=24 or not given[["Light","Zone","Species"]].equals(
                                    now[["Light","Zone","Species"]]):
            raise ValueError("ORIGINAL_XUE_SOURCE_ZONE_GROUP_IDS_DRIFT")
        if not np.allclose(given.biomass_g,now.biomass_g,rtol=1e-11,atol=1e-11) or not (
            given.n_surviving_cylinders==now.n_surviving_cylinders).all():
            raise ValueError("ORIGINAL_XUE_SOURCE_ZONE_BIOMASS_MEAN_DRIFT")
        report.update({
          "status":"ORIGINAL_WORKBOOK_77_CYLINDERS_24_ZONE_CELLS_IDENTICAL",
          "original_surviving_cylinder_records_validated_this_run":True,
          "source_numeric_zone_means_identical_to_original_this_run":True,
          "source_original_77_records":77,
          "source_original_24_zone_plant_cells":24,
          "source_original_3_ALAN_and_3_dark_zones_confirmed":True,
        })
    except Exception as e:
        report["status"]="ORIGINAL_WORKBOOK_SECOND_FETCH_STOP" if not report[
          "author_original_workbook_SHA_verified_this_run"] else "SOURCE_AFTER_SHA_SCHEMA_STOP"
        report["error_class"]=type(e).__name__
        report["http_code"]=getattr(e,"code",None)
        report["no_original_bytes_or_uncorroborated_rows_committed"]=True
    return report

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    result=execute()
    path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("XUE2026_INDEPENDENT_ORIGINAL_XLSX_SHA_AND_ZONE_REPLAY "+
          json.dumps(result,sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
