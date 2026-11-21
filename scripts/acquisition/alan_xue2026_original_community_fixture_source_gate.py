#!/usr/bin/env python3
"""Xue 2026 actual mixed-community original ScienceDB workbook fixture audit.

Original DOI 10.57760/sciencedb.35715 Community.xlsx SHA256
24b6cd24c40403803cdf5fa9d9251c46383a585bbd8da168d8917a81be2644ad.
Original 57 measured mesocosm containers / 60 assigned; NO assumption
of 60 physical independent lighting fixtures. Original source has 3
missing biomass containers all ALAN+native snail (not zero).

Goal: independently download author original, verify original 57 rows and
all 12 published source group means, identify source documented light
zone/fixture keys, and STOP causal light×invasive-snail estimation unless
physical illumination allocation units can be certified. All original
individual plant biomass rows remain uncommitted; only source schemas,
counts, SHA and derived cell statistics appear in receipt.
"""
from __future__ import annotations
import argparse,hashlib,io,json,urllib.request
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA="alan_Xue2026_original_community_57containers_fixture_source_gate_v1"
URL="https://www.scidb.cn/api/sdb-download-service/downloadFileMole?fileId=c10ee136b42f743bcf291e328ed840ed"
ORIGINAL_SHA="24b6cd24c40403803cdf5fa9d9251c46383a585bbd8da168d8917a81be2644ad"
ROOT=Path(__file__).resolve().parents[2]
COMMUNITY=ROOT/"database/alan_experiments/xue_sources/community_descriptive_only.csv"
PIN="8ebf5148126ca2610532f1d2cae2a8f51e505165"
def chk(ok,why):
    if not ok: raise RuntimeError("XUE_ORIGINAL_57_COMMUNITY_SOURCE_STOP:"+why)
def git_blob(raw):
    return hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()

def execute():
    report={
      "schema":SCHEMA,
      "original_sciencedb_doi":"10.57760/sciencedb.35715",
      "original_community_xlsx_SHA256_expected":ORIGINAL_SHA,
      "author_original_xlsx_source_url":URL,
      "original_community_57_source_rows_verified":False,
      "original_community_2x3_group_means_verified":False,
      "true_independent_lighting_fixture_units_VERIFIED":False,
      "source_community_57_containers_not_57_independent_lights":True,
      "no_original_raw_container_measurements_committed":True,
      "no_new_inferential_light_x_invasive_snail_coefficient_fitted":True,
      "status":"SECONDARY_ORIGINAL_SOURCE_ACCESS_STOP",
    }
    data=COMMUNITY.read_bytes()
    chk(git_blob(data)==PIN,"existing 12-source-cell source pin changed")
    fixed=pd.read_csv(io.BytesIO(data))
    chk(len(fixed)==12,"source original community group cells changed")
    try:
        req=urllib.request.Request(URL,headers={
           "User-Agent":"ALAN-Xue2026-independent-community-source-audit/2026",
           "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"})
        with urllib.request.urlopen(req,timeout=75) as res:
            x=res.read(3<<20)
        digest=hashlib.sha256(x).hexdigest()
        report["retrieved_bytes"]=len(x)
        report["observed_original_SHA256"]=digest
        if digest!=ORIGINAL_SHA:
            report["status"]="SECONDARY_ORIGINAL_SOURCE_HASH_DRIFT_STOP"
            return report
        z=pd.read_excel(io.BytesIO(x),engine="openpyxl")
        required={"ID","Rand","Light","Herbivore",
                  "Native_total_biomass","Invasive_total_biomass",
                  "Invasive_proportional_biomass"}
        chk(required.issubset(z.columns),"source expected original native/invasive measures absent")
        chk(len(z)==57 and z.ID.nunique()==57 and z.Rand.nunique()==57,
            "source 57 survivor independent containers changed")
        chk(np.isfinite(z[["Native_total_biomass","Invasive_total_biomass"]]).all().all(),
            "source native/invasive biomass has nonfinite observations")
        chk((z[["Native_total_biomass","Invasive_total_biomass"]]>=0).all().all(),
            "source original native/invasive biomass negative")
        calc=z.Invasive_total_biomass/(z.Native_total_biomass+z.Invasive_total_biomass)
        chk(np.allclose(calc,z.Invasive_proportional_biomass,rtol=1e-9,atol=1e-9),
            "original author source invasive proportion identity changed")
        rows=[]
        for (light,herb),g in z.groupby(["Light","Herbivore"]):
            for endpoint in ("Native_total_biomass","Invasive_total_biomass"):
                k=fixed[(fixed.light.eq(light))&(fixed.herbivore.eq(herb))&
                        (fixed.endpoint.eq(endpoint))]
                chk(len(k)==1,"original 12 source group mapping changed")
                chk(len(g)==int(k.iloc[0].observed_buckets),
                    "source group original observed containers inconsistent")
                m=float(g[endpoint].mean())
                chk(abs(m-float(k.iloc[0].mean_g))<1e-11,
                    "original source actual biomass mean changed")
                rows.append({"original_source_light":light,
                  "original_source_herbivore":herb,
                  "original_endpoint":endpoint,
                  "original_individual_mesocosms_observed":len(g),
                  "original_mean_g":m})
        # Only source-named explicit light fixture/zone fields could help
        # identify independent experimental light installations.
        rawcols=[str(c) for c in z.columns]
        candidates=[c for c in rawcols if any(term in c.lower() for term in
              ("zone","fixture","lamp","station","cluster","site","block","plot"))]
        report.update({
          "original_community_57_source_rows_verified":True,
          "original_community_2x3_group_means_verified":True,
          "original_community_columns":rawcols,
          "explicit_original_lighting_unit_candidate_column_names":candidates,
          "original_group_means_exact_to_pinned_repo":rows,
          "original_light_herbivore_observed_container_sizes":{
             f"{light} | {herb}":int(len(g))
             for (light,herb),g in z.groupby(["Light","Herbivore"])},
          "original_metadata_source_DOI_fixed":True,
          "source_group_missingness_not_zero_imputed":True,
          "status":"SOURCE_57_CONTAINERS_AND_12_MEANS_VERIFIED_LIGHT_FIXTURE_UNRESOLVED"
        })
        if "Zone" in z.columns:
            distinct=z.groupby("Zone").Light.nunique()
            report["original_explicit_Zones"]=int(z.Zone.nunique())
            report["original_zone_has_unique_light_arm"]=bool(distinct.eq(1).all())
            # Source still needs author design confirmation.
            report["source_fixture_inference_requiring_author_methods"]=True
        # Do not promote inferred light fixture count from original ID, Rand,
        # or only six group-treatment combinations.
    except Exception as e:
        report["status"]="SECONDARY_ORIGINAL_SOURCE_ACCESS_OR_SCHEMA_STOP"
        report["error_class"]=type(e).__name__
        report["HTTP_status"]=getattr(e,"code",None)
    return report

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    result=execute()
    path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("XUE2026_ORIGINAL_57_COMMUNITY_LAMP_UNIT_SOURCE_GATE "+
          json.dumps(result,sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
