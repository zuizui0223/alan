#!/usr/bin/env python3
"""Audit 2026 GBIF old Backbone vs COL XR gecko taxon KEY compatibility.

This is a diagnostic gate: NEVER equate a 0 result from an incompatible
taxonomy key with gecko absence. Official docs 2026 taxon backend migration.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,time,urllib.parse,urllib.request
from pathlib import Path

LEGACY="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
COL="7ddf754f-d193-4cc9-b351-99906754a03b"
API="https://api.gbif.org/v1"
DOI="10.15468/dl.8yrbe7"
SCHEMA="ALAN_true_GBIF_gecko_dual_taxonomy_2026_old_COLXR_key_validation_v1"
NAMES=(("Gekkonidae","FAMILY"),("Hemidactylus","GENUS"),
       ("Gekko","GENUS"),("Phelsuma","GENUS"),
       ("Hemidactylus frenatus","SPECIES"),
       ("Gekko japonicus","SPECIES"))

def request(endpoint,args):
    url=API+endpoint+"?"+urllib.parse.urlencode(args,doseq=True)
    r=urllib.request.Request(url,headers={"Accept":"application/json",
        "User-Agent":"ALAN-gecko-taxonomy-source-validation/1.0"})
    for i in range(4):
        try:
            with urllib.request.urlopen(r,timeout=60) as z:
                data=z.read()
                if z.status!=200:raise RuntimeError("Unexpected GBIF HTTP")
            return json.loads(data),url,hashlib.sha256(data).hexdigest()
        except Exception:
            if i==3:raise
            time.sleep(i+1)
    raise AssertionError("No GBIF API")

def one_count(label,params):
    d,u,h=request("/occurrence/search",[("limit","0"),*params])
    if not isinstance(d.get("count"),int):raise ValueError("No real GBIF count")
    return {"label":label,"actual_GBIF_record_count":d["count"],
            "source_query":u,"source_response_sha256":h}

def run(out):
    out=Path(out)
    if out.exists():raise ValueError("Refuse overwrite original 2026 gecko audit")
    matches={};gecko=None;comparison={}
    for name,rank in NAMES:
        matches[name]={}
        for key in (LEGACY,COL):
            d,u,h=request("/species/match",{
                "name":name,"strict":"true","checklistKey":key})
            matches[name][key]={"usageKey":d.get("usageKey"),
                 "rank":d.get("rank"),"matchType":d.get("matchType"),
                 "source_url":u,"sha256":h}
        print("ALAN_GECKO_GBIF_DUAL_TAXON_MATCH",
              json.dumps({name:matches[name]},sort_keys=True),flush=True)
    base=[
        ("legacy_key_no_declared_checklist",[("taxonKey","5666")]),
        ("legacy_key_explicit_legacy_checklist",[
            ("taxonKey","5666"),("checklistKey",LEGACY)]),
        ("COL_correct_matching_gecko_key",[("taxonKey",
             str(matches["Gekkonidae"][COL].get("usageKey"))),
             ("checklistKey",COL)]),
        ("gecko_scientificName_legacy_checklist",[
             ("scientificName","Gekkonidae"),("checklistKey",LEGACY)]),
        ("gecko_scientificName_COL_checklist",[
             ("scientificName","Gekkonidae"),("checklistKey",COL)])]
    for label,filters in base:
        try:
            comparison[label]={
                "all_status_all_basis":one_count(
                    label+"_no_other_filters",filters),
                "2013_2025_human_georef_PRESENT":one_count(
                    label+"_source_2013_2025_human",
                    [*filters,
                     ("basisOfRecord","HUMAN_OBSERVATION"),
                     ("hasCoordinate","true"),
                     ("year","2013,2025"),("occurrenceStatus","PRESENT")])
            }
        except Exception as e:
            comparison[label]={"failed":str(e)[:350]}
        print("ALAN_GEK_XR_TRUE_COUNT_COMPARE",
              json.dumps({label:comparison[label]},sort_keys=True),flush=True)
    # Species controls independent of the family taxon key. Explicit
    # GEKKO JAP and invasive HEMIDACTYLUS; Phelsuma diurnal contrast.
    focus={}
    for name in ("Hemidactylus frenatus","Gekko japonicus","Phelsuma"):
        xr=matches[name][COL].get("usageKey")
        if xr:
            try:
                focus[name]=one_count(name+"_COL",
                    [("taxonKey",str(xr)),("checklistKey",COL),
                     ("year","2013,2025"),("basisOfRecord","HUMAN_OBSERVATION"),
                     ("hasCoordinate","true"),("occurrenceStatus","PRESENT")])
            except Exception as e:
                focus[name]={"source_error":str(e)[:200]}
    out.parent.mkdir(parents=True,exist_ok=True)
    data={"schema":SCHEMA,"original_GBIF_2026_preview_DOI":DOI,
          "public_GBIF_taxonomy_guidance":
            "https://techdocs.gbif.org/en/data-processing/taxonomy-interpretation",
          "dual_taxonomy_original_exact_name_match":matches,
          "Gekkonidae_source_count_query_variants":comparison,
          "curated_invasive_native_day_gecko_source_counts":focus,
          "gecko_no_records_not_true_absence":True,
          "legacy_taxon_key_zero_not_real_gecko_zero":True,
          "source_taxonomic_backbone_keys_are_not_safe_to_mix":True,
          "no_actual_moth_gecko_interaction_yet":True,
          "UTC_timestamp":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    return data

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.out)
