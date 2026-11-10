#!/usr/bin/env python3
"""Identify GBIF reptile source filter that gave an artificial zero.

Taxonomic Gekkonidae source family has ~0.6m actual records; **never**
interpret 0 for an over-restrictive search as biological gecko absence.
"""
from __future__ import annotations
import argparse,datetime,json
from pathlib import Path
from scripts.acquisition.alan_GBIF_gecko_taxonomic_backbone_source_diagnostic_20261110 import (
    one_count,LEGACY,COL)

SCHEMA="ALAN_GBIF_gecko_family_2026_occurrence_filter_intersection_audit_v1"
ROOT=[("taxonKey","5666"),("checklistKey",LEGACY)]
OPTIONS={
 "year":[("year","2013,2025")],
 "basis_HUMAN":[("basisOfRecord","HUMAN_OBSERVATION")],
 "hasCoordinates":[("hasCoordinate","true")],
 "PRESENT":[("occurrenceStatus","PRESENT")],
 "DETECTED":[("occurrenceStatus","DETECTED")],
}

def run(out):
    out=Path(out)
    if out.exists():raise ValueError("Refuse overwrite source filter audit")
    combos=[
      (),("year",),("basis_HUMAN",),("hasCoordinates",),("PRESENT",),
      ("DETECTED",),("year","basis_HUMAN"),("year","hasCoordinates"),
      ("basis_HUMAN","hasCoordinates"),
      ("year","basis_HUMAN","hasCoordinates"),
      ("year","basis_HUMAN","PRESENT"),
      ("year","basis_HUMAN","DETECTED"),
      ("year","hasCoordinates","PRESENT"),
      ("basis_HUMAN","hasCoordinates","PRESENT"),
      ("year","basis_HUMAN","hasCoordinates","PRESENT"),
      ("year","basis_HUMAN","hasCoordinates","DETECTED"),
    ]
    rows={}
    for combo in combos:
        params=list(ROOT)
        for name in combo:params.extend(OPTIONS[name])
        label="+".join(combo) or "unrestricted"
        rows[label]=one_count(label,params)
        print("ALAN_GEK_GEKKONID_TRUE_FILTER_SOURCE",
              json.dumps({label:rows[label]["actual_GBIF_record_count"]},
                         sort_keys=True),flush=True)
    summary={
      "schema":SCHEMA,"correct_original_family_GBIF_backbone_key":5666,
      "legacy_checklistKey":LEGACY,
      "true_original_Gekkonidae_source_filter_counts":rows,
      "no_gecko_absence_inference_from_gbif_filtered_zeros":True,
      "no_actual_predation_or_pollination_causal_effect":True,
      "source_timestamp_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.out)
