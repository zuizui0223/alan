#!/usr/bin/env python3
"""Source-locked schema-only audit of the EERES2021 drift taxon columns.

No source biological observations or individual unit rows are exported.
This checks whether the German EERES source has taxon resolution to compare
its in-stream drift *composition* with the distinct Henn2014 five-stream
experiment (different light conditions and field-vs-flume design).

Source DOI: 10.5061/dryad.rfj6q57k5. Physical experiment 6 header tanks.
"""
import argparse
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
GATE=ROOT/"scripts/acquisition/alan_eeres_2021_source_gate.py"
SHA="496dcfa07436619084dc677317525dc674775e9965a0fdfb482cfa1cb4d8977f"
spec=importlib.util.spec_from_file_location("source_gate_eeres",GATE)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def check_schema(source):
    if source.get("status")!="SOURCE_WORKBOOK_VERIFIED":
        raise ValueError("STOP: original source unavailable")
    workbook=source.get("original_workbook",{})
    if workbook.get("sha256")!=SHA:
        raise ValueError("STOP: original EERES workbook hash changed")
    meta=workbook.get("metadata",{})
    if meta.get("status")!="SOURCE_EXCEL_SHEETS_READ":
        raise ValueError("STOP: no original source metadata")
    sheets=meta.get("sheets",{})
    out={"physical_programme":"EERES2021_Disconnected18flumes_6headers",
         "dataset_doi":"10.5061/dryad.rfj6q57k5",
         "verified_sha256":SHA,
         "status":"SOURCE_COLUMNS_ONLY_NOT_TAXON_GUILD_INFERENCE",
         "source_observation_rows_shared":0,
         "no_meta_effect_or_pooled_contrast":True,
         "sheets":{}}
    for sheet, offset, expected_n in (("drift",8,144),("benthos",6,108)):
        block=sheets.get(sheet)
        if not block or block.get("n_nonempty_data_rows")!=expected_n:
            raise ValueError(f"STOP: EERES {sheet} count changed")
        fields=block.get("columns")
        if not isinstance(fields,list) or len(fields)<=offset:
            raise ValueError("STOP: no taxonomic columns")
        out["sheets"][sheet]={
            "n_original_rows":expected_n,
            "metadata_column_count":offset,
            "biological_measurement_column_count":len(fields)-offset,
            "biological_column_labels":[str(x) for x in fields[offset:]],
            "taxonomy_classification_status":"UNCLASSIFIED_SOURCE_LABELS_DO_NOT_ASSUME_ORDER",
        }
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=Path("results/eeres2021_taxonomy_schema_20261124.json"))
    args=ap.parse_args()
    obj=check_schema(mod.run())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("SCHEMA_ONLY",obj["physical_programme"])
    for name,v in obj["sheets"].items():
        print("ORIGINAL_COLUMNS",name,"N",v["biological_measurement_column_count"],
              "NAMES",json.dumps(v["biological_column_labels"],ensure_ascii=False))


if __name__=="__main__":
    main()
