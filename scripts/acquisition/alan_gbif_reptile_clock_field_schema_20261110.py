#!/usr/bin/env python3
"""Audit whether source-authenticated Squamata GBIF Parquet has clock metadata.

Only source footers are inspected, not individual animal observations.
The existence of an eventtime column does NOT establish gecko activity
timing or reliable local time. No field-level pollination data.
"""
from __future__ import annotations
import argparse,datetime,json
from pathlib import Path
import pyarrow.parquet as pq
from pyarrow.fs import S3FileSystem

DOI="10.15468/dl.8yrbe7"
GROUP="taxonpartition=Animalia_Chordata_Reptilia_Squamata"
SCHEMA="ALAN_GBIF_202608_Squamata_time_of_day_public_schema_source_audit_v1"
CHECK=("eventtime","eventdate","verbatimeventdate","day","month","year",
       "timeofday","localtime","timezone","verbatimeventtime",
       "samplingprotocol","samplesizevalue","individualcount",
       "datasetkey","basisofrecord","recordedby","gbifid")
def require(ok,msg):
    if not ok:raise ValueError("ALAN_GBIF_DIEL: "+str(msg))

def run(manifest,out):
    out=Path(out)
    require(not out.exists(),"Cannot overwrite original public source schema")
    m=json.loads(Path(manifest).read_text())
    require(m["schema"]=="ALAN_original_GBIF_202608_Gekkota_reptile_Parquet_partition_S3_audit_v1"
            and m["original_GBIF_snapshot_DOI"]==DOI,
            "Not original source-verified Squamata Parquet")
    files=m["matching_reptile_gecko_partition_folders"][GROUP]["original_source_object_metadata"]
    fs=S3FileSystem(anonymous=True,region="us-east-1")
    info=[]
    for v in files:
        f=fs.open_input_file("gbif-public-data/"+v["key"])
        p=pq.ParquetFile(f)
        columns=list(p.schema_arrow.names)
        lower={x.casefold():x for x in columns}
        found={name:lower.get(name) for name in CHECK}
        info.append({
          "original_source_object":v["key"],
          "actual_original_rows":p.metadata.num_rows,
          "columns_available":len(columns),
          "requested_source_clock_and_effort_fields":found})
        f.close()
    flags={name:all(j["requested_source_clock_and_effort_fields"][name]
                    for j in info) for name in CHECK}
    result={
       "schema":SCHEMA,"true_original_GBIF_DOI":DOI,
       "true_original_reptile_source_rows":
           sum(x["actual_original_rows"] for x in info),
       "real_original_file_field_schemas":info,
       "time_or_month_columns_in_all_original_source_files":flags,
       "no_source_observation_clock_values_read_or_inferred":True,
       "time_of_day_field_does_not_imply_known_timezone":True,
       "no_gecko_moth_predation_or_activity_observed":True,
       "source_time_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_GECKO_GBIF_GENUINE_TIME_FIELDS",
          json.dumps({"source_rows":result["true_original_reptile_source_rows"],
            "fields":flags},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--manifest",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    x=a.parse_args()
    run(x.manifest,x.out)
