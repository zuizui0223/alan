#!/usr/bin/env python3
"""Verify authentic frozen GBIF 2026-08 Squamata Parquet footers and schema.

No source gecko presence inference until actual year-filtered records and
GSHHG island polygon membership are measured.
"""
from __future__ import annotations
import argparse,json,datetime
from pathlib import Path
import pyarrow.parquet as pq
from pyarrow.fs import S3FileSystem
import pandas as pd
import numpy as np

BUCKET="gbif-public-data"
SOURCE_GROUP="taxonpartition=Animalia_Chordata_Reptilia_Squamata"
SOURCE_DOI="10.15468/dl.8yrbe7"
SCHEMA="ALAN_original_GBIF_2026_08_Squamata_gecko_source_columns_v1"

def require(b,msg):
    if not b:raise ValueError("ALAN_GENUINE_GECKO_PARQUET: "+str(msg))

def get_col(cols,name):
    hits=[c for c in cols if c.casefold()==name.casefold()]
    require(len(hits)==1,"Original Squamata source missing "+name)
    return hits[0]

def run(manifest,out):
    out=Path(out);require(not out.exists(),"No overwrite source measurement")
    m=json.loads(Path(manifest).read_text())
    require(m["schema"]=="ALAN_original_GBIF_202608_Gekkota_reptile_Parquet_partition_S3_audit_v1"
            and m["original_GBIF_snapshot_DOI"]==SOURCE_DOI,
            "Original audited 2026-08 source reptile folder manifest absent")
    group=m["matching_reptile_gecko_partition_folders"].get(SOURCE_GROUP)
    require(group is not None and group["source_parquet_files"]==3,
            "Genuine source Squamata 3-blob partition unavailable")
    fs=S3FileSystem(anonymous=True,region="us-east-1")
    receipts=[];sample={}
    for i,v in enumerate(group["original_source_object_metadata"]):
        f=fs.open_input_file(BUCKET+"/"+v["key"])
        parquet=pq.ParquetFile(f)
        fields=parquet.schema_arrow.names
        names=("gbifid","year","family","genus","species","specieskey",
               "basisofrecord","occurrencestatus",
               "decimallatitude","decimallongitude","datasetkey")
        selected={k:get_col(fields,k) for k in names}
        metadata=parquet.metadata
        require(metadata.num_rows>0 and metadata.num_row_groups>0,
                "Original public reptile source Parquet empty")
        row={"original_object_key":v["key"],
             "original_compressed_source_bytes":v["compressed_original_bytes"],
             "original_parquet_true_rows":metadata.num_rows,
             "original_parquet_row_groups":metadata.num_row_groups,
             "source_field_mapping":selected,
             "source_total_columns":len(fields)}
        receipts.append(row)
        if i==0:
            tab=parquet.read_row_group(0,columns=list(selected.values()))
            data=tab.slice(0,min(tab.num_rows,15000)).to_pandas()
            fam=data[selected["family"]].fillna("").astype(str)
            years=pd.to_numeric(data[selected["year"]],errors="coerce")
            sample={
                "original_first_rowgroup_rows":tab.num_rows,
                "first_source_15000_records_read":len(data),
                "sample_top_family_counts":{str(k):int(n)
                    for k,n in fam.value_counts().head(25).items()},
                "sample_gecko_Gekkonidae_family_rows":int(fam.eq("Gekkonidae").sum()),
                "sample_source_2013_to_2025_year_rows":
                    int(years.between(2013,2025).sum()),
                "nonreleased_original_individual_lizard_coordinates":True}
        f.close()
    result={"schema":SCHEMA,"source_DOI":SOURCE_DOI,"original_taxon_partition":SOURCE_GROUP,
        "original_total_reptile_source_rows_all_years_all_basis":
            sum(f["original_parquet_true_rows"] for f in receipts),
        "original_compressed_reptile_Parquet_GiB":
            sum(f["original_compressed_source_bytes"] for f in receipts)/1024**3,
        "original_parquet_file_receipts":receipts,
        "first_row_group_technical_sample_NOT_global_data_estimate":sample,
        "no_individual_records_exported_or_original_source_species_presence_inferred":True,
        "source_timestamp_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_TRUE_SQUAMATA_GECKO_ORIGINAL_PARQUET_SCHEMA",
          json.dumps({"true_reptile_source_rows":
              result["original_total_reptile_source_rows_all_years_all_basis"],
              "original_source_GiB":
              result["original_compressed_reptile_Parquet_GiB"],
              "first_rowgroup_family_counts":sample["sample_top_family_counts"]},
              sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.manifest,a.out)
