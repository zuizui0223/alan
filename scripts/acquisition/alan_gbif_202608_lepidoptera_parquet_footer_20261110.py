#!/usr/bin/env python3
"""Read authentic Parquet FOOTERS of GBIF 2026-08 Lepidoptera only.

No actual occurrence rows are transferred. Public S3 footer metadata gives
true row groups, exact schema, source coverage and source bytes. This is the
cost gate for streaming tens of millions of records rather than scraping
GBIF APIs thousands of times. Preview DOI 10.15468/dl.8yrbe7.
"""
from __future__ import annotations
import argparse,json,datetime,hashlib
from pathlib import Path
import pyarrow.parquet as pq
from pyarrow.fs import S3FileSystem

BUCKET="gbif-public-data"
SCHEMA="ALAN_GBIF_2026_08_original_lepidoptera_Parquet_frozen_footers_v1"

def require(cond,msg):
    if not cond:raise ValueError("ALAN_GLOBAL_PARQUET_FOOTER: "+str(msg))

def run(manifest,out):
    out=Path(out);require(not out.exists(),"No overwrite of true source footers")
    source=json.loads(Path(manifest).read_text())
    require(source["schema"]=="ALAN_2026_08_GBIF_actual_taxon_partition_public_S3_manifest_v1"
            and source["GBIF_source_DOI"]=="10.15468/dl.8yrbe7"
            and source["lepidoptera_original_source_partitions"]==6,
            "Prior original GBIF S3 manifest not source checked")
    fs=S3FileSystem(anonymous=True,region="us-east-1")
    results={}
    for name,source_files in source["genuine_lepidoptera_parquet_partition_manifests"].items():
        files=[]
        for s in source_files["object_keys"]:
            key=s["key"]
            # S3FileSystem reads only footer and schema metadata, not rows.
            reader=fs.open_input_file(BUCKET+"/"+key)
            pf=pq.ParquetFile(reader)
            md=pf.metadata
            fields=list(pf.schema_arrow.names)
            require(md.num_rows>0 and md.num_row_groups>0,
                    "Original published source Parquet contains no observed rows")
            need=("gbifid","decimallatitude","decimallongitude",
                  "year","family","basisofrecord")
            require(all(x.casefold() in {c.casefold() for c in fields} for x in need),
                    "Native original GBIF occurrence source schema missing fields")
            files.append({"source_key":key,"source_bytes":s["bytes"],
                 "original_parquet_rows":md.num_rows,"source_row_groups":md.num_row_groups,
                 "source_footer_serialized_bytes":md.serialized_size,
                 "source_original_columns":fields[:85],
                 "source_original_species_field_present":
                      "specieskey" in {c.casefold() for c in fields},
                 "source_original_recorder_field_present":
                      "recordedby" in {c.casefold() for c in fields}})
            reader.close()
        results[name]={"original_parquet_files":files,
             "total_original_rows":sum(v["original_parquet_rows"] for v in files),
             "total_original_bytes":sum(v["source_bytes"] for v in files)}
        print("ALAN_TRUE_WORLDWIDE_LEPIDOPTERA_PARQUET_FOOTER",
              json.dumps({name:{"real_original_rows":results[name]["total_original_rows"],
                                "source_bytes":results[name]["total_original_bytes"],
                                "n_files":len(files)}},sort_keys=True),flush=True)
    rep={
       "schema":SCHEMA,"published_GBIF_snapshot_DOI":"10.15468/dl.8yrbe7",
       "published_preview_S3_root":source["GBIF_experimental_source_snapshot"],
       "source_true_taxonomic_partitions":len(results),
       "source_Parquet_original_row_total":
          sum(v["total_original_rows"] for v in results.values()),
       "source_Parquet_original_compressed_GiB":
          sum(v["total_original_bytes"] for v in results.values())/1024**3,
       "original_partition_footers":results,
       "did_not_read_any_observed_occurrence_rows":True,
       "source_basis_of_record_filters_NOT_applied_to_footer_totals":True,
       "no_island_vs_continent_point_assignment_yet":True,
       "no_moth_night_activity_or_effective_pollination_measured":True,
       "source_date_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat()
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(rep,sort_keys=True,indent=2)+"\n")
    print("ALAN_GBIF_2026_08_TRUE_LEPIDOPTERA_PARQUET_ROW_SUPPORT",
          json.dumps({"original_global_lepidoptera_parquet_rows":
                    rep["source_Parquet_original_row_total"],
                    "original_global_lepidoptera_compressed_GiB":
                    rep["source_Parquet_original_compressed_GiB"]},sort_keys=True),flush=True)
    return rep

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.manifest,a.out)
