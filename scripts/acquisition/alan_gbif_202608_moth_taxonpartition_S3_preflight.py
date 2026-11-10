#!/usr/bin/env python3
"""Inspect actual 2026-08 GBIF taxon-partitioned Parquet S3 object manifests.

Zero occurrence records downloaded. This gives true source bytes, partition
names, and tractability for worldwide Lepidoptera, rather than promising
impossible GBIF API pagination. Experimental 2026-08 DWCA source may differ
from ordinary GBIF monthly snapshots; eventual citations/quality gates
must use its exact 2026-08 DOI. No biological claim.
"""
from __future__ import annotations
import argparse, collections, datetime, hashlib, json
from pathlib import Path
import boto3
from botocore import UNSIGNED
from botocore.config import Config

BUCKET="gbif-public-data"
ROOT="development/2026-08-01-DWCA/p_taxon/"
DOI="10.15468/dl.8yrbe7"
SCHEMA="ALAN_2026_08_GBIF_actual_taxon_partition_public_S3_manifest_v1"

def require(x,message):
    if not x:raise ValueError("ALAN_GBIF_S3_TAXON_PREFLIGHT: "+str(message))

def run(out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite a source audit manifest")
    client=boto3.client("s3",region_name="us-east-1",
          config=Config(signature_version=UNSIGNED,retries={"max_attempts":3},
                        connect_timeout=20,read_timeout=50))
    items=[]
    for page in client.get_paginator("list_objects_v2").paginate(
            Bucket=BUCKET,Prefix=ROOT,PaginationConfig={"PageSize":1000}):
        for v in page.get("Contents",[]):
            key=str(v["Key"])
            if not key.endswith(".parquet"):continue
            items.append({"s3_key":key,"size_bytes":int(v["Size"]),
                          "etag":str(v["ETag"]).strip('"'),
                          "last_modified_UTC":v["LastModified"].isoformat()})
    require(len(items)>0,"Official GBIF development snapshot not listable; no invented table")
    groups=collections.defaultdict(lambda:{"source_files":0,"source_bytes":0,
                                           "object_keys":[]})
    for x in items:
        key=x["s3_key"];bits=key.split("/")
        part=next((v.split("=",1)[1] for v in bits
                   if v.startswith("taxonPartition=")),None)
        require(part is not None,"GBIF listing lacks partition key")
        d=groups[part]
        d["source_files"]+=1
        d["source_bytes"]+=x["size_bytes"]
        d["object_keys"].append({"key":key,"bytes":x["size_bytes"],
                                 "etag_not_guaranteed_md5":x["etag"]})
    moth={k:v for k,v in groups.items()
          if "Lepidoptera" in k}
    require(len(moth)>0,"No genuine Lepidoptera source taxon partitions found")
    total_bytes=sum(v["source_bytes"] for v in groups.values())
    moth_bytes=sum(v["source_bytes"] for v in moth.values())
    rep={
        "schema":SCHEMA,"GBIF_experimental_source_snapshot":ROOT,
        "GBIF_source_DOI":DOI,"source_S3_bucket":BUCKET,
        "source_time_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "listable_real_parquet_source_objects":len(items),
        "all_original_source_taxon_partitions":len(groups),
        "total_original_snapshot_listing_GiB":total_bytes/1024**3,
        "lepidoptera_original_source_partitions":len(moth),
        "lepidoptera_total_source_GiB":moth_bytes/1024**3,
        "genuine_lepidoptera_parquet_partition_manifests":dict(sorted(moth.items())),
        "has_original_site_observer_effort_factors_unverified":True,
        "no_original_occurrence_records_downloaded_or_analyzed":True,
        "no_coarse_georeferenced_moth_island_counts_claimed":True,
        "GBIF_2026_08_DWCA_IS_experimental_preview":True,
        "S3_ETag_may_be_multipart_not_MD5_checksum":True,
        "2026_10_normal_AWS_GBIF_snapshots_reported_incomplete_in_GBIF_issue_609":True
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(rep,sort_keys=True,indent=2)+"\n")
    print("ALAN_GBIF_ACTUAL_WORLDWIDE_MOTH_PARTITION_SOURCE",
          json.dumps({k:rep[k] for k in (
              "listable_real_parquet_source_objects",
              "all_original_source_taxon_partitions",
              "total_original_snapshot_listing_GiB",
              "lepidoptera_original_source_partitions",
              "lepidoptera_total_source_GiB",
          )},sort_keys=True),flush=True)
    return rep

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    run(args.out)
