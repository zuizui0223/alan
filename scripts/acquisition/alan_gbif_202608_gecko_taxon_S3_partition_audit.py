#!/usr/bin/env python3
"""Find the real GBIF 2026-08 DWCA Parquet reptile/gecko source partition.

List only public S3 object METADATA (names, object size, ETag), no individual
gecko/moth coordinates, no wildlife capture records or gecko presence guesses.
A taxon not partitioned separately may reside in a broader Reptilia,
Squamata or PARTIAL partition; infer neither species absence nor data gap
from absence of a Gekkonidae-named partition.
"""
from __future__ import annotations
import argparse,collections,datetime,json
from pathlib import Path
import boto3
from botocore.config import Config
from botocore import UNSIGNED

BUCKET="gbif-public-data"
ROOT="development/2026-08-01-DWCA/p_taxon/"
SCHEMA="ALAN_original_GBIF_202608_Gekkota_reptile_Parquet_partition_S3_audit_v1"
DOI="10.15468/dl.8yrbe7"
KEYWORDS=("reptilia","squamata","gekk","hemidactyl","phelsuma",
          "tarentola","sauria","lizard","scincidae","amphibia")

def require(x,reason):
    if not x: raise ValueError("ALAN_GECKO_S3_SOURCE: "+str(reason))

def run(out):
    out=Path(out)
    require(not out.exists(),"No overwrite genuine source manifest")
    s3=boto3.client("s3",region_name="us-east-1",config=Config(
        signature_version=UNSIGNED,connect_timeout=20,read_timeout=60,
        retries={"max_attempts":3}))
    bins=collections.defaultdict(list)
    n=0
    for p in s3.get_paginator("list_objects_v2").paginate(
        Bucket=BUCKET,Prefix=ROOT,PaginationConfig={"PageSize":1000}):
        for o in p.get("Contents",[]):
            key=o["Key"]
            if not key.endswith(".parquet"):continue
            suffix=key[len(ROOT):]
            group=suffix.split("/")[0]
            require(bool(group) and "/" in suffix,
                    "Unknown original taxon partition naming")
            n+=1
            bins[group].append({
                "key":key,
                "compressed_original_bytes":int(o["Size"]),
                "ETag_not_guaranteed_MD5":str(o["ETag"]).strip('"')})
    matching={name:{"source_parquet_files":len(obj),
       "source_GiB":sum(x["compressed_original_bytes"] for x in obj)/1024**3,
       "original_source_object_metadata":obj}
       for name,obj in sorted(bins.items())
       if any(w in name.lower() for w in KEYWORDS)}
    require(n>100 and len(bins)>10,"Published original GBIF source S3 unavailable")
    report={"schema":SCHEMA,"original_GBIF_snapshot_DOI":DOI,
       "source_bucket":BUCKET,"source_root":ROOT,
       "original_S3_objects_listed":n,
       "original_taxon_folders_total":len(bins),
       "matching_reptile_gecko_partition_folders":matching,
       "total_matching_reptile_gecko_GiB":sum(
            v["source_GiB"] for v in matching.values()),
       "no_2026_08_occurrence_rows_downloaded":True,
       "separate_gecko_family_folder_not_required_for_gecko_source_presence":True,
       "S3_ETags_not_original_MD5_guarantees":True,
       "source_time_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_GBIF_TRUE_GECKO_S3_PARTITION_SOURCE",json.dumps({
        "original_taxon_folders_total":report["original_taxon_folders_total"],
        "matching_partition_folders":{k:{"files":v["source_parquet_files"],
                "source_GiB":v["source_GiB"]}
         for k,v in matching.items()}},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.out)
