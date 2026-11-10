#!/usr/bin/env python3
"""Global gecko source and predator-moth taxon audit (not ecological effects).

Source-count true GBIF georeferenced 2013-25 HUMAN_OBSERVATION gecko records,
plus genuinely indexed 2026-08 frozen preview Parquet gecko/reptile partitions.
Don't mistake all Gekkonidae (diurnal genera exist) for obligate nocturnal
moth-eaters, a GBIF gecko occurrence for measured predation, or missing records
for predator absence. This is data-access/coverage, not a causal model.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import boto3
from botocore import UNSIGNED
from botocore.config import Config

GBIF="https://api.gbif.org/v1"
BUCKET="gbif-public-data"
ROOT="development/2026-08-01-DWCA/p_taxon/"
GBIF_DOI="10.15468/dl.8yrbe7"
YEAR_MIN,YEAR_MAX=2013,2025
TARGETS={
 "Gekkonidae":"FAMILY",
 "Hemidactylus":"GENUS",
 "Gekko":"GENUS",
 "Gehyra":"GENUS",
 "Lepidodactylus":"GENUS",
 "Phelsuma":"GENUS",
 "Hemidactylus frenatus":"SPECIES",
 "Gekko japonicus":"SPECIES",
 "Gehyra dubia":"SPECIES",
}
COUNTRIES=("JP","TW","PH","ID","AU","NZ","MG","MU","SC","FJ","NC",
           "US","CU","BS","GR","ES","PT")
SCHEMA="ALAN_real_Gekkonidae_nightlight_moth_global_source_opportunity_20261110_v1"


def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_GECKO_SOURCE: "+str(msg))


def api_json(endpoint,params):
    q=urllib.parse.urlencode(params,doseq=True)
    url=f"{GBIF}{endpoint}?{q}"
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-academic-global-gecko-moth-source-audit/1.0",
        "Accept":"application/json"})
    for trial in range(4):
        try:
            with urllib.request.urlopen(req,timeout=55) as r:
                require(r.status==200,"GBIF source not HTTP 200")
                raw=r.read()
            return json.loads(raw),url,hashlib.sha256(raw).hexdigest()
        except (OSError,TimeoutError,ValueError):
            if trial==3:raise
            time.sleep(trial+1)
    raise RuntimeError("No original public GBIF source")


def exact_taxon(name,rank):
    response,url,checksum=api_json("/species/match",{
        "name":name,"strict":"true"})
    if response.get("usageKey") and response.get("matchType")!="NONE" \
            and str(response.get("rank","")).upper()==rank:
        return {"taxonKey":int(response["usageKey"]),"source_name":name,
                "rank":rank,"original_query":url,
                "original_source_sha256":checksum,"resolution":"exact_rank_match"}
    # A genus/species may be absent or rank-ambiguous in some source catalogues.
    # Do not silently assign a parent/family key or impute predator behavior.
    return {"source_name":name,"rank":rank,"source_resolution":"not_exactly_resolved",
            "original_source_response":response,"original_query":url}


def count_source(key,country=None):
    query=[("taxonKey",int(key)),("basisOfRecord","HUMAN_OBSERVATION"),
           ("hasCoordinate","true"),("occurrenceStatus","PRESENT"),
           ("year",f"{YEAR_MIN},{YEAR_MAX}"),("limit","0")]
    if country:query.append(("country",country))
    source,url,checksum=api_json("/occurrence/search",query)
    require(isinstance(source.get("count"),int) and
            source["count"]>=0,"Missing true GBIF gecko occurrence count")
    return {"GBIF_original_georeferenced_human_records":int(source["count"]),
            "source_url":url,"source_response_sha256":checksum}


def s3_partition_census():
    cli=boto3.client("s3",region_name="us-east-1",
        config=Config(signature_version=UNSIGNED,
                      retries={"max_attempts":3},
                      connect_timeout=20,read_timeout=45))
    matches=[]
    listed=0
    parts=set()
    for page in cli.get_paginator("list_objects_v2").paginate(
            Bucket=BUCKET,Prefix=ROOT,PaginationConfig={"PageSize":1000}):
        for item in page.get("Contents",[]):
            key=item["Key"]
            if not key.endswith(".parquet"):continue
            listed+=1
            folder=key[len(ROOT):].split("/")[0]
            parts.add(folder)
            if any(x in folder.casefold() for x in (
                "gekko","gekkon","hemidactyl","gehyra","lepidodactyl",
                "phelsuma","reptilia","squamata")):
                matches.append({
                   "source_GBIF_parquet_object":key,
                   "source_taxon_partition":folder,
                   "source_compressed_bytes":int(item["Size"]),
                   "etag_not_guaranteed_md5":item["ETag"].strip('"')})
    return {"s3_gbif_snapshot":ROOT,"source_DOI":GBIF_DOI,
            "total_source_parquet_files_all_partitions":listed,
            "total_source_taxonomic_partitions":len(parts),
            "reptile_or_gecko_source_partition_objects":matches,
            "gecko_specific_original_source_files":sum(
               "gekkon" in x["source_taxon_partition"].casefold()
               for x in matches),
            "no_source_row_census_performed":True}


def run(output):
    output=Path(output)
    require(not output.exists(),"Do not overwrite authentic gecko source receipt")
    taxonomy={}
    source={}
    for name,rank in TARGETS.items():
        taxon=exact_taxon(name,rank)
        taxonomy[name]=taxon
        if taxon.get("taxonKey"):
            source[name]=count_source(taxon["taxonKey"])
        else:source[name]={"status":"taxon_not_unambiguously_resolved"}
        print("ALAN_ACTUAL_GBIF_GECKO_TAXON_SOURCE",json.dumps({
            "taxon":name,"rank":rank,
            "taxonomy":{k:v for k,v in taxon.items()
                        if k!="original_source_response"},
            "record_source":source[name]},sort_keys=True),flush=True)
    family=taxonomy["Gekkonidae"].get("taxonKey")
    require(bool(family),"Cannot source-verify Gekkonidae original family")
    countries={}
    for country in COUNTRIES:
        countries[country]=count_source(family,country)
        print("ALAN_TRUE_GECKONID_COUNTRY_SOURCE",
              json.dumps({country:countries[country]["GBIF_original_georeferenced_human_records"]},
                         sort_keys=True),flush=True)
    frozen=s3_partition_census()
    report={
      "schema":SCHEMA,"query_source":"public actual GBIF LIVE occurrence search 2026-11-10",
      "queried_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "original_taxon_source_provenance":taxonomy,
      "global_taxa_2013to2025_georeferenced_human_record_counts":source,
      "countries_are_NOT_true_GSHHG_island_records":countries,
      "frozen_source_original_2026_08_DWCA_parquet":frozen,
      "Gekkonidae_includes_day_active_species":True,
      "species_activity_and_native_invasive_status_REQUIRE_independent_curated_traits":True,
      "Hemidactylus_frenatus_is_not_all_geckos":True,
      "past_moth_record_observer_effort_not_automatically_gecko_effort":True,
      "GBIF_gecko_presence_record_not_direct_predation_observation":True,
      "gecko_source_no_record_NOT_predator_absence":True,
      "countryside_and_cities_confounded_with_ALAN":True,
      "no_gecko_x_moth_x_ALAN_statistical_association_claimed":True,
      "no_causal_predation_or_night_pollination_inferred":True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_REAL_GLOBAL_GECKO_SOURCE_SUMMARY",json.dumps({
      "global_Gekkonidae_source_records":
        source["Gekkonidae"].get("GBIF_original_georeferenced_human_records"),
      "snapshot_gecko_specific_parquet_files":
        frozen["gecko_specific_original_source_files"],
      "snapshot_reptile_gecko_related_parquet_files":
        len(frozen["reptile_or_gecko_source_partition_objects"]),
      "source_countries_queried":len(countries)
    },sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.out)
