#!/usr/bin/env python3
"""Check actual eventTime availability in original 2026-08 GBIF Squamata rows.

Same three original 11,634,053 Squamata Parquet objects used in ALAN PR 39.
Counts are **source-field completeness**, not nocturnal behavior. No precise
coordinates, timestamps, individual records or observer names are exported.
Original Darwin Core eventTime clock values may lack time-zone semantics;
ISO eventDate midnight may be generated from date-only reports.
"""
from __future__ import annotations
import argparse,collections,datetime,json,re
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pyarrow.fs import S3FileSystem
from scripts.analysis.alan_global_GBIF_gecko_island_original_stream_20261110 import GECKO_FAMILIES

DOI="10.15468/dl.8yrbe7"
PARTITION="taxonpartition=Animalia_Chordata_Reptilia_Squamata"
SCHEMA="ALAN_Squamata_2013_25_original_eventTime_photo_record_completeness_v1"
FIELDS=("family","genus","year","basisofrecord","occurrencestatus",
        "decimallatitude","decimallongitude","eventtime","eventdate")
def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_GECKO_ORIGINAL_EVENTCLOCK: "+str(msg))
def actual_field(allnames,key):
    found=[k for k in allnames if k.casefold()==key.casefold()]
    require(len(found)==1,"Missing original GBIF time-source field "+key)
    return found[0]

def run(manifest,out):
    out=Path(out);require(not out.exists(),"Do not replace original true source clock count")
    m=json.loads(Path(manifest).read_text())
    require(m["schema"]=="ALAN_original_GBIF_202608_Gekkota_reptile_Parquet_partition_S3_audit_v1"
            and m["original_GBIF_snapshot_DOI"]==DOI,
            "Original 2026-08 source Squamata not authenticated")
    items=m["matching_reptile_gecko_partition_folders"][PARTITION]["original_source_object_metadata"]
    fstore=S3FileSystem(anonymous=True,region="us-east-1")
    global_counts=collections.Counter()
    groups=collections.defaultdict(collections.Counter)
    for original in items:
        opened=fstore.open_input_file("gbif-public-data/"+original["key"])
        parquet=pq.ParquetFile(opened)
        cols=list(parquet.schema_arrow.names)
        mapping={name:actual_field(cols,name) for name in FIELDS}
        for rg in range(parquet.metadata.num_row_groups):
            a=parquet.read_row_group(rg,columns=list(mapping.values())).to_pandas()
            global_counts["all_reptile_source_rows"]+=len(a)
            global_counts["original_source_row_groups_read"]+=1
            year=pd.to_numeric(a[mapping["year"]],errors="coerce")
            family=a[mapping["family"]].fillna("").astype(str)
            is_species=family.isin(tuple(GECKO_FAMILIES)+("Scincidae",))
            good=(is_species.to_numpy()&
                year.between(2013,2025).to_numpy()&
                a[mapping["basisofrecord"]].fillna("").astype(str).str.upper()
                 .eq("HUMAN_OBSERVATION").to_numpy()&
                a[mapping["occurrencestatus"]].fillna("").astype(str).str.upper()
                 .eq("PRESENT").to_numpy())
            x=a.loc[good].copy()
            global_counts["2026GBIF_source_gecko_scink_photo_candidate_rows"]+=len(x)
            if x.empty:continue
            latitude=pd.to_numeric(x[mapping["decimallatitude"]],errors="coerce")
            longitude=pd.to_numeric(x[mapping["decimallongitude"]],errors="coerce")
            geo=latitude.between(-90,90)&longitude.between(-180,180)
            x=x.loc[geo].copy()
            global_counts["gecko_scink_valid_original_georeferenced_human_records"]+=len(x)
            if x.empty:continue
            fam=x[mapping["family"]].fillna("").astype(str)
            genus=x[mapping["genus"]].fillna("").astype(str)
            event=x[mapping["eventtime"]].fillna("").astype(str).str.strip()
            date=x[mapping["eventdate"]].fillna("").astype(str).str.strip()
            useful=event.ne("")&~event.str.lower().isin(("nan","none","null","n/a"))
            hours=useful & event.str.match(
                r"^(?:[01]\d|2[0-3]):[0-5]\d(?::[0-5]\d(?:\.\d+)?)?(?:Z|[+-][0-2]\d:?\d{2})?$",
                na=False)
            explicit_tz=hours & event.str.contains(
                r"(?:Z|[+-][0-2]\d:?\d{2})$",regex=True,na=False)
            date_time=date.str.contains(r"T\d{2}:\d{2}",regex=True,na=False)
            date_has_zone=date.str.contains(
                r"(?:Z|[+-][0-2]\d:?\d{2})$",regex=True,na=False)
            date_00=date.str.contains(r"T00:00:00(?:Z|[+-]\d\d:?\d\d)?$",regex=True,na=False)
            # Count only source fields, no timestamp extraction or local-time assumptions.
            for cat,mask in (
                ("gecko",fam.isin(GECKO_FAMILIES)),
                ("scink",fam.eq("Scincidae")),
                ("Hemidactylus",genus.eq("Hemidactylus")),
                ("Gekko",genus.eq("Gekko")),
                ("Phelsuma",genus.eq("Phelsuma"))):
                entry=groups[cat]
                entry["geo_source_records"]+=int(mask.sum())
                entry["original_eventtime_nonempty"]+=int((mask&useful).sum())
                entry["eventtime_valid_HHMM_format"]+=int((mask&hours).sum())
                entry["eventtime_explicit_zone"]+=int((mask&explicit_tz).sum())
                entry["eventdate_has_T_clock"]+=int((mask&date_time).sum())
                entry["eventdate_has_timezone_suffix"]+=int((mask&date_has_zone).sum())
                entry["eventdate_midnight_literal"]+=int((mask&date_00).sum())
        opened.close()
    require(global_counts["all_reptile_source_rows"]==11634053,
            "Not all genuine 2026-08 source reptile occurrence rows were inspected")
    result={
        "schema":SCHEMA,"original_GBIF_snapshot_DOI":DOI,
        "source_Squamata_all_basis_all_year_rows_read":
            global_counts["all_reptile_source_rows"],
        "original_source_row_groups_read":
            global_counts["original_source_row_groups_read"],
        "original_source_selection_quality":dict(global_counts),
        "source_family_and_genus_photo_eventtime_completeness":
            {k:dict(v) for k,v in groups.items()},
        "source_counts_are_georeferenced_record_types_not_exact_islands":True,
        "field_completeness_is_NOT_gecko_lamp_visitation_or_predation":True,
        "event_date_T_midnight_can_be_generated_default":True,
        "event_time_hhmm_without_timezone_is_not_verified_local_solar_time":True,
        "no_original_individual_times_locations_or_observers_exported":True,
        "source_UTC_audit":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_GECKO_SOURCE_EVENTTIME_COVERAGE",
      json.dumps(result["source_family_and_genus_photo_eventtime_completeness"],
       sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.manifest,a.out)
