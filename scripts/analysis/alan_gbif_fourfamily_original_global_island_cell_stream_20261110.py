#!/usr/bin/env python3
"""Stream true GBIF 2026-08 original moth/butterfly occurrences onto GSHHG islands.

Families: Noctuidae, Geometridae, Nymphalidae, Pieridae, four disjoint source
partitions. Each worker reads Parquet row groups lazily from genuine public
GBIF source (2026-08 DOI 10.15468/dl.8yrbe7); never downloads images.
Restrict HUMAN_OBSERVATION, 2013-2025 actual year, valid coordinates,
PRESENT, exact source family, GSHHG 2.3.7 islands area >=20km2, using actual
point-in-polygon before 0.05° cell/year aggregation. Effort/detection-biased
RECORD counts, not insect population abundance or nocturnal pollen service.
"""
from __future__ import annotations
import argparse,collections,datetime,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import pyarrow.parquet as pq
from pyarrow.fs import S3FileSystem

from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land

FAMILIES=("Noctuidae","Geometridae","Nymphalidae","Pieridae")
SCHEMA="ALAN_GBIF_202608_true_fourfamily_GSHHG_island_cell_year_occurrences_v1"
DOI="10.15468/dl.8yrbe7"
YEAR_LOW,YEAR_HIGH=2013,2025
MIN_ISLAND_AREA_KM2=20.
GRID=0.05
BUCKET="gbif-public-data"

def require(ok,reason):
    if not ok:raise ValueError("ALAN_GBIF_GLOBAL_ISLAND_STREAM: "+str(reason))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):
            h.update(block)
    return h.hexdigest()

def mask_for_island_bbox(land):
    mask=np.zeros((180,360),dtype=np.bool_)
    for bounds in land.geometry.bounds.itertuples(index=False,name=None):
        a,b,c,d=map(float,bounds)
        if not np.isfinite([a,b,c,d]).all():continue
        lo=max(0,min(359,int(math.floor(a+180))))
        hi=max(0,min(359,int(math.floor(c+180))))
        bottom=max(0,min(179,int(math.floor(b+90))))
        top=max(0,min(179,int(math.floor(d+90))))
        if lo<=hi and bottom<=top:
            mask[bottom:top+1,lo:hi+1]=True
    return mask

def original_source_files(receipt,family):
    raw=json.loads(Path(receipt).read_text())
    require(raw.get("schema")=="ALAN_GBIF_2026_08_original_lepidoptera_Parquet_frozen_footers_v1"
            and raw.get("published_GBIF_snapshot_DOI")==DOI,
            "Original worldwide taxonomy+source source lineage not signed")
    candidates=[(n,d) for n,d in raw["original_partition_footers"].items()
                if n.endswith("_"+family)]
    require(len(candidates)==1,"Expected one verified taxon partition "+family)
    name,meta=candidates[0]
    require(len(meta["original_parquet_files"])>=1,"No source parquet")
    return name,meta["original_parquet_files"],meta["total_original_rows"]

def required_column(names,name):
    candidates=[v for v in names if v.casefold()==name.casefold()]
    require(len(candidates)==1,"Original GBIF schema missing expected column: "+name)
    return candidates[0]

def scan(family,receipt,coast,out,max_row_groups=None):
    out=Path(out);require(not out.exists(),"Original source receipt already exists")
    require(family in FAMILIES,"Unreviewed taxonomic guild or family")
    partition,files,expected_rows=original_source_files(receipt,family)
    land=load_signed_land(coast)
    island=land.loc[land.land_class.eq("island") &
                    land.area_km2.ge(MIN_ISLAND_AREA_KM2),
                    ["land_id","area_km2","geometry"]].copy()
    require(island.land_id.is_unique and len(island)>1000,
            "Original GSHHG 2.3.7 island source not complete")
    island["land_id"]=island.land_id.astype(str)
    mask=mask_for_island_bbox(island)
    counters=collections.Counter()
    unique_species=collections.defaultdict(set)
    tallies=collections.Counter()
    s3=S3FileSystem(anonymous=True,region="us-east-1")
    actual_file_meta=[]
    stop=False
    for path in files:
        if stop:break
        key=path["source_key"]
        reader=s3.open_input_file(BUCKET+"/"+key)
        parquet=pq.ParquetFile(reader)
        require(parquet.metadata.num_rows>0 and
                parquet.metadata.num_rows==path["original_parquet_rows"],
                "Parquet file counts drift vs verified source footer")
        names=list(parquet.schema_arrow.names)
        cols={v:required_column(names,v) for v in (
            "gbifid","year","decimallatitude","decimallongitude",
            "basisofrecord","occurrencestatus","family","specieskey")}
        for rg in range(parquet.metadata.num_row_groups):
            if max_row_groups is not None and tallies["source_original_row_groups_read"]>=max_row_groups:
                stop=True;break
            tab=parquet.read_row_group(rg,columns=list(cols.values()))
            actual=tab.to_pandas()
            tallies["source_original_records_read"]+=len(actual)
            tallies["source_original_row_groups_read"]+=1
            # Exact source eligibility; never turn missing year/site into a zero count.
            year=pd.to_numeric(actual[cols["year"]],errors="coerce")
            mask1=(year.between(YEAR_LOW,YEAR_HIGH).to_numpy() &
                actual[cols["basisofrecord"]].fillna("").astype(str)
                      .str.upper().eq("HUMAN_OBSERVATION").to_numpy() &
                actual[cols["occurrencestatus"]].fillna("").astype(str)
                      .str.upper().eq("PRESENT").to_numpy() &
                actual[cols["family"]].fillna("").astype(str)
                      .str.casefold().eq(family.casefold()).to_numpy())
            d=actual.loc[mask1]
            tallies["human_observation_in_range_exact_family"]+=len(d)
            if d.empty:continue
            lat=pd.to_numeric(d[cols["decimallatitude"]],errors="coerce").to_numpy()
            lon=pd.to_numeric(d[cols["decimallongitude"]],errors="coerce").to_numpy()
            ok=np.isfinite(lat)&np.isfinite(lon)&(lat>=-90)&(lat<=90)&(lon>=-180)&(lon<=180)
            d=d.loc[ok].copy()
            lat=lat[ok];lon=lon[ok]
            tallies["with_valid_original_coordinates"]+=len(d)
            if d.empty:continue
            ix=np.minimum(359,np.maximum(0,np.floor(lon+180).astype(int)))
            iy=np.minimum(179,np.maximum(0,np.floor(lat+90).astype(int)))
            nearby=mask[iy,ix]
            d=d.loc[nearby].copy()
            lat=lat[nearby];lon=lon[nearby]
            tallies["candidate_locations_near_original_islands"]+=len(d)
            if d.empty:continue
            d["observed_year"]=year.loc[d.index].astype(int).to_numpy()
            d["species_token"]=d[cols["specieskey"]].fillna("").astype(str)
            d["grid_x"]=np.minimum(7199,np.floor((lon+180)/GRID).astype(int))
            d["grid_y"]=np.minimum(3599,np.floor((lat+90)/GRID).astype(int))
            loc=gpd.GeoDataFrame(
                d[["observed_year","species_token","grid_x","grid_y"]],
                geometry=gpd.points_from_xy(lon,lat),crs="EPSG:4326")
            hit=gpd.sjoin(loc,island[["land_id","geometry"]],how="inner",
                          predicate="within")
            # An original source point on overlapping GSHHG island polygons
            # must not count multiple times.
            dup=hit.index.duplicated(keep=False)
            tallies["ambiguous_overlapping_island_assignments"]+=int(dup.sum())
            hit=hit.loc[~dup].copy()
            tallies["real_2026GBIF_island_geo_observations"]+=len(hit)
            if not len(hit):continue
            keys=["land_id","observed_year","grid_x","grid_y"]
            n=hit.groupby(keys,sort=False).size()
            for k,v in n.items():counters[k]+=int(v)
            sampled=hit.loc[hit.species_token.ne(""),keys+["species_token"]].drop_duplicates()
            for k,g in sampled.groupby(keys,sort=False):
                unique_species[k].update(g["species_token"].tolist())
            if tallies["source_original_row_groups_read"]%25==0:
                print("ALAN_GBIF_TRUE_SOURCE_ISLAND_SCAN_PROGRESS",
                      json.dumps({"family":family,
                        "source_records":tallies["source_original_records_read"],
                        "actual_island_points":tallies["real_2026GBIF_island_geo_observations"],
                        "unique_cell_years":len(counters)},sort_keys=True),flush=True)
        actual_file_meta.append({"source_key":key,"original_rows":parquet.metadata.num_rows,
                                 "original_row_groups":parquet.metadata.num_row_groups,
                                 "actual_source_row_groups_read_to_date":tallies["source_original_row_groups_read"]})
        reader.close()
    complete=max_row_groups is None
    if complete:
        require(tallies["source_original_records_read"]==expected_rows,
                "Did not read every original occurrence in this taxon partition")
    rows=[]
    for (island_id,year,gx,gy),n in sorted(counters.items()):
        rows.append({"island":island_id,"year":year,"latbin":int(gy),
                     "lonbin":int(gx),"family":family,
                     "original_GBIF_human_occurrence_records":int(n),
                     "distinct_original_species_keys_recorded":len(unique_species[(island_id,year,gx,gy)])})
    frame=pd.DataFrame(rows,columns=["island","year","latbin","lonbin","family",
        "original_GBIF_human_occurrence_records","distinct_original_species_keys_recorded"])
    require(not frame.duplicated(["island","year","latbin","lonbin","family"]).any(),
            "Original island source cell/year pseudo-replication")
    out.mkdir(parents=True)
    path=out/(family+"_GBIF_202608_original_island_photo_cell_year.csv.gz")
    frame.to_csv(path,index=False,compression="gzip")
    audit={
        "schema":SCHEMA,"GBIF_source_2026_08_DOI":DOI,
        "source_original_family":family,"source_taxon_partition":partition,
        "source_expected_family_original_all_year_all_basis_records":expected_rows,
        "source_original_parquet_files":actual_file_meta,
        "source_scan_complete_all_files_all_row_groups":complete,
        "source_scan_max_row_groups_only":max_row_groups,
        "source_actual_scan":dict(tallies),
        "original_GSHHG_islands_ge_20km2":len(island),
        "actual_source_island_cell_year_family_rows":len(frame),
        "actual_observed_source_islands":int(frame.island.nunique()),
        "actual_observed_source_years":sorted(frame.year.unique().astype(int).tolist()) if len(frame) else [],
        "source_island_level_photo_observation_record_total":int(
            frame.original_GBIF_human_occurrence_records.sum()),
        "source_aggregated_cell_CSV_SHA256":sha(path),
        "source_coordinates_not_exported_only_coarse_0p05bins":True,
        "actual_source_taxa_necessary_but_not_all_moth_or_butterfly_families":True,
        "different_GBIF_records_not_independent_insect_population_density":True,
        "source_original_moth_observations_not_verified_night_activity_or_floral_visits":True,
        "photo_observer_effort_and_GBIF_study_source_not_adjusted":True,
        "no_population_decline_or_causal_ALAN_effect_estimated":True,
        "UTC_executed":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (out/(family+"_GBIF_202608_real_island_source_audit.json")).write_text(
        json.dumps(audit,sort_keys=True,indent=2)+"\n")
    print("ALAN_GBIF_REAL_GLOBAL_ISLAND_FAMILY_SOURCE",
          json.dumps({k:audit[k] for k in (
              "source_original_family","source_scan_complete_all_files_all_row_groups",
              "source_actual_scan","original_GSHHG_islands_ge_20km2",
              "actual_source_island_cell_year_family_rows",
              "actual_observed_source_islands",
              "source_island_level_photo_observation_record_total")},
              sort_keys=True),flush=True)
    return audit

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--family",choices=FAMILIES,required=True)
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--coast",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--max-row-groups",type=int,default=None)
    a=p.parse_args()
    scan(a.family,a.manifest,a.coast,a.out,a.max_row_groups)
