#!/usr/bin/env python3
"""Globally matched *prior-year* PCNL at real GBIF island photo-grid cells.

Input: 2026-08 GBIF source SHA-locked four-family 2013–25 observations,
22,713,063 actual source records on GSHHG island polygons >=20km2.
Observation unit: original island × 0.05deg cell × observed year, NOT
independent moth individual, species absence or insect density.

Satellite cell anchors use literal 0.05deg cell CENTERS derived from original
source grid IDs. Centers on sea or a different island are explicitly NA:
this is a strict source-land center eligibility, not a claim that actual GBIF
occurrences were outside land. Ocean, missing rasters and future years never
turn into fake dark-light exposure or zero moth counts.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
import shapely
from rasterio.warp import transform

from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land
from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one,manifest_read)
SCHEMA="ALAN_GBIF_global_fourfamily_2013_2025_real_island_prior_12_24_PCNL_v1"
SOURCE_SCHEMA="ALAN_GBIF_2026_08_4family_global_GSHHG_island_yearcell_real_counts_v1"
SOURCE_SHA="a993c0908600211e3e0ff66b88920ad0d47124e122f373bbce0e78ba9902ee6d"
OBS_YEARS=list(range(2013,2026))
NROWS=228361
NOBS=22713063
SHARDS=4
SOURCE_PATHNAME="GBIF_202608_actual_global_island_0p05deg_year_4family_occurrence_counts.csv.gz"
SOURCE_RECEIPT="ALAN_2026_08_real_global_fourfamily_island_counts_receipt.json"
MANIFEST="config/alan_pcnl_1992_2024_zenodo_v3_md5.json"

def require(x,message):
    if not bool(x):raise ValueError("ALAN_GLOBAL_PRIORLIGHT: "+str(message))

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def years(shard):
    require(0<=shard<SHARDS,"Unrecognized original historical year shard")
    return OBS_YEARS[shard::SHARDS]

def load_original_global(csv,receipt):
    p=json.loads(Path(receipt).read_text())
    require(p["schema"]==SOURCE_SCHEMA and
            p["source_grid_csv_SHA256"]==SOURCE_SHA and
            p["observed_original_island_0p05_cell_years"]==NROWS and
            p["original_fourfamily_record_count_on_actual_islands"]==NOBS and
            p["observed_actual_distinct_islands"]==1217 and
            p["at_least_5_records_and_both_guilds_cell_year"]==70219 and
            sha256(csv)==SOURCE_SHA,
            "Source is NOT exactly original global four-family SHA signed record")
    d=pd.read_csv(csv,dtype={"island":str},low_memory=False)
    keys=["island","year","latbin","lonbin"]
    required=set(keys+["moth_family_records","butterfly_family_records",
                       "all_four_family_records",
                       "relative_moth_family_record_fraction",
                       "has_both_moth_and_butterfly_family_records",
                       "has_at_least_5_four_family_records"])
    require(len(d)==NROWS and required.issubset(d.columns) and
            not d.duplicated(keys).any() and
            set(d.year)==set(OBS_YEARS) and
            int(d.all_four_family_records.sum())==NOBS,
            "Year/cell/source/GloBI-free GBIF record counts drift")
    require((d.moth_family_records+d.butterfly_family_records==
             d.all_four_family_records).all(),
            "Not a true moth-plus-butterfly original source composition")
    d["longitude_center"]=(d.lonbin.to_numpy(float)+.5)*.05-180.
    d["latitude_center"]=(d.latbin.to_numpy(float)+.5)*.05-90.
    require(d.longitude_center.between(-180,180).all() and
            d.latitude_center.between(-90,90).all(),
            "Photo grid center formula changed")
    return d

def same_island_photo_grid_centers(d,coast):
    land=load_signed_land(coast)
    land=land.set_index("land_id")
    islands=d.island.unique()
    require(set(islands).issubset(set(land.index)) and
            land.loc[islands,"land_class"].eq("island").all() and
            land.loc[islands,"area_km2"].ge(20).all(),
            "GSHHG source family source-island identity/area mismatched")
    d=d.copy()
    d["area_km2"]=d.island.map(land.area_km2)
    inside=np.zeros(len(d),dtype=bool)
    # Vectorized per original polygon; no country centroid, ocean or other
    # nearby island is silently promoted to a same-island source cell.
    for island,idx in d.groupby("island",sort=False).indices.items():
        part=d.iloc[idx]
        polygon=land.loc[island,"geometry"]
        within=shapely.contains_xy(
            polygon,
            part.longitude_center.to_numpy(float),
            part.latitude_center.to_numpy(float))
        inside[idx]=within
    d["center_inside_same_original_GSHHG_island"]=inside
    return d

def genuine_values(source,d):
    require(source.crs is not None and source.count==1,
            "Original PCNL source band and CRS invalid")
    x=d.longitude_center.to_numpy(float)
    y=d.latitude_center.to_numpy(float)
    if source.crs.to_epsg()!=4326:
        x,y=transform("EPSG:4326",source.crs,x.tolist(),y.tolist())
    data=np.full(len(d),np.nan,dtype=float)
    for start in range(0,len(d),10000):
        samples=source.sample(list(zip(x[start:start+10000],y[start:start+10000])),
                              indexes=1,masked=True)
        for j,arr in enumerate(samples,start=start):
            value=arr[0]
            if np.ma.is_masked(value):continue
            vv=float(value)
            if np.isfinite(vv) and vv>=0:data[j]=vv
    return data

def worker(source,receipt,coast,shard,out):
    out=Path(out);require(not out.exists(),"Refuse overwrite year-locked original raster source")
    all_source=load_original_global(source,receipt)
    selected=all_source.loc[all_source.year.isin(years(shard))].copy()
    require(len(selected)>0 and set(selected.year)==set(years(shard)),
            "Original GBIF source year shard has missing year")
    selected=same_island_photo_grid_centers(selected,coast)
    selected["prior_PCNl_value"]=np.nan
    m=manifest_read(MANIFEST)["files"]
    logs=[]
    out.mkdir(parents=True)
    for y in years(shard):
        raster_year=y-1
        file=f"PCNL{raster_year}.tif"
        r=fetch_one(file,m[file],out/file)
        eligible=selected.year.eq(y)&selected.center_inside_same_original_GSHHG_island
        d=selected.loc[eligible]
        with rasterio.open(out/file) as src:
            vals=genuine_values(src,d)
            original_geometry={"CRS":str(src.crs),
                               "original_res":list(src.res),
                               "nodata":src.nodata,
                               "original_band_scale":list(src.scales)}
        selected.loc[eligible,"prior_PCNl_value"]=vals
        logs.append({
            "original_2026GBIF_observation_year":y,
            "actual_prior_PCNL_raster_year":raster_year,
            "original_source_year_cell_rows":int(selected.year.eq(y).sum()),
            "original_source_year_cells_center_on_same_island":int(eligible.sum()),
            "original_source_year_cells_raster_numeric_success":int(np.isfinite(vals).sum()),
            "original_source_year_cells_original_center_off_island":int(
                 (selected.year.eq(y)&~selected.center_inside_same_original_GSHHG_island).sum()),
            "source_PCNl_original_MD5":m[file],
            "source_raster_original_receipt":r,
            "original_PCNl_raster_geometry":original_geometry})
        print("ALAN_GLOBAL_GBIF_TRUE_PRIOR_YEAR_PCNL",
              json.dumps(logs[-1],sort_keys=True),flush=True)
        (out/file).unlink()
    selected["PCNL_year"]=selected.year-1
    selected["prior_log1p_PCNL"]=np.log1p(selected.prior_PCNl_value)
    target=out/f"global_GBIF_island_fourfamily_prior_PCNL_shard_{shard}.csv.gz"
    selected.to_csv(target,index=False,compression="gzip")
    audit={
        "schema":SCHEMA,"worker_shard":shard,
        "original_source_GBIF_2026_08_SHA256":SOURCE_SHA,
        "original_fourfamily_GBIF_years":years(shard),
        "real_source_cell_year_rows":len(selected),
        "original_yearly_raster_coverage":logs,
        "site_GSHHG_center_inside_same_island_rows":int(
              selected.center_inside_same_original_GSHHG_island.sum()),
        "source_valid_prior_PCNL_cell_years":int(selected.prior_log1p_PCNL.notna().sum()),
        "source_year_shard_csv_SHA256":sha256(target),
        "NO_future_light_for_old_survey_year":True,
        "no_raw_moth_photo_coordinates_republished":True,
        "original_record_source_not_moth_abundance_or_pollination":True}
    (out/f"global_GBIF_island_fourfamily_prior_PCNL_shard_{shard}.json").write_text(
       json.dumps(audit,sort_keys=True,indent=2)+"\n")
    return audit

def collect(source,receipt,parts,out):
    out=Path(out);require(not out.exists(),"Refuse overwrite signed global PCNL source")
    old=load_original_global(source,receipt)
    inputs=[];seen=[]
    for shard in range(SHARDS):
        p=Path(parts)/f"global_GBIF_island_fourfamily_prior_PCNL_shard_{shard}.csv.gz"
        q=Path(parts)/f"global_GBIF_island_fourfamily_prior_PCNL_shard_{shard}.json"
        a=json.loads(q.read_text())
        require(a["schema"]==SCHEMA and
                a["worker_shard"]==shard and
                a["original_fourfamily_GBIF_years"]==years(shard) and
                a["original_source_GBIF_2026_08_SHA256"]==SOURCE_SHA and
                sha256(p)==a["source_year_shard_csv_SHA256"],
                "Historical source worker missed original MD5 or source SHA")
        v=pd.read_csv(p,dtype={"island":str})
        require(len(v)==a["real_source_cell_year_rows"] and
                (v.PCNL_year==v.year-1).all(),
                "Past-year light leaks or original cell-year source rows missing")
        inputs.append(v)
        seen.extend(a["original_fourfamily_GBIF_years"])
    require(sorted(seen)==OBS_YEARS,
            "Not all 13 original moth-source observation years joined with real PCNL")
    full=pd.concat(inputs,ignore_index=True)
    keys=["island","year","latbin","lonbin"]
    require(len(full)==NROWS and not full.duplicated(keys).any() and
            int(full.all_four_family_records.sum())==NOBS,
            "Original 22.7m insect observation rows changed in raster join")
    joined=old.merge(full[keys+["prior_PCNl_value","prior_log1p_PCNL",
                               "PCNL_year","center_inside_same_original_GSHHG_island",
                               "area_km2"]],
                     on=keys,how="left",validate="one_to_one",indicator=True)
    require(joined._merge.eq("both").all() and (joined.PCNL_year==joined.year-1).all(),
            "Exact observed cell-year to authentic prior-year raster join failed")
    out.mkdir(parents=True)
    csv=out/"ALAN_worldwide_GBIF_228361_island_cell_years_original_prior_PCNL_2012_2024.csv.gz"
    joined.to_csv(csv,index=False,compression="gzip")
    valid=joined.loc[joined.prior_log1p_PCNL.notna()].copy()
    supported=valid.loc[valid.moth_family_records.gt(0) &
                        valid.butterfly_family_records.gt(0) &
                        valid.all_four_family_records.ge(5)]
    study_groups=supported.groupby(["island","year"]).agg(
       n_distinct_sites=("lonbin","size"),
       distinct_light=("prior_log1p_PCNL","nunique"),
       light_min=("prior_log1p_PCNL","min"),
       light_max=("prior_log1p_PCNL","max"))
    spatial=study_groups.loc[(study_groups.n_distinct_sites>=2)&
          (study_groups.distinct_light>=2)&
          ((study_groups.light_max-study_groups.light_min)>=.02)]
    d={
       "schema":SCHEMA,"source_GBIF_original_CSV_SHA256":SOURCE_SHA,
       "true_global_original_source_cells":NROWS,
       "true_global_fourfamily_original_human_observation_records":NOBS,
       "original_source_islands":int(joined.island.nunique()),
       "true_13_observation_years":OBS_YEARS,
       "all_prior_raster_years":[y-1 for y in OBS_YEARS],
       "source_GSHHG_center_inside_same_island_cells":int(
           joined.center_inside_same_original_GSHHG_island.astype(str).str.lower().eq("true").sum()),
       "real_complete_prioryear_PCNL_sites":len(valid),
       "real_prioryear_PCNL_supported_islands":int(valid.island.nunique()),
       "original_source_both_guild_5plus_record_cells":70219,
       "real_PCNL_supported_5plus_both_guild_cells":len(supported),
       "real_PCNL_supported_5plus_both_guild_islands":int(supported.island.nunique()),
       "same_island_year_2plus_distinct_PCNl_cells_span0p02_groups":len(spatial),
       "same_island_year_2plus_distinct_PCNl_cells_span0p02_islands":int(
           spatial.reset_index().island.nunique()),
       "real_complete_joined_global_source_CSV_SHA256":sha256(csv),
       "center_on_island_is_strict_0p05_pixel_center_proxy_not_original_GBIF_point":True,
       "GBIF_insect_counts_not_population_density_or_true_moth_absence":True,
       "only_pre_observation_PCNL_images_not_sameyear_or_future":True,
       "zero_trap_nights_and_actual_observer_effort_not_provided":True,
       "no_insect_plant_pollen_or_population_effect_estimated":True}
    (out/"ALAN_2026_11_global_fourfamily_prior_PCNL_real_source_receipt.json").write_text(
       json.dumps(d,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_WORLDWIDE_GBIF_PRIOR_PCNL_FULL_JOIN",
          json.dumps(d,sort_keys=True),flush=True)
    return d

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="command",required=True)
    w=sub.add_parser("worker")
    w.add_argument("--source",required=True,type=Path)
    w.add_argument("--source-receipt",required=True,type=Path)
    w.add_argument("--coast",required=True,type=Path)
    w.add_argument("--shard",required=True,type=int)
    w.add_argument("--out",required=True,type=Path)
    c=sub.add_parser("collect")
    c.add_argument("--source",required=True,type=Path)
    c.add_argument("--source-receipt",required=True,type=Path)
    c.add_argument("--parts",required=True,type=Path)
    c.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    if a.command=="worker":
        worker(a.source,a.source_receipt,a.coast,a.shard,a.out)
    else:
        collect(a.source,a.source_receipt,a.parts,a.out)

if __name__=="__main__":
    main()
