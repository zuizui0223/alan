#!/usr/bin/env python3
"""Gecko-vs-skink source PHOTO records × REAL PREYEAR global PCNL, not occupancy.

Actual GBIF 2026-08 GSHHG source audited 66,038 gecko and 35,361
Scincidae HUMAN_OBSERVATION source records 2013-2025. This independently
retains geographic sampled cells with gecko OR skink source record, with
the exact 2026-08 source release MD5 SHA256 plus 13 genuine annually
MD5-verified source PCNL rasters (year minus one). Zero in a category means
'not documented among the audited source photo records', NOT species absent.

Restricted original 0.05° cell *centres* must be on the SAME authentic
GSHHG island land, not ocean. No 2024 light assigned to earlier visits.
No extrapolated island occupancy, actual gecko consumption of moths,
effective pollination or true species abundance.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from scripts.analysis.alan_global_GBIF_true4family_priorPCNL_20261110 import (
    same_island_photo_grid_centers,genuine_values)
from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one,manifest_read)
from scripts.analysis.alan_global_GBIF_gecko_island_original_stream_20261110 import (
    GECKO_FAMILIES,SCHEMA as GECKO_SCHEMA,sha256)

SCHEMA="ALAN_true_202608_gecko_v_scink_photo_cellyear_priorPCNL_year12to24_v1"
ORIGINAL_GECKO_SHA="c2e0c823ee9fcf61182b920724e488d9130b1d1ca9d53266511c2da8873983c9"
SHARDS=4
YEARS=tuple(range(2013,2026))

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_GECKO_SOURCE_ANNUAL_PCNL: "+str(msg))

def shard_years(shard):
    require(shard in range(SHARDS),"Unknown source year shard")
    return list(YEARS[shard::SHARDS])

def original_reptile_source(path,receipt):
    d=json.loads(Path(receipt).read_text())
    require(d["schema"]==GECKO_SCHEMA and
            d["full_original_Squamata_all_year_and_basis_rows"]==11634053 and
            d["source_gecko_family_cell_year_CSV_SHA256"]==ORIGINAL_GECKO_SHA and
            sha256(path)==ORIGINAL_GECKO_SHA and
            d["source_gecko_photo_observation_records"]==66038 and
            d["source_Scincidae_photo_observation_records"]==35361,
            "Not original authentic 101399 reptile GBIF source photo occurrences")
    z=pd.read_csv(path,dtype={"island":str})
    k=["island","year","latbin","lonbin"]
    require(not z.duplicated(k+["actual_source_taxon_family"]).any() and
            z.year.between(2013,2025).all() and
            set(z.actual_source_taxon_family).issubset(
                set(GECKO_FAMILIES)|{"Scincidae"}) and
            int(z.actual_GBIF_human_observation_records.sum())==101399,
            "Original reptile photo taxonomic source rows cannot be reconciled")
    z["gecko_source_photo_count"]=np.where(
        z.actual_source_taxon_family.isin(GECKO_FAMILIES),
        z.actual_GBIF_human_observation_records,0)
    z["scink_source_photo_count"]=np.where(
        z.actual_source_taxon_family.eq("Scincidae"),
        z.actual_GBIF_human_observation_records,0)
    all_cells=z.groupby(k,sort=True)[[
        "gecko_source_photo_count","scink_source_photo_count"]].sum().reset_index()
    require(not all_cells.duplicated(k).any() and
            int(all_cells.gecko_source_photo_count.sum())==66038 and
            int(all_cells.scink_source_photo_count.sum())==35361,
            "Original source taxon photo counts lost during true site-year aggregation")
    all_cells["original_total_reptile_photo_source_events"]=(
        all_cells.gecko_source_photo_count+all_cells.scink_source_photo_count)
    all_cells["original_gecko_photo_fraction"]=(
        all_cells.gecko_source_photo_count/
        all_cells.original_total_reptile_photo_source_events)
    all_cells["longitude_center"]=(all_cells.lonbin.astype(float)+.5)*.05-180.
    all_cells["latitude_center"]=(all_cells.latbin.astype(float)+.5)*.05-90.
    return all_cells

def worker(source,receipt,coast,shard,out):
    out=Path(out)
    require(not out.exists(),"Never overwrite source original year-matched light")
    d=original_reptile_source(source,receipt)
    selected=d.loc[d.year.isin(shard_years(shard))].copy()
    require(len(selected)>0 and set(selected.year)==set(shard_years(shard)),
            "Actual source year photos missing from shard")
    selected=same_island_photo_grid_centers(selected,coast)
    require(selected.center_inside_same_original_GSHHG_island.notna().all(),
            "Source island land centre membership missing")
    selected["prior_PCNL_native"]=np.nan
    source_md5=manifest_read("config/alan_pcnl_1992_2024_zenodo_v3_md5.json")["files"]
    out.mkdir(parents=True)
    rec=[]
    for year in shard_years(shard):
        file=f"PCNL{year-1}.tif"
        provenance=fetch_one(file,source_md5[file],out/file)
        valid=selected.year.eq(year)&selected.center_inside_same_original_GSHHG_island
        part=selected.loc[valid]
        with rasterio.open(out/file) as src:
            samples=genuine_values(src,part)
            geom={"native_CRS":str(src.crs),"native_pixel_size_degrees":list(src.res),
                  "original_nodata":src.nodata}
        selected.loc[valid,"prior_PCNL_native"]=samples
        rec.append({
           "source_year_observed_GBIF_reptile_photo":year,
           "actual_PCNL_raster_year":year-1,
           "source_reptile_photo_cells_observed":int(selected.year.eq(year).sum()),
           "original_photo_source_cells_center_on_island_land":int(valid.sum()),
           "original_PCNL_raster_valid_numeric_cells":int(np.isfinite(samples).sum()),
           "original_raster_MD5_provenance":provenance,
           "native_raster_geometry":geom})
        print("ALAN_TRUE_GECKO_SCINK_PRIOR_YEAR_PCNL",json.dumps({
            "year":year,"island_photo_site_years":rec[-1]["source_reptile_photo_cells_observed"],
            "original_PCNL_valid":rec[-1]["original_PCNL_raster_valid_numeric_cells"]},
            sort_keys=True),flush=True)
        (out/file).unlink()
    selected["prior_PCNL_year"]=selected.year-1
    selected["log1p_prior_PCNL"]=np.log1p(selected.prior_PCNL_native)
    selected.sort_values(["island","year","latbin","lonbin"],inplace=True)
    csv=out/f"gecko_scink_prior_PCNL_original_reptile_shard{shard}.csv.gz"
    selected.to_csv(csv,index=False,compression="gzip")
    audit={
      "schema":SCHEMA,
      "shard":shard,
      "original_reptile_source_snapshot_SHA256":ORIGINAL_GECKO_SHA,
      "source_observed_years":shard_years(shard),
      "true_original_reptile_site_year_rows":len(selected),
      "actual_source_PCNL_valid_cells":int(selected.log1p_prior_PCNL.notna().sum()),
      "source_observed_GSHHG_islands":int(selected.island.nunique()),
      "source_all_reptile_photo_events":int(
          selected.original_total_reptile_photo_source_events.sum()),
      "source_Gecko_photo_events":int(selected.gecko_source_photo_count.sum()),
      "source_Scink_photo_events":int(selected.scink_source_photo_count.sum()),
      "original_raster_source_MD5_year_receipts":rec,
      "original_4way_source_CSV_SHA256":sha256(csv),
      "photodocumentation_not_true_reptile_density_or_invasive_presence":True,
      "every_satellite_source_year_precedes_observed_reptile_photo":True}
    (out/f"gecko_scink_prior_PCNL_original_reptile_shard{shard}.json").write_text(
         json.dumps(audit,sort_keys=True,indent=2)+"\n")
    return audit

def collect(source,receipt,parts,out):
    out=Path(out);require(not out.exists(),"Do not overwrite original source model")
    base=original_reptile_source(source,receipt)
    all_shards=[];years=[]
    for shard in range(SHARDS):
        file=Path(parts)/f"gecko_scink_prior_PCNL_original_reptile_shard{shard}.csv.gz"
        ledger=Path(parts)/f"gecko_scink_prior_PCNL_original_reptile_shard{shard}.json"
        a=json.loads(ledger.read_text())
        require(a["schema"]==SCHEMA and a["shard"]==shard and
                a["original_reptile_source_snapshot_SHA256"]==ORIGINAL_GECKO_SHA and
                a["source_observed_years"]==shard_years(shard) and
                sha256(file)==a["original_4way_source_CSV_SHA256"],
                "Missing or altered one original annual nightlight/reptile source")
        frame=pd.read_csv(file,dtype={"island":str})
        require(not frame.duplicated(["island","year","latbin","lonbin"]).any() and
                len(frame)==a["true_original_reptile_site_year_rows"] and
                (frame.prior_PCNL_year==frame.year-1).all(),
                "Original observed photo cell/year not unique or future pixel assigned")
        all_shards.append(frame);years.extend(shard_years(shard))
    require(sorted(years)==list(YEARS),
            "All 13 original GBIF reptile observed years not source-checked")
    data=pd.concat(all_shards,ignore_index=True)
    k=["island","year","latbin","lonbin"]
    require(len(data)==len(base) and not data.duplicated(k).any() and
            int(data.gecko_source_photo_count.sum())==66038 and
            int(data.scink_source_photo_count.sum())==35361,
            "Original reptile source counts lost before historical PCNL")
    merged=base.merge(data[k+["prior_PCNL_native","log1p_prior_PCNL",
                              "prior_PCNL_year","center_inside_same_original_GSHHG_island",
                              "area_km2"]],
                on=k,how="left",validate="one_to_one",indicator=True)
    require(merged._merge.eq("both").all() and
            merged.log1p_prior_PCNL.notna().sum()>0,
            "Source gecko photo site-year or geographic centre lost")
    out.mkdir(parents=True)
    outfile=out/"gecko_scink_true_2013_2025_original_island_photo_cells_priorPCNL.csv.gz"
    merged.to_csv(outfile,index=False,compression="gzip")
    valid=merged.loc[merged.log1p_prior_PCNL.notna()].copy()
    within=valid.groupby(["island","year"]).agg(
        actual_photo_cells=("lonbin","size"),
        actual_photo_light_values=("log1p_prior_PCNL","nunique"),
        light_min=("log1p_prior_PCNL","min"),
        light_max=("log1p_prior_PCNL","max"),
        gecko_photos=("gecko_source_photo_count","sum"),
        scink_photos=("scink_source_photo_count","sum"))
    supported=within.loc[(within.actual_photo_cells>=2)&
                          (within.actual_photo_light_values>=2)&
                          (within.light_max-within.light_min>=.02)&
                          (within.gecko_photos>=1)&(within.scink_photos>=1)]
    support2=within.loc[(within.actual_photo_cells>=2)&
                          (within.actual_photo_light_values>=2)&
                          (within.light_max-within.light_min>=.02)&
                          (within.gecko_photos>=5)&(within.scink_photos>=5)]
    audit={
       "schema":SCHEMA,
       "original_GBIF_snap_DOI":"10.15468/dl.8yrbe7",
       "historically_prior_PCNL_original_DOI":"10.5281/zenodo.17013414",
       "all_original_gecko_and_scink_photo_record_events":101399,
       "original_Gekkonidae_plus_other5photo_event_count":66038,
       "original_Scincidae_photo_event_count":35361,
       "real_source_0p05_reptile_photo_cell_years":len(merged),
       "original_record_bearing_GSHHG_islands":int(merged.island.nunique()),
       "real_source_valid_prior_year_PCNL_site_years":len(valid),
       "real_source_valid_prior_year_PCNL_islands":int(valid.island.nunique()),
       "same_island_year_source_dual_gecko_scink_min1_light_groups":len(supported),
       "same_island_year_source_dual_gecko_scink_min1_independent_islands":
             int(supported.reset_index().island.nunique()),
       "same_island_year_source_dual_gecko_scink_min5_light_groups":len(support2),
       "same_island_year_source_dual_gecko_scink_min5_independent_islands":
             int(support2.reset_index().island.nunique()),
       "source_prior_year_lighting_not_physical_lux_at_gecko":True,
       "original_reptile_photo_cell_CSV_SHA256":sha256(outfile),
       "source_zeros_mean_no_photo_NOT_biological_absence":True,
       "no_gecko_moth_feeding_or_population_impact_measured":True}
    (out/"gecko_scink_original_13prior_PCNL_source_receipt.json").write_text(
        json.dumps(audit,sort_keys=True,indent=2)+"\n")
    print("ALAN_GECKO_SCINK_GENUINE_PRIOR_PCNL_GLOBAL_SITE_SUPPORT",
          json.dumps(audit,sort_keys=True),flush=True)
    return audit

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="command",required=True)
    w=sub.add_parser("worker")
    for name in ("source","receipt","coast","out"):
        w.add_argument("--"+name,required=True,type=Path)
    w.add_argument("--shard",required=True,type=int)
    c=sub.add_parser("collect")
    for name in ("source","receipt","parts","out"):
        c.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    if a.command=="worker":worker(a.source,a.receipt,a.coast,a.shard,a.out)
    else:collect(a.source,a.receipt,a.parts,a.out)

if __name__=="__main__":main()
