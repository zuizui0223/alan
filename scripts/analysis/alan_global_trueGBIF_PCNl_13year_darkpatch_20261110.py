#!/usr/bin/env python3
"""True annual 2012-2024 PCNL island dark mosaic for 2013-2025 GBIF moth sites.

Only exact original verified 62,688 co-observed moth+butterfly photo source
site-years on 196 GSHHG islands supporting within-island same-year light
contrasts. Measure genuine pre-observation annual ~1km PCNL dark fractions
and four-neighbour patch shape at 5km geodesic radius, original island land
minus fresh waters. No images, future light, source absence imputation, pollen
delivery, moth density or causal light-effect claims.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely import union_all

from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land
from scripts.analysis.alan_dark_refugia_local_pcnl_20261110 import dark_refuge_metrics
from scripts.analysis.alan_global_GBIF_moth_record_fraction_within_island_FE_20261110 import (
    source as signed_source,frame)
from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one,manifest_read)
SCHEMA="ALAN_true_62688_global_moth_photo_cells_2012_2024_PCNl_5km_mosaic_v1"
YEARS=tuple(range(2013,2026))
SHARDS=4
THRESHOLDS=(0.,3.)
EXPECTED_SITES=62688
EXPECTED_ISLANDS=196
SOURCE_SHA="64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_GLOBAL_DARKPATCH: "+str(msg))
def digest(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()
def years(shard):
    require(shard in range(SHARDS),"Unregistered actual source year shard")
    return list(YEARS[shard::SHARDS])

def source_eligible(source,receipt):
    require(digest(source)==SOURCE_SHA,
            "Historical 13year original source GBIF/PCNL SHA changed")
    d,info=signed_source(source,receipt)
    x,support=frame(d,"recorded_moth_family_share_pp",
                    include_zero_category=False,obs_threshold=5)
    require(len(x)==EXPECTED_SITES and x.island.nunique()==EXPECTED_ISLANDS and
            x.island_year.nunique()==1083 and
            x.prior_log1p_PCNL.notna().all() and
            x.moth_family_records.gt(0).all() and
            x.butterfly_family_records.gt(0).all() and
            x.all_four_family_records.ge(5).all(),
            "Genuine global 196-island within-year photo-source cohort drift")
    require(x.center_inside_same_original_GSHHG_island.astype(str).str.lower().eq("true").all(),
            "Original source photo-cell center is not on original GSHHG island land")
    return x

def worker(source,receipt,coast,shard,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwrite real 5km dark-patch source")
    x=source_eligible(source,receipt)
    years_this=years(shard)
    chosen=x.loc[x.year.isin(years_this)].copy()
    require(chosen.year.nunique()==len(years_this) and
            not chosen.duplicated(["island","year","latbin","lonbin"]).any(),
            "Original source photos duplicated or year subset not complete")
    land=load_signed_land(coast).set_index("land_id",drop=False)
    require(set(chosen.island).issubset(set(land.index)),
            "Original source island not in SHA-pinned GSHHG")
    lakes=gpd.read_file(Path(coast)/"GSHHS_h_L2.shp")
    require(lakes.crs is not None,"Original GSHHG inland lakes have no source CRS")
    lakes=lakes.to_crs(4326)
    li=lakes.sindex
    polygons={}
    for island in chosen.island.unique():
        v=land.loc[island,"geometry"]
        nearest=li.query(v,predicate="intersects")
        if len(nearest):
            w=[q for q in lakes.geometry.iloc[nearest] if q.intersects(v)]
            if w:v=v.difference(union_all(w))
        polygons[island]=v if not v.is_empty and v.is_valid else None
    out.mkdir(parents=True)
    rows=[];raster_audit=[];errors=[]
    published=manifest_read("config/alan_pcnl_1992_2024_zenodo_v3_md5.json")["files"]
    for observed_year in years_this:
        raster_year=observed_year-1
        raster_name=f"PCNL{raster_year}.tif"
        status=fetch_one(raster_name,published[raster_name],out/raster_name)
        group=chosen.loc[chosen.year.eq(observed_year)]
        with rasterio.open(out/raster_name) as src:
            require(src.count==1 and src.crs is not None,
                    "Original PCNL source raster source invalid")
            for island,g in group.groupby("island",sort=True):
                poly=polygons[island]
                for row in g.itertuples(index=False):
                    src_id={"island":str(island),
                            "year":int(observed_year),
                            "latbin":int(row.latbin),
                            "lonbin":int(row.lonbin),
                            "original_moth_family_record_count":int(row.moth_family_records),
                            "original_butterfly_family_record_count":int(row.butterfly_family_records),
                            "original_total_fourfamily_source_records":int(row.all_four_family_records),
                            "prior_year_PCNL_center_log1p":float(row.prior_log1p_PCNL)}
                    if poly is None:
                        for t in THRESHOLDS:
                            rows.append({**src_id,"threshold_PCNL_native":t,
                                         "status":"invalid_original_island_after_lake_mask"})
                        continue
                    try:
                        metrics=dark_refuge_metrics(
                           src,poly,float(row.longitude_center),float(row.latitude_center),
                           radius_km=5.,thresholds=THRESHOLDS,
                           min_valid_land_pixels=16,min_coverage=.8)
                    except ValueError as e:
                        errors.append({"island":str(island),"year":observed_year,
                                       "latbin":int(row.latbin),"lonbin":int(row.lonbin),
                                       "reason":str(e)})
                        for t in THRESHOLDS:
                            rows.append({**src_id,"threshold_PCNL_native":t,
                                         "status":"source_geometry_unsupported"})
                        continue
                    for m in metrics:
                        rows.append({**src_id,**m})
        raster_audit.append({"source_raster_original":raster_name,
                             "actual_photo_observation_year":observed_year,
                             "original_MD5":published[raster_name],
                             "source_raster_receipt":status,
                             "actual_original_cog_photo_cell_years":len(group)})
        print("ALAN_WORLDWIDE_GBIF_DARK_PATCH_GENUINE_YEAR",
              json.dumps({"year":observed_year,"source_pcNl_year":raster_year,
                          "source_original_photo_cells":len(group),
                          "accumulated_source_photo_cells_2thresholds":len(rows)},
                          sort_keys=True),flush=True)
        (out/raster_name).unlink()
    data=pd.DataFrame(rows)
    require(len(data)==2*len(chosen) and
            not data.duplicated(["island","year","latbin","lonbin",
                                 "threshold_PCNL_native"]).any(),
            "True original global moth source site-year dark-patch cells duplicated")
    target=out/f"alan_true_worldwide_GBIF_PCNl_darkpatch_yearshard_{shard}.csv.gz"
    data.to_csv(target,index=False,compression="gzip")
    g=data.loc[data.threshold_PCNL_native.eq(0.)]
    audit={"schema":SCHEMA,"shard":shard,
           "source_actual_photo_observation_years":years_this,
           "original_genuine_62688_global_GBIF_source_SHA256":SOURCE_SHA,
           "source_site_years_processed":len(chosen),
           "source_5km_landpatch_metric_success":int(
                g.status.eq("observed_satellite_geometry_only").sum()),
           "source_5km_landpatch_metric_fail_or_unobserved":int(
                (~g.status.eq("observed_satellite_geometry_only")).sum()),
           "source_5km_landpatch_unique_true_islands":int(chosen.island.nunique()),
           "source_2threshold_original_darkpatch_csv_SHA256":digest(target),
           "source_year_MD5_receipts":raster_audit,
           "source_geometry_failure_examples":errors[:20],
           "no_500m_2024_future_light_backprojection":True,
           "zeros_in_radiance_are_NOT_true_insect_absence":True,
           "real_moth_photo_source_ratio_NOT_moth_density":True}
    (out/f"alan_true_worldwide_GBIF_PCNl_darkpatch_yearshard_{shard}.json").write_text(
          json.dumps(audit,sort_keys=True,indent=2)+"\n")
    return audit

def collect(source,receipt,parts,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite globally signed source patch table")
    original=source_eligible(source,receipt)
    frames=[];seen=[]
    for shard in range(SHARDS):
        p=Path(parts)/f"alan_true_worldwide_GBIF_PCNl_darkpatch_yearshard_{shard}.csv.gz"
        q=Path(parts)/f"alan_true_worldwide_GBIF_PCNl_darkpatch_yearshard_{shard}.json"
        metadata=json.loads(q.read_text())
        require(metadata["schema"]==SCHEMA and
                metadata["shard"]==shard and
                metadata["source_actual_photo_observation_years"]==years(shard) and
                metadata["original_genuine_62688_global_GBIF_source_SHA256"]==SOURCE_SHA and
                digest(p)==metadata["source_2threshold_original_darkpatch_csv_SHA256"],
                "One original satellite source patch shard was not successful")
        d=pd.read_csv(p,dtype={"island":str})
        require(len(d)==2*metadata["source_site_years_processed"] and
                set(d.threshold_PCNL_native)==set(THRESHOLDS),
                "Original source raster-year threshold/source identity drift")
        frames.append(d);seen.extend(metadata["source_actual_photo_observation_years"])
    require(sorted(seen)==list(YEARS),"Original 13 observation years not complete")
    all_exposure=pd.concat(frames,ignore_index=True)
    keys=["island","year","latbin","lonbin","threshold_PCNL_native"]
    require(not all_exposure.duplicated(keys).any() and
            len(all_exposure)==2*EXPECTED_SITES,
            "Original 62688×2 true satellite source patch rows do not reconcile")
    expected=original[["island","year","latbin","lonbin"]]
    unique=all_exposure.loc[all_exposure.threshold_PCNL_native.eq(0.)]
    check=expected.merge(unique[["island","year","latbin","lonbin"]],
                         on=["island","year","latbin","lonbin"],
                         how="left",indicator=True,validate="one_to_one")
    require(check._merge.eq("both").all(),
            "One original GSHHG photo-grid-year source cell was omitted")
    out.mkdir(parents=True)
    target=out/"alan_real_global_62688_GBIF_moth_source_annual_prior_PCNl_5km_darkpatch.csv.gz"
    all_exposure.to_csv(target,index=False,compression="gzip")
    support={}
    for t in THRESHOLDS:
        x=all_exposure.loc[all_exposure.threshold_PCNL_native.eq(t)].copy()
        good=x.loc[x.status.eq("observed_satellite_geometry_only")]
        known=good.loc[good.largest_dark_component_fraction_of_dark.notna()]
        support[str(t)]={
            "genuine_source_valid_photo_grid_years":len(good),
            "genuine_source_valid_islands":int(good.island.nunique()),
            "actual_dark_region_present_photo_grid_years":len(known),
            "median_original_land_dark_fraction":
                float(good.dark_land_fraction.median()) if len(good) else None,
            "median_largest_connected_dark_fraction":
                float(known.largest_dark_component_fraction_of_dark.median()) if len(known) else None,
            "genuine_source_no_dark_pixels_detected":int(good.n_dark_raster_pixels.eq(0).sum()),
            "photo_grid_years_uncovered":int(len(x)-len(good))}
    audit={
        "schema":SCHEMA,
        "source_actual_global_GBIF_photo_cells":EXPECTED_SITES,
        "source_original_islands":EXPECTED_ISLANDS,
        "year_matched_corrected_PCNL_original_rasters_2012_2024":True,
        "source_2threshold_global_patch_CSV_SHA256":digest(target),
        "real_prior_year_5km_darkpatch_support":support,
        "sea_pixels_and_inland_water_masked":True,
        "no_insect_density_pollen_deposition_or_effective_moth_visitation":True,
        "source_photographic_GBIF_record_share_not_population_response":True,
        "2thresholds_same_photo_site_year_NOT_independent_replicates":True}
    (out/"alan_actual_global_GBIF_moth_darkpatch_all13years_receipt.json").write_text(
         json.dumps(audit,sort_keys=True,indent=2)+"\n")
    print("ALAN_GLOBAL_TRUE_13YEAR_DARK_PATCH_SITE_MATCH",
          json.dumps(audit,sort_keys=True),flush=True)
    return audit

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="mode",required=True)
    w=sub.add_parser("worker")
    w.add_argument("--source",type=Path,required=True)
    w.add_argument("--receipt",type=Path,required=True)
    w.add_argument("--coast",type=Path,required=True)
    w.add_argument("--shard",type=int,required=True)
    w.add_argument("--out",type=Path,required=True)
    c=sub.add_parser("collect")
    c.add_argument("--source",type=Path,required=True)
    c.add_argument("--receipt",type=Path,required=True)
    c.add_argument("--parts",type=Path,required=True)
    c.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.mode=="worker":worker(a.source,a.receipt,a.coast,a.shard,a.out)
    else:collect(a.source,a.receipt,a.parts,a.out)

if __name__=="__main__":
    main()
