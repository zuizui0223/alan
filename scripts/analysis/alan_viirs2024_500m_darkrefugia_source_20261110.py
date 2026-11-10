#!/usr/bin/env python3
"""True ~500m 2024 VIIRS vs ~1km PCNL 2024 dark-mosaic resolution audit.

Only original 2025 photographed plant site coordinates, frozen original
2024 PCNL source exposure and original SHA-pinned GSHHG polygons are used.
This is NOT moth abundance, nocturnal pollination, light at a flower or a
causal assessment of night light. A source change is confounded with pixel
resolution: we cannot attribute every difference to resolution alone.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely import union_all

from scripts.analysis.alan_dark_refugia_local_pcnl_20261110 import (
    dark_refuge_metrics, sha256 as source_sha256)
from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land

DATASET = "https://zenodo.org/records/17294744"
RASTER_NAME = "nightlights.average_viirs.v21_m_500m_s_20240101_20241231_go_epsg4326_v20250904.tif"
RASTER_MD5 = "0b976892d40c2d0d99177c727289eaf2"
SOURCE_SHA = "32327abb37ee85f3d1d9e7446d9d742fb5adf1457dab60fae1aae791c29f7837"
THRESHOLDS_NW = (0., 1., 3.)
# The source publisher states VIIRS nW/cm2/sr was multiplied by 10 for COG storage.
PUBLISHED_MULTIPLIER = 10.
SCHEMA = "ALAN_original_2025_photo_sites_VIIRS2024_500m_vs_PCNL2024_1km_dark_refugia_v1"

def require(ok, why):
    if not bool(ok): raise ValueError("ALAN_VIIRS500_SOURCE: " + str(why))

def file_md5(path):
    h=hashlib.md5()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):
            h.update(b)
    return h.hexdigest()

def valid_old_source(source, receipt):
    info=json.loads(Path(receipt).read_text())
    require(info["schema"]=="alan_dark_refugia_original_2025_sites_prior_2024_PCNL_GSHHG_v1",
            "Wrong original satellite cohort")
    require(info["source_output_csv_sha256"]==SOURCE_SHA and
            source_sha256(source)==SOURCE_SHA, "Original 2025 site PCNL archive changed")
    old=pd.read_csv(source,dtype={"island":str},low_memory=False)
    keys=["island","latbin","lonbin","photo_observation_year"]
    need=set(keys+["threshold_PCNL_native","status","lat","lon","dark_land_fraction",
                  "largest_dark_component_fraction_of_dark","local_mean_PCNL_native"])
    require(need.issubset(old.columns),"Missing original real PCNL columns")
    require(set(old.threshold_PCNL_native)=={0.,1.,3.} and
            old.photo_observation_year.eq(2025).all(), "Wrong source years/thresholds")
    require(not old.duplicated(keys+["threshold_PCNL_native"]).any(),
            "Repeated original photo site/threshold")
    photo_sites=old.loc[old.threshold_PCNL_native.eq(0.)].copy()
    require(len(photo_sites)==4728 and
            photo_sites.island.nunique()==1758 and
            photo_sites.status.eq("observed_satellite_geometry_only").sum()==3543,
            "Frozen genuine 4728/1758 and 3543 site support changed")
    return photo_sites,old,info

def quantiles(x):
    x=pd.to_numeric(x,errors="coerce").dropna()
    return {str(p):float(x.quantile(p)) for p in (.05,.25,.50,.75,.95)} if len(x) else {}

def run(coarse,coarse_receipt,coast,raster,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwrite previously measured 500m satellite record")
    photo_sites,old,old_receipt=valid_old_source(coarse,coarse_receipt)
    require(Path(raster).name==RASTER_NAME and file_md5(raster)==RASTER_MD5,
            "Actual published Zenodo 2024 COG filename or MD5 differs")
    land=load_signed_land(coast).set_index("land_id",drop=False)
    require(photo_sites.island.isin(land.index).all() and
            land.loc[photo_sites.island,"land_class"].eq("island").all(),
            "Original GSHHG land/photograph identity mismatches")
    lakes=gpd.read_file(Path(coast)/"GSHHS_h_L2.shp")
    require(lakes.crs is not None,"No original GSHHG inland lake geometry CRS")
    lakes=lakes.to_crs(4326)
    lakes_index=lakes.sindex
    result=[]
    failures=[]
    raster_details={}
    with rasterio.open(raster) as src:
        require(src.count==1 and src.crs is not None and src.crs.to_epsg()==4326,
                "VIIRS file is not one EPSG:4326 original-raster band")
        require(0 < abs(src.res[0]) <= .007 and 0 < abs(src.res[1]) <= .007,
                "Not published sub-km VIIRS source pixel size")
        require(len(src.scales)==1 and np.isclose(src.scales[0],1.),
                "GeoTIFF band has unexpected internal scaling; cannot multiply published factor")
        raster_details={
            "height":src.height,"width":src.width,"crs":str(src.crs),
            "native_res_degrees":list(src.res),"nodata":src.nodata,
            "raster_scale":src.scales[0],"raster_offset":src.offsets[0],
            "band_tags":src.tags(1),"file_tags":src.tags(),
            "value_conversion_to_nW_cm2_sr":"original stored numeric raster / 10",
            "file_MD5_verified":RASTER_MD5,
        }
        for island,g in photo_sites.groupby("island",sort=True):
            poly=land.loc[island,"geometry"]
            nearby=lakes_index.query(poly,predicate="intersects")
            if len(nearby):
                relevant=[v for v in lakes.geometry.iloc[nearby] if v.intersects(poly)]
                if relevant: poly=poly.difference(union_all(relevant))
            if poly.is_empty or not poly.is_valid:
                failures.extend({"island":island,"reason":"source_island_invalid_post_lake_subtraction",
                                  "latbin":r.latbin,"lonbin":r.lonbin}
                                for r in g.itertuples(index=False))
                continue
            for row in g.itertuples(index=False):
                try:
                    # Thresholds are defined in original physical VIIRS units.
                    # This product stores values x10; pass native 0/10/30.
                    metrics=dark_refuge_metrics(
                        src,poly,float(row.lon),float(row.lat),radius_km=5.,
                        thresholds=tuple(t*PUBLISHED_MULTIPLIER for t in THRESHOLDS_NW),
                        min_valid_land_pixels=64,min_coverage=.80)
                except ValueError as err:
                    failures.append({"island":island,"latbin":row.latbin,
                                     "lonbin":row.lonbin,"reason":str(err)})
                    continue
                for a,threshold in zip(metrics,THRESHOLDS_NW):
                    b=dict(a)
                    b["local_mean_viirs_nW_cm2_sr"]=(
                        b.pop("local_mean_PCNL_native")/PUBLISHED_MULTIPLIER
                        if "local_mean_PCNL_native" in b else None)
                    b.pop("threshold_PCNL_native")
                    b.update({"island":island,
                              "latbin":row.latbin,"lonbin":row.lonbin,
                              "lat":row.lat,"lon":row.lon,
                              "observation_year":2025,"satellite_year":2024,
                              "VIIRS_dark_threshold_nW_cm2_sr":threshold,
                              "source_product":"OpenGeoHub annual VIIRS VNL V2.2 2024 native ~500m COG",
                              "original_photo_site_PCNL_2024_status":row.status})
                    result.append(b)
    data=pd.DataFrame(result)
    require(len(data)>0,"No real original fine-scale island image processed")
    keys=["island","latbin","lonbin","observation_year",
          "VIIRS_dark_threshold_nW_cm2_sr"]
    require(not data.duplicated(keys).any() and
            set(data.VIIRS_dark_threshold_nW_cm2_sr)==set(THRESHOLDS_NW),
            "500m satellite source repetitions or wrong thresholds")
    out.mkdir(parents=True)
    path=out/"original_2025_photo_site_VIIRS2024_500m_landmasked_dark_refugia.csv.gz"
    data.to_csv(path,index=False,compression="gzip")
    summary={}
    for t in THRESHOLDS_NW:
        x=data.loc[data.VIIRS_dark_threshold_nW_cm2_sr.eq(t)]
        ok=x.loc[x.status.eq("observed_satellite_geometry_only")].copy()
        d={"n_attempted_original_sites":len(x),
           "n_valid_original_500m_sites":len(ok),
           "n_valid_islands":int(ok.island.nunique()),
           "median_dark_land_fraction":float(ok.dark_land_fraction.median()) if len(ok) else None,
           "dark_land_fraction_quantiles":quantiles(ok.dark_land_fraction),
           "median_fraction_of_dark_in_largest_connected_component":
               float(ok.largest_dark_component_fraction_of_dark.median()) if len(ok) else None,
           "conditional_connectivity_quantiles":
               quantiles(ok.largest_dark_component_fraction_of_dark),
           "n_sites_no_detected_dark_pixels":int(ok.n_dark_raster_pixels.eq(0).sum()) if len(ok) else 0,
           "n_new_valid_sites_where_original_PCNL1km_failed":
               int(ok.original_photo_site_PCNL_2024_status.ne("observed_satellite_geometry_only").sum()) if len(ok) else 0}
        if t in (0.,3.) and len(ok):
            former=old.loc[old.threshold_PCNL_native.eq(t),[
                "island","latbin","lonbin","dark_land_fraction",
                "largest_dark_component_fraction_of_dark","status"]]
            a=ok.merge(former,on=["island","latbin","lonbin"],validate="one_to_one",
                       suffixes=("_VIIRS500m","_PCNL1km"))
            both=a.loc[a.status_PCNL1km.eq("observed_satellite_geometry_only")].copy()
            delta=both.dark_land_fraction_VIIRS500m-both.dark_land_fraction_PCNL1km
            d["co_supported_source_photo_sites"]=len(both)
            d["co_supported_islands"]=int(both.island.nunique())
            d["VIIRS_minus_PCNL_dark_fraction_quantiles"]=quantiles(delta)
            d["same_sites_spearman_dark_fraction"]=(
                float(both.dark_land_fraction_VIIRS500m.corr(
                    both.dark_land_fraction_PCNL1km,method="spearman"))
                if both.dark_land_fraction_VIIRS500m.nunique()>1 and
                both.dark_land_fraction_PCNL1km.nunique()>1 else None)
            cont=both.loc[both.largest_dark_component_fraction_of_dark_PCNL1km.notna() &
                          both.largest_dark_component_fraction_of_dark_VIIRS500m.notna()]
            d["both_have_dark_regions"]=len(cont)
            d["n_coarse_largest_dark_0p95_plus_but_fine_below_0p90"]=int((
                cont.largest_dark_component_fraction_of_dark_PCNL1km.ge(.95) &
                cont.largest_dark_component_fraction_of_dark_VIIRS500m.lt(.90)).sum())
        summary[str(t)]=d
    # Threshold 0 vs 1 must be rechecked for true 500m raster; 1km
    # was identical at these settings. Do not infer if source masking differs.
    i0=data.loc[data.VIIRS_dark_threshold_nW_cm2_sr.eq(0.)].set_index(
        ["island","latbin","lonbin"])
    i1=data.loc[data.VIIRS_dark_threshold_nW_cm2_sr.eq(1.)].set_index(
        ["island","latbin","lonbin"])
    common=i0.index.intersection(i1.index)
    common_valid=i0.loc[common].status.eq("observed_satellite_geometry_only") & (
        i1.loc[common].status.eq("observed_satellite_geometry_only"))
    ident=(i0.loc[common,"dark_land_fraction"][common_valid].to_numpy() ==
           i1.loc[common,"dark_land_fraction"][common_valid].to_numpy())
    report={
        "schema":SCHEMA,
        "source_raster_Zenodo_record":DATASET,
        "raster_name":RASTER_NAME,
        "source_raster_MD5":RASTER_MD5,
        "previous_PCNl_original_source_CSV_SHA256":SOURCE_SHA,
        "raster_details":raster_details,
        "source_original_photo_sites_attempted":len(photo_sites),
        "source_islands_attempted":int(photo_sites.island.nunique()),
        "source_original_2024_PCNL_valid_sites":int(
            photo_sites.status.eq("observed_satellite_geometry_only").sum()),
        "source_original_2024_PCNL_valid_islands":int(photo_sites.loc[
            photo_sites.status.eq("observed_satellite_geometry_only"),"island"].nunique()),
        "2024_VIIRS500m_threshold_results":summary,
        "n_orig_valid_sites_with_both_500m_threshold_0_and_1":int(common_valid.sum()),
        "fine500m_fraction_threshold0_equal_threshold1":float(ident.mean()) if len(ident) else None,
        "n_failures_geometry_or_antimeridian":len(failures),
        "geometry_failure_examples":failures[:12],
        "output_CSV_SHA256":source_sha256(path),
        "observed_VIIRS2024_not_historical_extrapolation":True,
        "raster_value_scale_converted_10_to_1":True,
        "coarse_and_fine_different_satellite_products_not_resolution_only":True,
        "500m_pixels_are_not_40m_nor_single_streetlights":True,
        "no_genuine_moth_abundance_pollen_delivery_plant_seed_set_effect_tested":True,
        "source_GSHHG_2p3p7_original_land_minus_original_inland_lakes":True,
        "frozen_2025_photo_site_cohort_not_randomized_global_flora_sample":True,
        "minimum_land_pixel_rule":"64 finer pixels versus 16 coarse pixels",
        "min_real_land_raster_coverage_fraction":.80}
    (out/"viirs500m_vs_pcnl1km_actual_source_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2,default=str)+"\n")
    print("ALAN_REAL_500M_VIIRS_VS_1KM_PCNL_SOURCE",json.dumps(report,sort_keys=True,default=str),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--coarse",type=Path,required=True)
    p.add_argument("--coarse-receipt",type=Path,required=True)
    p.add_argument("--coast",type=Path,required=True)
    p.add_argument("--raster",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.coarse,a.coarse_receipt,a.coast,a.raster,a.out)

if __name__=="__main__":main()
