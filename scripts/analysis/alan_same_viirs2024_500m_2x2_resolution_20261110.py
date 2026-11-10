#!/usr/bin/env python3
"""Isolate pixel-grid resolution within *one actual* 2024 VIIRS VNL V2.2 COG.

Use same land-only geodesic 5 km 2025 photographed sites, compare original
15-arc-second (~500m) binary dark patch with geographically aligned 2x2
(30-arc-second ~1km) quartets made from those exact same authentic pixels.
Only *all-land, all-inside-radius, all-valid* quartets enter either side of
a pair: no coastline ocean mixing or imputed nodata. This is satellite
geometry only, not measured moths or a causal resolution effect on biology.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from rasterio.windows import Window, from_bounds
from scipy import ndimage
from shapely import contains_xy, union_all

from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land

RASTER_MD5="0b976892d40c2d0d99177c727289eaf2"
FINE_SHA="68ed09378a8bfbcb36572ff0b0578334e4f1795efe22209f437960dbb01d5896"
SCHEMA="ALAN_2024_true_same_VIIRS_product_500m_vs_2x2_1km_dark_mosaic_v1"
GEOD=Geod(ellps="WGS84")
T=(0.,1.,3.)
STRUCTURE=np.array([[0,1,0],[1,1,1],[0,1,0]],dtype=np.uint8)

def require(condition,message):
    if not bool(condition):raise ValueError("ALAN_SAME_VIIRS_RES: "+str(message))
def digest(path,name="sha256"):
    h=hashlib.new(name)
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def patch_fraction(mask,weights):
    """Measured largest 4-connected component divided by all measured dark."""
    labels,n=ndimage.label(mask,STRUCTURE)
    if not n: return None
    sizes=np.bincount(labels[mask],weights=weights[mask],minlength=n+1)
    return float(sizes[1:].max()/sum(weights[mask]))

def quantile(x):
    d=pd.to_numeric(x,errors="coerce").dropna()
    return {str(p):float(d.quantile(p)) for p in (.05,.25,.50,.75,.95)} if len(d) else {}

def subset_and_coarsen(src,polygon,lon,lat,radius_km=5.):
    require(src.crs.to_epsg()==4326 and src.count==1 and
            np.isclose(src.res[0],1/240,atol=1e-10) and
            np.isclose(src.res[1],1/240,atol=1e-10) and
            np.isclose(src.scales[0],10.) and np.isclose(src.offsets[0],0.),
            "Wrong authentic original VIIRS 15-arcsecond COG and scale")
    require(np.isfinite([lon,lat]).all() and abs(lon)<179.5 and abs(lat)<80,
            "Site requires explicit antimeridian/polar raster treatment")
    dy=radius_km/110.57
    dx=radius_km/(111.32*np.cos(np.deg2rad(lat)))
    w=from_bounds(lon-dx,lat-dy,lon+dx,lat+dy,src.transform)
    # Align the source raster window to the GLOBAL even-row/even-column
    # origin, so all geographically aligned quartets are reproducible.
    c0=int(np.floor(w.col_off))-2
    r0=int(np.floor(w.row_off))-2
    c1=int(np.ceil(w.col_off+w.width))+2
    r1=int(np.ceil(w.row_off+w.height))+2
    left=2*(c0//2)
    top=2*(r0//2)
    right=2*((c1+1)//2)
    bottom=2*((r1+1)//2)
    win=Window(left,top,right-left,bottom-top)
    require(win.width>0 and win.height>0 and int(win.width)%2==0 and
            int(win.height)%2==0 and win.width*win.height<100000,
            "Not a safe 2x2 aligned geodesic source window")
    raw=src.read(1,window=win,boundless=True,masked=True).astype(float).filled(np.nan)
    ny,nx=raw.shape
    trans=src.window_transform(win)
    rr,cc=np.indices(raw.shape)
    xx=trans.c+(cc+.5)*trans.a+(rr+.5)*trans.b
    yy=trans.f+(rr+.5)*trans.d+(cc+.5)*trans.e
    _,_,dist=GEOD.inv(np.full(raw.shape,lon),np.full(raw.shape,lat),xx,yy)
    inside=dist<=radius_km*1000.
    land=inside & contains_xy(polygon,xx,yy)
    valid=land & np.isfinite(raw) & (raw>=0)
    small=(ny//2,nx//2)
    block_valid=valid.reshape(small[0],2,small[1],2).all(axis=(1,3))
    # The full land-only quartet exactly tiles a source raster cell group;
    # partial land/coastline/radius quartets are excluded on BOTH grids.
    eligible=np.repeat(np.repeat(block_valid,2,axis=0),2,axis=1)
    nblock=int(block_valid.sum())
    if nblock<16:
        return [{"status":"insufficient_16_full_source_land_quartets",
                 "threshold_nW_cm2_sr":t,"n_complete_land_quartets":nblock}
                for t in T]
    weights=np.maximum(np.cos(np.deg2rad(yy)),0.)
    fine_weights=weights*eligible
    block_weight=weights.reshape(small[0],2,small[1],2).sum(axis=(1,3))*block_valid
    denom=float(np.sum(block_weight))
    require(denom>0 and np.isfinite(denom),"Invalid true original geodetic pixel weights")
    raw_coarse=(np.where(eligible,raw,0.)*weights).reshape(
        small[0],2,small[1],2).sum(axis=(1,3))/np.maximum(block_weight,1e-12)
    out=[]
    for t in T:
        fine_dark=eligible & (raw<=t)
        dark_count=fine_dark.reshape(small[0],2,small[1],2).sum(axis=(1,3))
        coarse_majority=block_valid & (dark_count>=3)
        coarse_mean=block_valid & (raw_coarse<=t)
        fine_dark_fraction=float(fine_weights[fine_dark].sum()/denom)
        coarse_majority_fraction=float(block_weight[coarse_majority].sum()/denom)
        coarse_mean_fraction=float(block_weight[coarse_mean].sum()/denom)
        out.append({
            "status":"same_VIIRS_paired_original_pixels",
            "threshold_nW_cm2_sr":t,
            "n_complete_land_quartets":nblock,
            "n_fine_500m_original_land_pixels":nblock*4,
            "fine_dark_land_fraction_on_complete_quartets":fine_dark_fraction,
            "coarse_1km_majority_dark_land_fraction":coarse_majority_fraction,
            "coarse_1km_mean_radiance_dark_land_fraction":coarse_mean_fraction,
            "fine_500m_largest_dark_component_fraction_of_dark":
                patch_fraction(fine_dark,fine_weights),
            "coarse_1km_majority_largest_dark_component_fraction_of_dark":
                patch_fraction(coarse_majority,block_weight),
            "coarse_1km_mean_largest_dark_component_fraction_of_dark":
                patch_fraction(coarse_mean,block_weight),
            "fine_500m_n_dark_original_pixels":int(fine_dark.sum()),
            "coarse_1km_n_majority_dark_blocks":int(coarse_majority.sum()),
            "coarse_1km_n_mean_dark_blocks":int(coarse_mean.sum()),
        })
    return out

def run(fine_csv,fine_receipt,coast,raster,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite true satellite same-product result")
    metadata=json.loads(Path(fine_receipt).read_text())
    require(metadata["output_CSV_SHA256"]==FINE_SHA and
            digest(fine_csv)==FINE_SHA and digest(raster,"md5")==RASTER_MD5,
            "2024 authentic source COG or 2025 photo-site subset drift")
    f=pd.read_csv(fine_csv,dtype={"island":str})
    sites=f.loc[f.VIIRS_dark_threshold_nW_cm2_sr.eq(0.) &
                f.status.eq("observed_satellite_geometry_only")].copy()
    require(len(sites)==3513 and sites.island.nunique()==957 and
            not sites.duplicated(["island","latbin","lonbin"]).any(),
            "True 2025 VIIRS source site identity not the original 3513")
    land=load_signed_land(coast).set_index("land_id",drop=False)
    require(sites.island.isin(land.index).all(),"GSHHG island ID drift")
    lakes=gpd.read_file(Path(coast)/"GSHHS_h_L2.shp").to_crs(4326)
    index=lakes.sindex
    records=[];failures=[]
    with rasterio.open(raster) as src:
        for island,g in sites.groupby("island",sort=True):
            poly=land.loc[island,"geometry"]
            isects=index.query(poly,predicate="intersects")
            if len(isects):
                rings=[p for p in lakes.geometry.iloc[isects] if p.intersects(poly)]
                if rings:poly=poly.difference(union_all(rings))
            if poly.is_empty or not poly.is_valid:
                failures.append({"island":island,"reason":"original island geometry invalid after lake mask"})
                continue
            for s in g.itertuples(index=False):
                try:
                    a=subset_and_coarsen(src,poly,float(s.lon),float(s.lat))
                except (ValueError,TypeError) as err:
                    failures.append({"island":island,"latbin":s.latbin,
                                     "lonbin":s.lonbin,"reason":str(err)})
                    continue
                for v in a:
                    records.append({"island":island,"latbin":s.latbin,
                                    "lonbin":s.lonbin,"photo_year":2025,"VIIRS_year":2024,
                                    **v})
    d=pd.DataFrame(records)
    require(len(d)>0 and not d.duplicated(["island","latbin","lonbin",
                                           "threshold_nW_cm2_sr"]).any(),
            "Fake or duplicate 500m/2x2 same-source satellite row")
    out.mkdir(parents=True)
    csv=out/"genuine_same_VIIRS2024_500m_vs_1km_matched_complete_land_quartets.csv.gz"
    d.to_csv(csv,index=False,compression="gzip")
    stats={}
    for t in T:
        x=d.loc[d.threshold_nW_cm2_sr.eq(t)]
        x=x.loc[x.status.eq("same_VIIRS_paired_original_pixels")]
        cond=x.loc[x.fine_500m_largest_dark_component_fraction_of_dark.notna() &
                   x.coarse_1km_majority_largest_dark_component_fraction_of_dark.notna()]
        stats[str(t)]={
            "same_source_support_sites":len(x),
            "same_source_support_islands":int(x.island.nunique()),
            "n_both_grids_have_dark_connected_regions":len(cond),
            "fine_dark_fraction_median":float(x.fine_dark_land_fraction_on_complete_quartets.median())
                if len(x) else None,
            "coarse_majority_dark_fraction_median":float(x.coarse_1km_majority_dark_land_fraction.median())
                if len(x) else None,
            "coarse_mean_radiance_dark_fraction_median":float(x.coarse_1km_mean_radiance_dark_land_fraction.median())
                if len(x) else None,
            "fine_conditional_connectivity_median":
                float(cond.fine_500m_largest_dark_component_fraction_of_dark.median()) if len(cond) else None,
            "coarse_majority_conditional_connectivity_median":
                float(cond.coarse_1km_majority_largest_dark_component_fraction_of_dark.median()) if len(cond) else None,
            "coarse_mean_conditional_connectivity_median":
                float(x.coarse_1km_mean_largest_dark_component_fraction_of_dark.median())
                if x.coarse_1km_mean_largest_dark_component_fraction_of_dark.notna().any() else None,
            "paired_fine_minus_coarse_majority_dark_fraction_quantiles":
                quantile(x.fine_dark_land_fraction_on_complete_quartets-
                         x.coarse_1km_majority_dark_land_fraction),
            "paired_fine_minus_coarse_majority_conditional_connectivity_quantiles":
                quantile(cond.fine_500m_largest_dark_component_fraction_of_dark-
                         cond.coarse_1km_majority_largest_dark_component_fraction_of_dark),
            "n_coarse_majority_0p95_plus_but_fine_below_0p90":int(
                (cond.coarse_1km_majority_largest_dark_component_fraction_of_dark.ge(.95) &
                 cond.fine_500m_largest_dark_component_fraction_of_dark.lt(.90)).sum())
        }
    receipt={
        "schema":SCHEMA,"source_original_500m_COG_MD5":RASTER_MD5,
        "source_original_2025_photographed_site_exposure_SHA256":FINE_SHA,
        "source_original_500m_valid_photo_sites":len(sites),
        "source_original_500m_valid_islands":int(sites.island.nunique()),
        "source_actual_fine_and_2x2_comparison":stats,
        "source_CSV_SHA256":digest(csv),
        "n_geometry_or_window_failures":len(failures),
        "failure_examples":failures[:12],
        "land_quartet_rule":"all 4 original 500m pixel centers on same island land AND inside 5km AND valid; exclude coastal/missing/radius partial quartets",
        "coarse_dark_rules":">=3 of 4 true fine dark pixels (strict majority) and alternative mean raw radiance <= same threshold",
        "sample_selection":"only genuine 2025 photo sites with valid 500m original source; same cells for both resolutions",
        "same_sensor_year_and_product_resolution_comparison":True,
        "different_pixel_aggregation_binary_rule_changes_estimand":True,
        "no_new_true_moth_insect_records_or_effective_pollination":True
    }
    (out/"same_VIIRS_true_original_quartet_resolution_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("ALAN_TRUE_VIIRS_SAME_PRODUCT_500M_VS_1KM",json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser()
    for key in ("fine","fine_receipt","coast","raster","out"):
        p.add_argument("--"+key.replace("_","-"),required=True,type=Path)
    a=p.parse_args()
    run(a.fine,a.fine_receipt,a.coast,a.raster,a.out)

if __name__=="__main__":main()
