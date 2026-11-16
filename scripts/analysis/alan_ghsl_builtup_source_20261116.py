#!/usr/bin/env python3
"""Official JRC GHSL GHS-BUILT-S R2023A, 30arcsec WGS84 at ALAN photo sites.

Built-up surface is square metres per ORIGINAL raster cell. Convert to
built fraction by native source pixel geodesic area, not by WorldPop person
density, 0.05° photographic bin area, or a hardcoded km² value.

Epoch 2015 predates every original 2019–2025 photographic year.
Epoch 2020 is *sensitivity only* (post-2019, not a 2019 baseline).
A source release is authenticated by public S3 key, source byte count
and complete TIFF SHA256. Snapshot into digest-pinned artifacts.

No flower, private botanical species coordinate or raw species data export.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from rasterio.windows import Window, bounds
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    load_sites, require
)

SOURCE={
  2015:("GHS_BUILT_S_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif",170932198),
  2020:("GHS_BUILT_S_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif",173085842),
}
S3_PREFIX="s3://jrc-ghsl/ghs-built-s/r2023a/4326/30asec/"
ORIGINAL_SITES_SHA="440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c"

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as stream:
        for part in iter(lambda:stream.read(4<<20),b""):
            h.update(part)
    return h.hexdigest()

def geodesic_cell_area(ds,row,geod):
    x0,y0,x1,y1=bounds(Window(0,int(row),1,1),ds.transform)
    area,_=geod.polygon_area_perimeter(
        [x0,x1,x1,x0],[y0,y0,y1,y1])
    square_m=float(abs(area))
    require(square_m>0 and math.isfinite(square_m),
            "Original GHSL native geodesic cell area invalid")
    return square_m

def sample_ghsl(ds,lon,lat):
    require(ds.crs is not None and ds.crs.to_epsg()==4326 and
            ds.count==1 and ds.width==43200 and ds.height==21600,
            "Official GHS-BUILT-S must be global WGS84 30arcsec")
    require(abs(ds.res[0]-1/120)<1e-8 and abs(ds.res[1]-1/120)<1e-8,
            "Not original GHSL 30arcsecond resolution")
    n=len(lon)
    lon=np.asarray(lon,float);lat=np.asarray(lat,float)
    valid=(np.isfinite(lon)&np.isfinite(lat)&
           (lon>=ds.bounds.left)&(lon<ds.bounds.right)&
           (lat>ds.bounds.bottom)&(lat<=ds.bounds.top))
    rc=np.asarray([ds.index(float(lon[i]),float(lat[i]))
                   if valid[i] else (-1,-1)
                   for i in range(n)],dtype=np.int32)
    valid=valid&(rc[:,0]>=0)&(rc[:,0]<ds.height)&(
           rc[:,1]>=0)&(rc[:,1]<ds.width)
    out=np.full(n,np.nan,float)
    status=np.full(n,"OUTSIDE",object)
    raw=np.full(n,np.nan,float)
    ids=np.flatnonzero(valid)
    bh,bw=ds.block_shapes[0]
    # Batch original COG TIFF raster blocks for efficient random access.
    order=np.lexsort((rc[ids,1]//bw,rc[ids,0]//bh))
    ids=ids[order]
    geod=Geod(ellps="WGS84")
    areas={}
    for i,val in zip(ids,ds.sample(
            ((float(lon[i]),float(lat[i])) for i in ids),
             indexes=1,masked=True)):
        if np.ma.is_masked(val[0]):
            status[i]="NODATA";continue
        v=float(val[0])
        if not np.isfinite(v) or (
              ds.nodata is not None and v==float(ds.nodata)):
            status[i]="NODATA";continue
        if v<0:
            status[i]="INVALID_NEGATIVE";continue
        rr=int(rc[i,0])
        if rr not in areas:
            areas[rr]=geodesic_cell_area(ds,rr,geod)
        raw[i]=v
        out[i]=v/areas[rr]
        status[i]="OK"
    require((np.isfinite(out)==(status=="OK")).all(),
            "GHSL nodata misinterpreted as built-up zero")
    require(np.nanmin(out)>=0 and np.nanmax(out)<1.1,
            "Impossible GHSL built-up area fraction; check raster cell units")
    return out,raw,status,{
       "valid":int((status=="OK").sum()),
       "zero_built_area":int(((status=="OK")&(raw==0)).sum()),
       "nodata":int((status=="NODATA").sum()),
       "outside":int((status=="OUTSIDE").sum()),
       "invalid":int((status=="INVALID_NEGATIVE").sum()),
       "fraction_over_1":int(np.sum(out>1)),
       "max_native_fraction":float(np.nanmax(out)),
       "raster_block_shape":[int(x) for x in ds.block_shapes[0]]
    }

def run(site_file,root,out):
    require(not out.exists(),"Never overwrite signed source output")
    original=load_sites(site_file).reset_index(drop=True)
    require(len(original)==77647 and original.island.nunique()==562,
            "Original 562-island light sample drift")
    result=pd.DataFrame({"site_id":np.arange(len(original),dtype=np.int32)})
    source={}
    for epoch,(name,size) in SOURCE.items():
        raster=Path(root)/name
        require(raster.exists() and raster.stat().st_size==size,
                f"Wrong GHSL {epoch} official raster path/byte size")
        digest=sha256(raster)
        with rasterio.open(raster) as ds:
            val,raw,status,meta=sample_ghsl(
                ds,original.lon.to_numpy(float),original.lat.to_numpy(float))
            other={
              "crs":str(ds.crs),"width":ds.width,"height":ds.height,
              "resolution":[float(x) for x in ds.res],
              "nodata":None if ds.nodata is None else str(ds.nodata)
            }
        result[f"ghsl_built_{epoch}_fraction"]=val
        result[f"ghsl_built_{epoch}_status"]=status
        source[str(epoch)]={
           "official_s3_key":S3_PREFIX+str(epoch)+"/"+name,
           "original_source_byte_count":size,
           "sha256_original_unmodified_tif":digest,
           "native_units":"square metres built-up surface per WGS84 30arcsec pixel",
           "fraction_calculation":"built area m² / exact native 30arcsec WGS84 geodesic pixel m²",
           "status":meta,"raster":other,
           "epoch_2020_is_future_for_2019_photo_year":epoch==2020
        }
    require(result.site_id.is_unique and len(result)==77647,"ALAN site row loss")
    require((result.ghsl_built_2015_fraction>=0).fillna(True).all() and
            (result.ghsl_built_2020_fraction>=0).fillna(True).all(),
            "GHSL negative building fractions")
    # Sanity: source built-up after 2015 should not systematically vanish,
    # but do not enforce monotonicity because epochs are estimated and
    # land reclamation, demolitions and mapping uncertainty exist.
    both=result.dropna(subset=[
        "ghsl_built_2015_fraction","ghsl_built_2020_fraction"])
    source["temporal_epoch_audit"]={
      "site_years_with_both_epochs":len(both),
      "2020_lower_built_fraction_than_2015":int(
          both.ghsl_built_2020_fraction.lt(both.ghsl_built_2015_fraction).sum()),
      "2015_and_2020_source_values_do_not_imply_actual_light_year":True
    }
    out.mkdir(parents=True)
    artifact=out/"alan_77647_original_ghsl_build_2015_2020_source_sitevalues.csv.gz"
    result.to_csv(artifact,index=False,compression={"method":"gzip","mtime":0},
                  float_format="%.10g",lineterminator="\n")
    receipt={
      "schema":"ALAN_77k_photo_original_GHSL_BUILT_S_R2023A_30asec_epochs2015_2020_v1",
      "status":"INDEPENDENT_JRC_GHSL_ORIGINAL_GEOGRAPHIC_SOURCE_AUDIT",
      "original_ALAN_photo_site_years":len(original),
      "original_ALAN_islands":int(original.island.nunique()),
      "source_2015_independent_prior_urban_builtup":True,
      "source_2020_sensitivity_only":True,
      "source_release_S3_public_open_data":True,
      "source":source,
      "output_site_values_sha256":sha256(artifact),
      "no_botanical_species_x_location_rows_exported":True,
      "no_builtup_source_nodata_imputed_zero":True,
      "not_landcover_plant_habitat_or_lamp_illuminance":True
    }
    (out/"source_receipt.json").write_text(
        json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print("ALAN_GHSL_BUILTUP_AUTHENTIC_SOURCE",json.dumps({
      "ALAN_sites":len(original),
      "2015":source["2015"],"2020":source["2020"],
      "combined_epochs":source["temporal_epoch_audit"],
      "output_SHA":receipt["output_site_values_sha256"]
    },allow_nan=False),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("sites","ghsl_root","out"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.sites,a.ghsl_root,a.out)
