#!/usr/bin/env python3
"""Extract same-calendar-year WorldPop R2024B unconstrained population at ALL
frozen ALAN source photographic site coordinates, without MOIS photo overlap.

Input ALAN source is SHA-locked 85,318 site/year panel, reduced by its frozen
original eligible function to 77,647, 562 islands. Original WorldPop 30 arcsec
raster bytes for 2019,2021,2023,2024,2025 are individually SHA-pinned from
previously successful independent MOIS source acquisition.

WorldPop values are people PER 1km source raster pixel, converted to people
per geodesic square kilometre at that exact pixel. Missing != zero. No
biological species, taxon/trait, photographs or raw location data exported.
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
from rasterio.windows import Window, bounds
from pyproj import Geod
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, load_sites, require
)

RASTERS = {
    2019: ("16195d62c18a0a524d5742c641f4c3f4abf363750e35d3bbc4eac23ea26f2ae2",870798761,8257),
    2021: ("4a6bfc000230c1f81ab898304af4803f967e147c53b2373df37a65a48448eb81",873971195,12542),
    2023: ("4a182c22d8523691d9dcbc5c03bca164a8a05bff1445ec738e8c5357dba3a74b",869736358,17025),
    2024: ("e0ff729885f5449eb5bdd7347f9b337d99a8ee9bbbaadf2b3e679db8d8f53538",868021365,18890),
    2025: ("4ef9014dbe9f20462b719abb2987eb9efca2b4517bfdcd7f0f8d0a2a8cdd191f",865211479,20933),
}
BASE_URL="https://worldpop-public-data.soton.ac.uk/GIS/Population/Global_2015_2030/R2024B"
EXPECTED_N=77647
EXPECTED_ISLANDS=562

def file_sha256(path:Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4<<20),b""):
            h.update(block)
    return h.hexdigest()

def geodesic_pixel_km2(ds,row:int,geod:Geod) -> float:
    west,south,east,north=bounds(Window(0,int(row),1,1),ds.transform)
    area,_=geod.polygon_area_perimeter(
        [west,east,east,west],[south,south,north,north])
    km2=abs(area)/1.e6
    require(np.isfinite(km2) and km2>0,
            "Invalid WorldPop geodesic source pixel area")
    return float(km2)

def sample_source(ds, original:pd.DataFrame):
    require(len(original)>0,"No ALAN original source photo sites")
    require(ds.count==1 and ds.crs.to_epsg()==4326 and
            ds.width==43200 and ds.height==17280,
            "Wrong official WorldPop raster shape/CRS")
    require(abs(ds.res[0]-1/120)<1.e-7 and abs(ds.res[1]-1/120)<1.e-7,
            "WorldPop source not 30 arcseconds")
    xs=original.lon.to_numpy(float)
    ys=original.lat.to_numpy(float)
    valid_coords=np.isfinite(xs)&np.isfinite(ys)&(
        xs>=ds.bounds.left)&(xs<ds.bounds.right)&(
        ys>ds.bounds.bottom)&(ys<=ds.bounds.top)
    rowcol=np.array([ds.index(float(x),float(y))
                     if ok else (-1,-1)
                     for x,y,ok in zip(xs,ys,valid_coords)],dtype=int)
    within=valid_coords & (rowcol[:,0]>=0)&(rowcol[:,0]<ds.height)&(
        rowcol[:,1]>=0)&(rowcol[:,1]<ds.width)
    n=len(original)
    count=np.full(n,np.nan,float)
    density=np.full(n,np.nan,float)
    status=np.full(n,"OUTSIDE_RASTER",object)
    validids=np.flatnonzero(within)
    if len(validids):
        bh,bw=ds.block_shapes[0]
        # Sort actual raster block row/col to minimize repeated tile reads.
        order=np.lexsort((rowcol[validids,1]//bw,rowcol[validids,0]//bh))
        validids=validids[order]
        geod=Geod(ellps="WGS84")
        cache={}
        coords=((float(xs[i]),float(ys[i])) for i in validids)
        for i,sampled in zip(validids,ds.sample(coords,indexes=1,masked=True)):
            if np.ma.is_masked(sampled[0]):
                status[i]="NODATA"
                continue
            v=float(sampled[0])
            if not np.isfinite(v) or (
                ds.nodata is not None and v==float(ds.nodata)):
                status[i]="NODATA"
                continue
            if v<0:
                status[i]="NEGATIVE_INVALID"
                continue
            rr=int(rowcol[i,0])
            if rr not in cache:
                cache[rr]=geodesic_pixel_km2(ds,rr,geod)
            count[i]=v
            density[i]=v/cache[rr]
            status[i]="OK"
    require((np.isfinite(density)==(status=="OK")).all(),
            "Population invalid or missing source status mixed with valid data")
    require(np.all(density[status=="OK"]>=0),
            "True population cannot be negative")
    return count,density,status,within

def run(year:int,site_file:Path,raster:Path,outdir:Path):
    require(year in RASTERS,"No authenticated year")
    require(not outdir.exists(),"Do not overwrite verified source outputs")
    expected_sha, expected_bytes, expected_year_rows=RASTERS[year]
    expect_name=f"global_pop_{year}_UC_1km_R2024B_v1.tif"
    require(raster.name==expect_name and raster.stat().st_size==expected_bytes,
            "Wrong original file or byte size")
    original_raster_sha=file_sha256(raster)
    require(original_raster_sha==expected_sha,
            "WorldPop source SHA does not match frozen original public raster")
    all_sites=load_sites(site_file)
    require((len(all_sites),all_sites.island.nunique())==
            (EXPECTED_N,EXPECTED_ISLANDS),
            "ALAN original photo site population drift")
    all_sites=all_sites.copy().reset_index(drop=True)
    all_sites["site_id"]=np.arange(len(all_sites),dtype=np.int32)
    s=all_sites.loc[all_sites.year.eq(year)].copy()
    require(len(s)==expected_year_rows,"Frozen ALAN year source row drift")
    require(not s.duplicated(KEY).any(),"Duplicate site-year keys")
    with rasterio.open(raster) as ds:
        counts,density,status,within=sample_source(ds,s)
        metadata={"crs":str(ds.crs),"width":ds.width,"height":ds.height,
                  "res":[float(x) for x in ds.res],
                  "nodata":float(ds.nodata) if ds.nodata is not None else None,
                  "bounds":[float(x) for x in ds.bounds],
                  "source_original_sha256":original_raster_sha}
    # Stable site_id is a row index in the original checksum-locked table.
    # No source lat/lon/accepted taxon names included in saved artifact.
    out=pd.DataFrame({
        "site_id":s.site_id.to_numpy(np.int32),
        "original_observation_year":np.full(len(s),year,dtype=np.int16),
        "worldpop_year":np.full(len(s),year,dtype=np.int16),
        "population_status":status,
        "population_people_per_pixel":counts,
        "population_density_per_km2":density,
    })
    require(out.site_id.is_unique and len(out)==len(s),
            "Original ALAN site rows must be preserved one-to-one")
    csv=outdir/f"worldpop_R2024B_UC_{year}_original_ALAN_photo_sites.csv.gz"
    outdir.mkdir(parents=True,exist_ok=True)
    out.to_csv(csv,index=False,compression={"method":"gzip","mtime":0},
               float_format="%.11g",lineterminator="\n")
    receipt={
        "schema":"ALAN_original_full_photo_island_WorldPop_G2_R2024B_UC_v1",
        "status":"ORIGINAL_SOURCE_ANNUAL_GEOSAMPLED_NO_FILL",
        "photo_site_source_SHA256":
           "440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c",
        "population_year":year,
        "original_global_raster":expect_name,
        "raster_release":"WorldPop Global2 R2024B unconstrained v1",
        "original_raster_bytes":expected_bytes,
        "raster_metadata":metadata,
        "same_year_only":True,
        "original_alansites":len(s),
        "original_islands_with_photo_sites":int(s.island.nunique()),
        "sampled_valid_population":int((status=="OK").sum()),
        "sampled_nodata":int((status=="NODATA").sum()),
        "sampled_outside_raster":int((status=="OUTSIDE_RASTER").sum()),
        "sampled_negative_invalid":int((status=="NEGATIVE_INVALID").sum()),
        "zero_people_pixels":int(((status=="OK")&(counts==0)).sum()),
        "new_population_source_values_sha256":file_sha256(csv),
        "no_source_geography_species_rows_exported":True,
        "no_missing_population_imputation":True,
        "source_photo_ids_not_biological_absences":True
    }
    (outdir/"source_receipt.json").write_text(
        json.dumps(receipt,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_FULL_WORLDPOP_AUTHENTIC_YEAR",
          json.dumps({k:receipt[k] for k in (
               "population_year","original_alansites",
               "sampled_valid_population","sampled_nodata",
               "sampled_outside_raster","zero_people_pixels",
               "original_raster_bytes","new_population_source_values_sha256"
          )},allow_nan=False),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites",type=Path,required=True)
    p.add_argument("--raster",type=Path,required=True)
    p.add_argument("--year",type=int,required=True,choices=sorted(RASTERS))
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.year,a.sites,a.raster,a.out)
