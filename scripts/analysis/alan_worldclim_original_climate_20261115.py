#!/usr/bin/env python3
"""Source-verified WorldClim 2.1 10m BIO climate normal at actual ALAN sites.

BIO1 mean annual temperature, BIO12 annual precipitation and BIO15 rain
seasonality are 1970–2000 spatial CLIMATOLOGICAL normals, not contemporaneous
weather at photographic observation year. Original raster zip SHA is checked
before this workflow; pixel NO_DATA is never assigned as zero.

Non-geographic aggregate output only; use the return arrays in memory in
subsequent paired population-adjusted observational analyses.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY,load_sites,require
)

BIO=("bio1","bio12","bio15")
CLIMATE_SHA="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"

def original_raster_values(tiff,lon,lat):
    with rasterio.open(tiff) as ds:
        require(ds.crs is not None and ds.crs.to_epsg()==4326,
                "WorldClim source must be WGS84")
        require(ds.width==2160 and ds.height==1080 and ds.count==1,
                "WorldClim original 10 arcminute raster shape changed")
        require(abs(abs(ds.res[0])-1/6)<1e-7 and
                abs(abs(ds.res[1])-1/6)<1e-7,
                "WorldClim source not 10 arcminute scale")
        x=np.asarray(lon,float)
        y=np.asarray(lat,float)
        valid_coord=np.isfinite(x)&np.isfinite(y)
        rr=np.zeros(len(x),dtype=np.int64)
        cc=np.zeros(len(x),dtype=np.int64)
        xx=np.where(valid_coord,x,0)
        yy=np.where(valid_coord,y,0)
        rr[:],cc[:]=rasterio.transform.rowcol(ds.transform,xx,yy,op=np.floor)
        valid=valid_coord&(rr>=0)&(rr<ds.height)&(cc>=0)&(cc<ds.width)
        arr=ds.read(1)
        vals=np.full(len(x),np.nan,dtype=np.float64)
        rows=rr[valid];cols=cc[valid]
        vv=arr[rows,cols].astype(float)
        keep=np.isfinite(vv)
        if ds.nodata is not None and np.isfinite(ds.nodata):
            keep&=vv!=float(ds.nodata)
        iv=np.flatnonzero(valid)
        vals[iv[keep]]=vv[keep]
        meta={"crs":str(ds.crs),"width":ds.width,"height":ds.height,
              "resolution":[float(q) for q in ds.res],"nodata":str(ds.nodata),
              "valid":int(np.isfinite(vals).sum()),
              "missing":int(np.isnan(vals).sum())}
        return vals,meta

def sample_sites(sites,raster_root):
    require(not sites.duplicated(KEY).any(),"Original photographic source site duplicates")
    z=sites.copy().reset_index(drop=True)
    m={}
    for k,i in (("bio1",1),("bio12",12),("bio15",15)):
        p=Path(raster_root)/f"wc2.1_10m_bio_{i}.tif"
        require(p.exists(),"Missing independently sourced WorldClim BIO raster")
        z[k],m[k]=original_raster_values(p,z.lon,z.lat)
    z["climate_complete"]=np.isfinite(z[list(BIO)]).all(axis=1)
    require((z.loc[z.climate_complete,"bio12"]>=0).all(),
            "Impossible negative annual precipitation in source")
    require((z.loc[z.climate_complete,"bio15"]>=0).all(),
            "Impossible negative precipitation seasonality in source")
    for key in BIO:
        # Static climate not contemporaneous: value same for repeated grid cell
        same=z.groupby(["island","latbin","lonbin"],observed=True)[key].nunique(dropna=True)
        require(same.max()<=1,"Static climate shifted across years at fixed original site")
    return z,m

def audit(z,m):
    valid=z.climate_complete
    count=len(z)
    island=z.island.nunique()
    year=z.groupby("year",observed=True).agg(
        original_source_site_years=("island","size"),
        worldclim_allthree_valid=("climate_complete","sum"))
    site=z.loc[valid].copy()
    year_island_g=site.groupby(["island","year"],observed=True)
    iy=year_island_g.agg(
       n=("island","size"),
       bio1_range=("bio1",lambda x:float(x.max()-x.min())),
       bio12_range=("bio12",lambda x:float(x.max()-x.min())),
       bio15_range=("bio15",lambda x:float(x.max()-x.min())))
    climate_cov={
        "schema":"ALAN_WorldClim2p1_original10arcminute_spatial_normal_coverage_v1",
        "status":"TRUE_ORIGINAL_RASTER_SOURCE_SAMPLE_NOT_A_TEMPORAL_CLIMATE_MODEL",
        "source":"WorldClim v2.1 10m BIO1/BIO12/BIO15, static baseline normals 1970–2000",
        "original_alan_site_years":count,
        "original_alan_islands":island,
        "climate_complete_site_years":int(valid.sum()),
        "climate_complete_islands":int(site.island.nunique()),
        "climate_missing_site_years":int((~valid).sum()),
        "annual_photo_year_source_coverage":[
          {"year":int(y),"total":int(row.original_source_site_years),
           "complete":int(row.worldclim_allthree_valid)}
          for y,row in year.iterrows()],
        "site_climate_within_island_year_support":{
          "groups_with_two_sites":int(iy.n.ge(2).sum()),
          "groups_with_temperature_range_at_least_0p1C":int(
               (iy.n.ge(2)&iy.bio1_range.ge(.1)).sum()),
          "groups_with_precipitation_range_at_least_10mm":int(
               (iy.n.ge(2)&iy.bio12_range.ge(10)).sum()),
          "groups_with_both_temp_and_precip_variation":int(
               (iy.n.ge(2)&iy.bio1_range.ge(.1)&iy.bio12_range.ge(10)).sum()),
          "islands_with_two_site_climate_gradient":int(
               iy.loc[iy.n.ge(2)].index.get_level_values("island").nunique())
        },
        "bioclimatic_variable_grid_metadata":m,
        "climate_NODATA_retained":True,
        "source_is_10arcmin_static_1970_2000_normal_not_2019_2025_actual_weather":True,
        "no_raw_species_or_location_export":True,
        "not_locally_measured_flower_temperature":True
    }
    return climate_cov

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites",type=Path,required=True)
    p.add_argument("--raster_dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    require(not a.out.exists(),"Do not overwrite source receipt")
    site=load_sites(a.sites)
    require((len(site),site.island.nunique())==(77647,562),
            "Frozen ALAN photo-site source changed")
    z,meta=sample_sites(site,a.raster_dir)
    obj=audit(z,meta)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,indent=2,allow_nan=False,ensure_ascii=False)+"\n")
    print("ALAN_WORLDCLIM_10M_BIO_ORIGINAL_SOURCE_AUDIT",
          json.dumps(obj,ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":main()
