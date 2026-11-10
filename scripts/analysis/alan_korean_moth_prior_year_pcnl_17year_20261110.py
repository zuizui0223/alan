#!/usr/bin/env python3
"""True original Korean insular moth sampling site-years × PRE-survey PCNL.

Four independent year shards acquire original MD5-verified PCNL (2006-2024)
for the 17 *actually surveyed* years from the independently published
2007-2025 island moth data (DOI 10.5281/zenodo.21438943).
Never use 2024 VIIRS as the 2007-2024 historical field exposure.
No false moth absence, zero sampling effort or effective pollination inferred.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform

from scripts.acquisition.alan_korean_island_moth_real_source_preflight_20261110 import (
    get_original, SOURCES)
from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    manifest_read,fetch_one)
SURVEY_YEARS=(2007,2008,2009,2010,2011,2012,2014,2015,2016,
              2017,2018,2020,2021,2022,2023,2024,2025)
SHARDS=4
SCHEMA="ALAN_true_Korean_30island_moth_yearminus1_PCNL_original_17years_v1"
MANIFEST=Path("config/alan_pcnl_1992_2024_zenodo_v3_md5.json")
SOURCE_ROWS=377

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_KOREAN_PCNL_TRUE: "+str(msg))
def checksum(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for c in iter(lambda:f.read(4<<20),b""):h.update(c)
    return h.hexdigest()

def years_for_shard(shard):
    require(shard in range(SHARDS),"Unknown actual source year shard")
    return [y for i,y in enumerate(SURVEY_YEARS) if i%SHARDS==shard]

def genuine_site_years():
    m=pd.read_excel(io.BytesIO(get_original(*SOURCES["occurrences"])),
                    usecols=["island","year","sites"])
    loc=pd.read_csv(io.BytesIO(get_original(*SOURCES["sites"])),low_memory=False)
    require(len(m)==SOURCE_ROWS and len(loc)==332 and
            set(m.year.astype(int))==set(SURVEY_YEARS) and
            not m.duplicated(["island","year","sites"]).any(),
            "Original 2007-2025 Korean moth site-year metadata changed")
    require(not loc.duplicated(["island","site_id"]).any() and
            {"Latitude","Longitude","Elevation"}.issubset(loc.columns),
            "Authentic site geolocation missing or duplicated")
    x=m.merge(loc,left_on=["island","sites"],right_on=["island","site_id"],
              how="left",validate="many_to_one",indicator=True)
    require(x._merge.eq("both").all() and x.Latitude.between(30,40).all() and
            x.Longitude.between(120,135).all(),"Original Korean moth site coordinates changed")
    require(x[["island","sites","year"]].notna().all().all(),
            "Source site-year identity unknown")
    x["year"]=x.year.astype(int)
    x["raster_year"]=x.year-1
    return x[["island","year","sites","Latitude","Longitude",
              "Elevation","raster_year"]].copy()

def sample_original_raster(src,site):
    require(src.crs is not None and src.count==1 and src.crs.is_geographic,
            "Unsupported original PCNL raster; fail closed")
    xx=site.Longitude.to_numpy(float)
    yy=site.Latitude.to_numpy(float)
    if src.crs.to_epsg()!=4326:
        xx,yy=transform("EPSG:4326",src.crs,xx.tolist(),yy.tolist())
    vals=[]
    for v in src.sample(list(zip(xx,yy)),indexes=1,masked=True):
        if np.ma.is_masked(v[0]):
            vals.append(np.nan)
        else:
            n=float(v[0]);vals.append(n if np.isfinite(n) and n>=0 else np.nan)
    return np.asarray(vals)

def run_shard(shard,out):
    out=Path(out)
    require(not out.exists(),"Do not replace signed original PCNL source shard")
    manifest=manifest_read(MANIFEST)["files"]
    original=genuine_site_years()
    years=years_for_shard(shard)
    selected=original.loc[original.year.isin(years)].copy()
    require(len(selected)>0 and set(selected.year)==set(years),"Shard lost original moth sites")
    out.mkdir(parents=True)
    stats=[]
    parts=[]
    for y in years:
        raster_year=y-1
        raster_name=f"PCNL{raster_year}.tif"
        b=fetch_one(raster_name,manifest[raster_name],out/raster_name)
        g=selected.loc[selected.year.eq(y)].copy()
        with rasterio.open(out/raster_name) as raster:
            v=sample_original_raster(raster,g)
            grid={"crs":str(raster.crs),"res":list(raster.res),
                  "nodata":raster.nodata,"scale":list(raster.scales),
                  "offset":list(raster.offsets)}
        g["prior_year_PCNL_native"]=v
        g["prior_year_log1p_PCNL"]=np.log1p(v)
        parts.append(g)
        stats.append({"original_moth_survey_year":y,
                      "pre_survey_raster_year":raster_year,
                      "source_moth_site_years":len(g),
                      "n_true_PCNl_source_cells":int(np.isfinite(v).sum()),
                      "n_missing_PCNl_source_cells":int(np.isnan(v).sum()),
                      "raster_source":b,"original_raster_geometry":grid})
        print("ALAN_KOREAN_MOTH_GENUINE_PRIOR_PCNL_YEAR",
              json.dumps(stats[-1],sort_keys=True),flush=True)
        (out/raster_name).unlink() # archive only source proof, not rasters
    data=pd.concat(parts,ignore_index=True).sort_values(
        ["island","year","sites"]).reset_index(drop=True)
    require(not data.duplicated(["island","year","sites"]).any(),"Site×year duplicated")
    csv=out/f"korean_moth_true_prior_year_pcnl_shard{shard}.csv.gz"
    data.to_csv(csv,index=False,compression="gzip")
    report={
        "schema":SCHEMA,"shard":shard,"source_moth_DOI":"10.5281/zenodo.21438943",
        "source_PCNL_DOI":"10.5281/zenodo.17013414",
        "actual_source_survey_years":years,
        "real_original_moth_site_years":len(data),
        "original_yearminus1_raster_MD5_audit":stats,
        "real_source_CSV_SHA256":checksum(csv),
        "prior_light_not_future_2024_exposure_for_older_moths":True,
        "no_true_insect_population_density_or_pollen_deposition_estimated":True,
    }
    (out/f"korean_moth_true_prior_year_pcnl_shard{shard}.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

def collect(root,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite actual collected 17-year source")
    root=Path(root)
    original=genuine_site_years()
    sheets=[];seen=[]
    for shard in range(SHARDS):
        data=root/f"korean_moth_true_prior_year_pcnl_shard{shard}.csv.gz"
        rec=root/f"korean_moth_true_prior_year_pcnl_shard{shard}.json"
        report=json.loads(rec.read_text())
        require(report["schema"]==SCHEMA and report["shard"]==shard and
                report["actual_source_survey_years"]==years_for_shard(shard) and
                checksum(data)==report["real_source_CSV_SHA256"] and
                len(report["original_yearminus1_raster_MD5_audit"])==len(years_for_shard(shard)),
                "An actual independent prior-year PCNL source shard was not verified")
        sheet=pd.read_csv(data)
        require(not sheet.duplicated(["island","year","sites"]).any(),
                "Duplicated original survey year/position in joined raster")
        sheets.append(sheet);seen.extend(report["actual_source_survey_years"])
    require(sorted(seen)==sorted(SURVEY_YEARS),"Original 17 survey years are incomplete")
    data=pd.concat(sheets,ignore_index=True)
    require(len(data)==SOURCE_ROWS and
            not data.duplicated(["island","year","sites"]).any() and
            (data.raster_year==data.year-1).all(),"Original site-year full join missing or leaked future year")
    key=["island","year","sites"]
    x=original.merge(data[key+["prior_year_PCNL_native","prior_year_log1p_PCNL"]],
                     on=key,how="left",validate="one_to_one",indicator=True)
    require(x._merge.eq("both").all(),"A true original Korean moth sample was silently omitted")
    out.mkdir(parents=True)
    path=out/"korean_30island_377_actual_moth_siteyear_prior_PCNL_yearminus1.csv.gz"
    x.to_csv(path,index=False,compression="gzip")
    report={"schema":SCHEMA,"all_17_survey_years":list(SURVEY_YEARS),
            "original_survey_site_years":len(x),"original_distinct_islands":int(x.island.nunique()),
            "real_source_PCNL_valid_site_years":int(x.prior_year_log1p_PCNL.notna().sum()),
            "PCNL_missing_site_years":int(x.prior_year_log1p_PCNL.isna().sum()),
            "PCNL_sites_within_survey_year_at_least_2_real_light_levels":int(
                sum(g.prior_year_log1p_PCNL.nunique()>=2 and len(g)>=2
                    for _,g in x.groupby(["island","year"]))),
            "joined_full_source_CSV_SHA256":checksum(path),
            "cannot_infer_moth_persistence_from_unobserved_years":True,
            "not_sameyear_future_exposure":True,
            "not_measured_moth_pollination_effect":True}
    (out/"korean_30island_prior_PCNL_complete_source_receipt.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_30ISLAND_TRUE_MOTH_PCNL17YEAR_COMPLETE",json.dumps(report,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="cmd",required=True)
    w=sub.add_parser("worker")
    w.add_argument("--shard",type=int,required=True)
    w.add_argument("--out",type=Path,required=True)
    c=sub.add_parser("collect")
    c.add_argument("--root",type=Path,required=True)
    c.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.cmd=="worker":run_shard(a.shard,a.out)
    else:collect(a.root,a.out)

if __name__=="__main__":
    main()
