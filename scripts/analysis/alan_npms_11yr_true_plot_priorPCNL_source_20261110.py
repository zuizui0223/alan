#!/usr/bin/env python3
"""NPMS 2015-2025 genuine repeatedly visited UK plant plots × annual prior PCNL.

This produces a time-varying SATELLITE EXPOSURE frame *independent of plant
cover/outcomes*: true NPMS sample year 2015–2025 at stable source plot ID,
original Li&Cao annual PCNL 2014–2024 at source plot's fixed coordinate,
strictly one year BEFORE the original plant visit year. Publishes only
coarse site-ID/year exposure and audited source support. No species
nonrecord imputation or GloPL PL transferred to individual plot.

The original NPMS EIDC CSV member CRC/source size pins avoid nondeterministic
outer ZIP SHA. A plot with >500m coordinate drift across original visits
is UNSUPPORTED, not an invented independently revisitable plot.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform
from pyproj import Geod

from scripts.analysis.alan_npms_glopl_actual_species_repeat_support_20261110 import NPMS_MEMBERS
from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one,manifest_read)

DOI="10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8"
PCNL_DOI="10.5281/zenodo.17013414"
SCHEMA="ALAN_true_2015_2025_NPMS_repeated_plot_previous_PCNL_11years_v1"
YEAR_RANGE=list(range(2015,2026))
SHARDS=3
GEOD=Geod(ellps="WGS84")
MAX_SITE_DRIFT_M=500.

def require(x,reason):
    if not bool(x):raise ValueError("ALAN_NPMS_PREALAN: "+str(reason))

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def shard_years(shard):
    require(shard in range(SHARDS),"Unknown genuine NPMS observation year shard")
    return YEAR_RANGE[shard::SHARDS]

def site_year_source(package):
    with zipfile.ZipFile(package) as z:
        require(z.testzip() is None,"Original NPMS EIDC archive CRC corrupted")
        for member in ("data/sampleinfowithlatlong_2015to2025.csv",
                       "data/locations_2015to2025.csv"):
            src=z.getinfo(member)
            expected_size,expected_crc=NPMS_MEMBERS[member]
            require(src.file_size==expected_size and
                    format(src.CRC,"08x")==expected_crc,
                    "Original NPMS measured site/year or geolocation source changed: "+member)
        with z.open("data/sampleinfowithlatlong_2015to2025.csv") as stream:
            sample=pd.read_csv(stream,dtype=str,keep_default_na=False)
        with z.open("data/locations_2015to2025.csv") as stream:
            locations=pd.read_csv(stream,dtype=str,keep_default_na=False)
    require(len(sample)==27611 and len(locations)==5360 and
            locations.location_id.is_unique and sample.id.is_unique,
            "Original NPMS known visit and plot numbers no longer match")
    sample["date"]=pd.to_datetime(sample.date_start,errors="coerce",utc=True)
    sample["year"]=sample.date.dt.year
    sample["lat"]=pd.to_numeric(sample.LATITUDE,errors="coerce")
    sample["lon"]=pd.to_numeric(sample.LONGITUDE,errors="coerce")
    valid=sample.location_id.isin(set(locations.location_id)) & sample.year.between(
        2015,2025)&sample.lat.between(-90,90)&sample.lon.between(-180,180)
    sample=sample.loc[valid].copy()
    require(len(sample)==27603 and sample.location_id.nunique()==5360,
            "Original NPMS georeferenced visit source changed")
    original_site=sample.groupby("location_id",sort=False).agg(
        fixed_lon=("lon","median"),
        fixed_lat=("lat","median"),
        n_source_visits=("id","size"),
        n_original_observation_years=("year","nunique"))
    coords=sample.join(original_site[["fixed_lon","fixed_lat"]],on="location_id")
    _,_,meters=GEOD.inv(coords.fixed_lon.to_numpy(),coords.fixed_lat.to_numpy(),
                        coords.lon.to_numpy(),coords.lat.to_numpy())
    coords["source_coordinate_deviation_metres"]=np.asarray(meters,dtype=float)
    drift=coords.groupby("location_id").source_coordinate_deviation_metres.max()
    original_site["largest_original_visit_coordinate_shift_metres"]=drift
    original_site["site_is_spatially_fixed_within_500metres"]=(
        drift.le(MAX_SITE_DRIFT_M))
    sample=sample.merge(locations[["location_id","location_type"]],
        on="location_id",how="left",validate="many_to_one")
    require(sample.location_type.notna().all(),"Plot type missing from EIDC original")
    sample["month"]=sample.date.dt.month.astype(int)
    sample["original_survey_ID"]=sample.survey_id.astype(str).str.strip()
    sample["season"]=np.select(
        [sample.month.isin([4,5,6]),sample.month.isin([7,8,9])],
        ["AprJun","JulSep"],default="Other")
    year=sample.groupby(["location_id","year"],sort=False).agg(
        n_actual_source_survey_visits=("id","nunique"),
        n_distinct_original_NPMS_survey_codes=("original_survey_ID","nunique"),
        n_distinct_observed_season_windows=("season","nunique"),
        first_observed_month=("month","min"),
        last_observed_month=("month","max"),
        original_source_plot_type=("location_type","first")).reset_index()
    year=year.join(original_site[["fixed_lon","fixed_lat",
        "n_original_observation_years",
        "largest_original_visit_coordinate_shift_metres",
        "site_is_spatially_fixed_within_500metres"]],on="location_id")
    require(len(year)==16960 and not year.duplicated(["location_id","year"]).any() and
            set(year.year)==set(YEAR_RANGE),
            "Source original NPMS plot×year 16,960 row frame not preserved")
    year["prior_PCNL_year"]=year.year-1
    return year

def raster_sample(raster,rows):
    require(raster.count==1 and raster.crs is not None and
            (raster.crs.is_geographic or raster.crs.is_projected),
            "Original PCNL source raster unsupported")
    xx=rows.fixed_lon.to_numpy(float)
    yy=rows.fixed_lat.to_numpy(float)
    if raster.crs.to_epsg()!=4326:
        xx,yy=transform("EPSG:4326",raster.crs,xx.tolist(),yy.tolist())
    values=[]
    for a in raster.sample(list(zip(xx,yy)),indexes=1,masked=True):
        if np.ma.is_masked(a[0]):values.append(np.nan)
        else:
            v=float(a[0])
            values.append(v if np.isfinite(v) and v>=0 else np.nan)
    return np.asarray(values)

def worker(npms,shard,out):
    out=Path(out);require(not out.exists(),"Do not overwrite original year exposure")
    src=site_year_source(npms)
    yrs=shard_years(shard)
    q=src.loc[src.year.isin(yrs)].copy()
    require(len(q)>0 and set(q.year)==set(yrs),"Original observation years missing in source shard")
    original=manifest_read(Path("config/alan_pcnl_1992_2024_zenodo_v3_md5.json"))["files"]
    out.mkdir(parents=True)
    info=[]
    q["source_true_prior_PCNl"]=np.nan
    for observation_year in yrs:
        pre=observation_year-1
        name=f"PCNL{pre}.tif"
        real=fetch_one(name,original[name],out/name)
        active=q.year.eq(observation_year)&q.site_is_spatially_fixed_within_500metres
        s=q.loc[active].copy()
        with rasterio.open(out/name) as raster:
            measured=raster_sample(raster,s)
            geo={"src_crs":str(raster.crs),"native_resolution":list(raster.res),
                 "original_nodata":raster.nodata}
        q.loc[active,"source_true_prior_PCNl"]=measured
        info.append({
            "real_original_NPMS_observation_year":observation_year,
            "real_original_survey_plot_year_rows":int(q.year.eq(observation_year).sum()),
            "actual_prior_PCNL_image_year":pre,
            "original_source_plot_years_measured":int(np.isfinite(measured).sum()),
            "original_site_coordinate_drift_excluded":int((q.year.eq(observation_year)&
                ~q.site_is_spatially_fixed_within_500metres).sum()),
            "original_PCNl_MD5_source_receipt":real,
            "original_raster_geometry":geo})
        print("ALAN_NPMS_ORIGINAL_PRIOR_LIGHT_YEAR",json.dumps(info[-1],sort_keys=True),
              flush=True)
        (out/name).unlink()
    q["source_prior_log1p_PCNL"]=np.log1p(q.source_true_prior_PCNl)
    # Do not export coordinates at individual NPMS plot resolution.
    q=q.drop(columns=["fixed_lon","fixed_lat"])
    csv=out/f"NPMS_2015_2025_actual_plot_previous_PCNL_shard_{shard}.csv.gz"
    q.to_csv(csv,index=False,compression="gzip")
    audit={"schema":SCHEMA,"shard":shard,
       "real_source_NPMS_visit_years":yrs,
       "source_original_NPMS_plot_years":len(q),
       "source_valid_prior_PCNL_plot_years":int(q.source_prior_log1p_PCNL.notna().sum()),
       "source_original_raster_MD5_and_geodata":info,
       "source_CSV_SHA256":sha256(csv),
       "source_fixed_plot_not_photo_detection_zeros":True,
       "no_future_radiance_for_observed_plants":True}
    (out/f"NPMS_2015_2025_actual_plot_previous_PCNL_shard_{shard}.json").write_text(
        json.dumps(audit,sort_keys=True,indent=2)+"\n")
    return audit

def collect(npms,partials,out):
    out=Path(out);require(not out.exists(),"No overwrite original NPMS source-light cohort")
    original=site_year_source(npms)
    frames=[];verified=[]
    for sh in range(SHARDS):
        p=Path(partials)/f"NPMS_2015_2025_actual_plot_previous_PCNL_shard_{sh}.csv.gz"
        q=Path(partials)/f"NPMS_2015_2025_actual_plot_previous_PCNL_shard_{sh}.json"
        j=json.loads(q.read_text())
        require(j["schema"]==SCHEMA and j["shard"]==sh and
                j["real_source_NPMS_visit_years"]==shard_years(sh) and
                sha256(p)==j["source_CSV_SHA256"],
                "One actual original NPMS historical light source shard is unsupported")
        d=pd.read_csv(p,dtype={"location_id":str})
        require(len(d)==j["source_original_NPMS_plot_years"] and
                (d.prior_PCNL_year==d.year-1).all(),
                "Year-before source exposure not matched to plot-year")
        frames.append(d);verified.extend(j["real_source_NPMS_visit_years"])
    require(sorted(verified)==YEAR_RANGE,"Not all 2015-25 NPMS survey years present")
    frame=pd.concat(frames,ignore_index=True)
    keys=["location_id","year"]
    require(len(frame)==16960 and not frame.duplicated(keys).any(),
            "Original 16,960 source plot-year count failed after union")
    check=original[keys].merge(frame[keys],on=keys,how="left",
                               validate="one_to_one",indicator=True)
    require(check._merge.eq("both").all(),
            "Original repeated NPMS plot/year lost in annual PCNL join")
    out.mkdir(parents=True)
    csv=out/"NPMS_16960_real_plotyear_PCNl_prior2014to2024.csv.gz"
    frame.to_csv(csv,index=False,compression="gzip")
    observed=frame.loc[frame.source_prior_log1p_PCNL.notna()].copy()
    site_history=observed.groupby("location_id").agg(
         distinct_observation_years=("year","nunique"),
         min_prior_light=("source_prior_log1p_PCNL","min"),
         max_prior_light=("source_prior_log1p_PCNL","max"),
         n_unique_light_values=("source_prior_log1p_PCNL","nunique"))
    active=site_history.loc[site_history.distinct_observation_years.ge(2)]
    span=active.max_prior_light-active.min_prior_light
    info={
      "schema":SCHEMA,"original_NPMS_DOI":DOI,
      "corrected_PCNL_annual_original_DOI":PCNL_DOI,
      "genuine_source_2015_2025_original_plot_years":16960,
      "genuine_source_distinct_original_NPMS_plots":5360,
      "genuine_source_plot_years_valid_prior_light":len(observed),
      "genuine_source_independent_plots_with_valid_prior_light":int(observed.location_id.nunique()),
      "genuine_source_plots_revisited_two_or_more_years_with_valid_annual_PCNL":len(active),
      "same_plot_observed_years_at_least_2_different_satellite_light_values":int(
          active.n_unique_light_values.ge(2).sum()),
      "same_plot_actual_log1p_prior_PCNL_range_at_least_0p02":int(span.ge(.02).sum()),
      "same_plot_actual_log1p_prior_PCNL_range_at_least_0p1":int(span.ge(.1).sum()),
      "same_plot_actual_log1p_prior_PCNL_range_at_least_0p3":int(span.ge(.3).sum()),
      "same_plot_true_light_range_quantiles":{
          str(p):float(span.quantile(p)) for p in (.1,.25,.5,.75,.9,.95)},
      "actual_survey_years":YEAR_RANGE,
      "actual_prior_satellite_years":[y-1 for y in YEAR_RANGE],
      "plot_fixed_location_coordinate_shift_gate_met_500m":True,
      "original_plotyear_radiance_source_CSV_SHA256":sha256(csv),
      "not_moth_night_pollination_or_plant_cover_response":True,
      "no_unreported_species_zeros_or_plot_absence_created":True,
      "NPMS_survey_method_and_season_not_yet_enforced_in_light_only_frame":True}
    (out/"NPMS_priorPCNL_2014_2024_actual_repeated_plot_source_receipt.json").write_text(
        json.dumps(info,sort_keys=True,indent=2)+"\n")
    print("ALAN_NPMS_TRUE_11YEAR_WITHIN_PLOT_HISTORICAL_LIGHT_VARIATION",
          json.dumps(info,sort_keys=True),flush=True)
    return info

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="action",required=True)
    w=sub.add_parser("worker")
    w.add_argument("--npms",type=Path,required=True)
    w.add_argument("--shard",type=int,required=True)
    w.add_argument("--out",type=Path,required=True)
    c=sub.add_parser("collect")
    c.add_argument("--npms",type=Path,required=True)
    c.add_argument("--partials",type=Path,required=True)
    c.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.action=="worker": worker(args.npms,args.shard,args.out)
    else: collect(args.npms,args.partials,args.out)

if __name__=="__main__":main()
