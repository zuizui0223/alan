#!/usr/bin/env python3
"""MD5-verified 30 annual Li/Cao PCNL rasters sampled at original day/night
fruit/seed experimental site coordinates with time locked to true experiment
year_start. 6 disjoint year groups, 1992–2021, 5 annual raster each.

Do NOT substitute study publication year, assume same site if source coordinates
differ, or use any post-experiment annual PCNL. Source experimental outcome is
not used to choose geographic exposure or match-pixel availability.

The source consists of 301 site/time/metric independently paired experiment
summaries, 118 original studies. These are NOT independent 301 studies.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one,manifest_read
)
from scripts.analysis.alan_pcnl_verified_sharded_join_20261109 import sample_one_raster

START=1992
LAST=2021
YEARS=list(range(START,LAST+1))
SHARDS=6
EXPECTED_SHA="21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af"
N_PAIRS=301
SCHEMA="ALAN_original_DvN_301_true_study_site_year_1992_2021_precursor_PCNL_6shards_v1"

def require(ok,message):
    if not ok:raise ValueError("ALAN_TRUE_FIELD_DIEL_PRIOR_LIGHT: "+str(message))
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def source_read(src,receipt):
    r=json.loads(Path(receipt).read_text())
    require(r["schema"]=="ALAN_original_DvN_source_site_year_pre_PCNLP_v1" and
            r["actual_original_site_year_reproductive_pairs"]==N_PAIRS and
            r["n_independent_original_studies"]==118 and
            r["original_outcome_blind_satellite_source_sha256"]==EXPECTED_SHA and
            sha(src)==EXPECTED_SHA,
            "Authenticated original 301 source field-site experiment cohort drift")
    x=pd.read_csv(src,dtype={"study_key":str,"species":str})
    require(len(x)==N_PAIRS and x.source_pair_id.is_unique and
            x.year_exp.between(1993,2022).all() and
            x[["lat","lon"]].notna().all().all() and
            x.night_share.between(0,1).all(),
            "Source field site/experimental outcome integrity failed")
    return x

def years(shard):
    require(shard in range(SHARDS),"Unknown source disjoint year shard")
    return [y for y in YEARS if (y-START)%SHARDS==shard]

def worker(source,receipt,manifest,shard,rasterdir,out):
    import rasterio
    x=source_read(source,receipt)
    files=manifest_read(manifest)["files"]
    out=Path(out)
    require(not out.exists(),"Refuse overwrite of source raster annual shard")
    out.mkdir(parents=True)
    year_list=years(shard)
    sums=np.zeros(len(x),dtype=float)
    counts=np.zeros(len(x),dtype=np.int16)
    raster_receipts=[]
    Path(rasterdir).mkdir(parents=True,exist_ok=True)
    for year in year_list:
        fn=f"PCNL{year}.tif"
        rr=fetch_one(fn,files[fn],Path(rasterdir)/fn)
        with rasterio.open(Path(rasterdir)/fn) as tif:
            require(tif.count==1 and tif.crs is not None,
                    "Raster source missing original geo-reference")
            val=sample_one_raster(tif,x)
        before=x.year_exp.to_numpy(int)>year
        valid=before & np.isfinite(val) & (val>=0)
        sums[valid]+=np.log1p(val[valid])
        counts[valid]+=1
        raster_receipts.append({
          "raster_year":year,**rr,
          "field_experiments_strictly_after_this_raster_year":int(before.sum()),
          "source_fieldpoints_with_valid_original_pixels":int(valid.sum())
        })
        print("ALAN_REAL_NIGHT_DIEL_YEAR_RASTER_SAMPLE",json.dumps(
            raster_receipts[-1],sort_keys=True),flush=True)
    dest=out/f"diel_annual_shard_{shard}.csv.gz"
    pd.DataFrame({
       "source_pair_id":x.source_pair_id,
       "sum_log1p_PCNL_before_source_experiment":sums,
       "n_years_original_PCNL_before_experiment":counts
    }).to_csv(dest,index=False,compression="gzip")
    (out/f"diel_annual_shard_{shard}.json").write_text(json.dumps({
        "schema":SCHEMA,"shard":shard,"raster_years":year_list,
        "source_sha":EXPECTED_SHA,"source_pairs":len(x),
        "source_file_sha256":sha(dest),"authenticated_PCNL_raster_years":raster_receipts,
        "PCNL_strictly_before_start_of_field_experiment":True
    },indent=2,sort_keys=True)+"\n")

def collect(source,receipt,partialdir,out):
    x=source_read(source,receipt)
    totals=np.zeros(len(x),float)
    counts=np.zeros(len(x),np.int16)
    all_years=[]
    all_receipts=[]
    for shard in range(SHARDS):
        meta=json.loads((Path(partialdir)/f"diel_annual_shard_{shard}.json").read_text())
        p=Path(partialdir)/f"diel_annual_shard_{shard}.csv.gz"
        require(meta["schema"]==SCHEMA and meta["shard"]==shard and
                meta["raster_years"]==years(shard) and
                meta["source_sha"]==EXPECTED_SHA and
                meta["source_file_sha256"]==sha(p),
                "Original PCNL year shard source identity mismatch")
        d=pd.read_csv(p)
        require(len(d)==len(x) and
                np.array_equal(d.source_pair_id.to_numpy(),x.source_pair_id.to_numpy()),
                "Original source 301 pair IDs not aligned across annual years")
        totals+=d.sum_log1p_PCNL_before_source_experiment.to_numpy(float)
        counts+=d.n_years_original_PCNL_before_experiment.to_numpy(np.int16)
        all_years.extend(meta["raster_years"])
        all_receipts.extend(meta["authenticated_PCNL_raster_years"])
    require(sorted(all_years)==YEARS and len(set(all_years))==30,
            "Not every original 1992–2021 annual source raster sampled exactly once")
    expected=x.year_exp.to_numpy(int)-START
    require(np.all(counts<=expected) and expected.min()>=1,
            "Incomplete chronology or satellite exposure sampled after study")
    full=counts==expected
    x["pcnl_historical_complete"]=full
    x["pcnl_prior_expected_years"]=expected
    x["pcnl_prior_actual_sampled_years"]=counts
    x["pcnl_prior_mean_log1p"]=np.where(full,totals/expected,np.nan)
    x["pcnl_prior_sum_log1p"]=np.where(full,totals,np.nan)
    out=Path(out)
    require(not out.exists(),"Refuse rewriting signed source observations")
    out.mkdir(parents=True)
    dest=out/"Diel_original_daynight_reproductive_1993_2022_priorPCNL1992_2021.csv.gz"
    x.to_csv(dest,index=False,compression="gzip")
    good=x.loc[full]
    meta={
       "schema":"ALAN_DvN_reproductive_real_priorPCNL_source_join_v1",
       "original_source_pairs":N_PAIRS,
       "original_source_sha256":EXPECTED_SHA,
       "actual_raster_years":YEARS,
       "annual_rasters_MD5_source_verified":len(all_receipts),
       "source_pairs_full_pre_experiment_PCNL":len(good),
       "independent_published_studies_full_PCNLP":int(good.study_key.nunique()),
       "original_field_site_georeferenced_points_PCNL":int(
           good[["lat","lon"]].drop_duplicates().shape[0]),
       "original_field_study_site_years_PCNL":int(
           good[["study_key","lat","lon","year_exp"]].drop_duplicates().shape[0]),
       "original_field_plant_species_PCNL":int(good.species.nunique()),
       "original_field_experiment_start_years":[int(x.year_exp.min()),int(x.year_exp.max())],
       "fraction_zero_prior_light":float(good.pcnl_prior_mean_log1p.eq(0).mean()) if len(good) else None,
       "sd_original_mean_prior_light":float(good.pcnl_prior_mean_log1p.std()) if len(good) else None,
       "joined_file_sha256":sha(dest),
       "true_source_geography_year_prior_light_only":True,
       "no_future_year_or_1992_study_without_prior_light":True,
       "outcome_direction_NOT_used_to_choose_raster_cohort":True,
       "ALAN_causation_of_real_pollen_deposition_NOT_established":True,
       "author_experiment_site_level_night_function_not_moth_specific":True
    }
    (out/"Diel_original_priorPCNL_source_receipt.json").write_text(
        json.dumps(meta,indent=2,sort_keys=True)+"\n")
    print("ALAN_AUTHENTIC_ORIGINAL_DIEL_DAYNIGHT_FRUITSEED_PCNLPrior_SITEYEAR",
          json.dumps(meta,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="task",required=True)
    w=sub.add_parser("worker")
    for key in ("source","receipt","manifest","rasterdir","out"):
        w.add_argument("--"+key,type=Path,required=True)
    w.add_argument("--shard",type=int,required=True)
    c=sub.add_parser("collect")
    for key in ("source","receipt","partialdir","out"):
        c.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    if a.task=="worker":
        worker(a.source,a.receipt,a.manifest,a.shard,a.rasterdir,a.out)
    else:
        collect(a.source,a.receipt,a.partialdir,a.out)
