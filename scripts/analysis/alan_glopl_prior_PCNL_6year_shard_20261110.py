#!/usr/bin/env python3
"""Link ORIGINAL GloPL experiment-year/coordinate source to PRIOR PCNL 1992–2014.

This uses 6 disjoint annual-pixel-shard workers and one single aggregator.
The original GloPL cohort is outcome-blind with respect to ALAN selection:
original GloPL source year_of_study, coordinates, species and independent
publication identity were frozen BEFORE any PCNL values or outcomes were read.

Original GloPL supplementation log response ratios are NOT plant abundance,
moth pollination identity, local NPMS plot measurements or a causal ALAN effect.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one, manifest_read
)
from scripts.analysis.alan_pcnl_verified_sharded_join_20261109 import sample_one_raster

SOURCE_SHA="ac8e8d191edbe175fa0f72cb61387eece986d31487ae2db9cdf3619c58d80930"
SOURCE_ROWS=1419
YEARS=list(range(1992,2015))
SHARDS=6
SCHEMA="ALAN_real_GloPL_experimental_pollen_limitation_priorPCNL_six_year_shards_v1"

def require(ok,msg):
    if not ok:raise ValueError("ALAN_GLOPL_PCNL: "+str(msg))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def source_read(path):
    require(sha(path)==SOURCE_SHA,"Wrong original outcome-blind GloPL experimental cohort")
    x=pd.read_csv(path,dtype={"study_key":str,"species":str})
    require(len(x)==SOURCE_ROWS and
           {"original_row_id","year_exp","lat","lon","study_key","PL_Effect_Size",
            "PL_Effect_Size_Type1","PL_Effect_Size_Type2",
            "Constant_added","Level_of_Supplementation"}<=set(x),
            "Original 1419 GloPL experimental rows or design fields changed")
    require(x.original_row_id.is_unique and
            x.year_exp.between(1993,2015).all() and
            x.lat.between(-90,90).all() and x.lon.between(-180,180).all() and
            x.study_key.notna().all() and x.study_key.astype(str).str.strip().ne("").all(),
            "Invalid year, coordinate, row key, or publication source")
    return x

def years(shard):
    require(shard in range(SHARDS),"Not one of six frozen PCNL year groups")
    return [y for y in YEARS if (y-1992)%SHARDS==shard]

def worker(source,manifest,shard,raster_dir,out):
    import rasterio
    x=source_read(source)
    files=manifest_read(manifest)["files"]
    out=Path(out)
    require(not out.exists(),"Previous original source shard must not be overwritten")
    out.mkdir(parents=True)
    year_indices=years(shard)
    light_sum=np.zeros(len(x),dtype=float)
    light_n=np.zeros(len(x),dtype=np.int16)
    audits=[]
    Path(raster_dir).mkdir(parents=True,exist_ok=True)
    for y in year_indices:
        name=f"PCNL{y}.tif"
        source_receipt=fetch_one(name,files[name],Path(raster_dir)/name)
        with rasterio.open(Path(raster_dir)/name) as tif:
            require(tif.count==1 and tif.crs is not None,
                    "Original PCNL single-band year geolocation unavailable")
            values=sample_one_raster(tif,x)
        need=x.year_exp.to_numpy(int)>y
        valid=need & np.isfinite(values)&(values>=0)
        light_sum[valid]+=np.log1p(values[valid])
        light_n[valid]+=1
        audits.append({"year":y,**source_receipt,
                       "n_experiments_with_prior_year":int(need.sum()),
                       "n_valid_actual_satellite_pixels":int(valid.sum())})
        print("ALAN_REAL_GLOPL_PCNL_EXPERIMENT_PRIOR_YEAR",
              json.dumps(audits[-1],sort_keys=True),flush=True)
    outfile=out/f"glopl_pcnl_shard_{shard}.csv.gz"
    pd.DataFrame({"original_row_id":x.original_row_id,
                  "partial_logp_sums":light_sum,
                  "partial_valid_years":light_n}).to_csv(
                       outfile,index=False,compression="gzip")
    (out/f"glopl_pcnl_shard_{shard}.json").write_text(json.dumps({
        "schema":SCHEMA,"shard":shard,"year_group":year_indices,
        "input_source_sha256":SOURCE_SHA,"source_experiments":len(x),
        "raster_originals_verified":audits,"partial_sha256":sha(outfile),
        "historic_light_only_before_experiment":True},sort_keys=True,indent=2)+"\n")

def collect(source,partial_dir,out):
    x=source_read(source)
    year_sum=np.zeros(len(x),dtype=float)
    year_count=np.zeros(len(x),dtype=np.int16)
    all_years=[]
    receipts=[]
    for shard in range(SHARDS):
        path=Path(partial_dir)/f"glopl_pcnl_shard_{shard}.csv.gz"
        receipt=json.loads((Path(partial_dir)/f"glopl_pcnl_shard_{shard}.json").read_text())
        require(receipt["schema"]==SCHEMA and receipt["shard"]==shard and
                receipt["year_group"]==years(shard) and
                receipt["input_source_sha256"]==SOURCE_SHA and
                receipt["source_experiments"]==SOURCE_ROWS and
                sha(path)==receipt["partial_sha256"],
                "One independently verified source-shard is missing or changed")
        d=pd.read_csv(path)
        require(len(d)==len(x) and
                np.array_equal(d.original_row_id.to_numpy(),x.original_row_id.to_numpy()),
                "Original GloPL row key shifted between annual PCNL shards")
        year_sum+=d.partial_logp_sums.to_numpy(float)
        year_count+=d.partial_valid_years.to_numpy(np.int16)
        all_years.extend(receipt["year_group"])
        receipts.append(receipt)
    require(sorted(all_years)==YEARS and len(set(all_years))==len(YEARS),
            "Incomplete or duplicate annual nightlight source years")
    expected=x.year_exp.to_numpy(int)-1992
    require(np.all(year_count<=expected),"Future/duplicate PCNL light year included")
    full=year_count==expected
    x["pcnl_pre_experiment_expected_years"]=expected
    x["pcnl_pre_experiment_observed_years"]=year_count
    x["pcnl_prior_mean_log1p"]=np.where(full,year_sum/expected,np.nan)
    x["pcnl_prior_sum_log1p"]=np.where(full,year_sum,np.nan)
    x["pcnl_historical_complete"]=full
    dest=Path(out)
    require(not dest.exists(),"Refuse replacing genuine GloPL×ALAN result")
    dest.mkdir(parents=True)
    outfile=dest/"GloPL_experimental_PL1993_2015_priorPCNL1992_2014.csv.gz"
    x.to_csv(outfile,index=False,compression="gzip")
    audited=x.loc[x.pcnl_historical_complete].copy()
    receipt={
       "schema":"ALAN_real_GloPL_PCNLPrior1992_2014_source_join_v1",
       "original_outcome_blind_source_sha256":SOURCE_SHA,
       "source_experiment_rows":len(x),
       "complete_PCNLPrior_experiment_rows":len(audited),
       "complete_PCNLPrior_sites":len(audited[["lat","lon"]].drop_duplicates()),
       "complete_PCNLPrior_site_years":len(audited[["lat","lon","year_exp"]].drop_duplicates()),
       "complete_PCNLPrior_publications":int(audited.study_key.nunique()),
       "n_complete_site_years_with_zero_prior_PCNLP":int(audited.pcnl_prior_sum_log1p.eq(0).sum()),
       "complete_light_iqr":(audited.pcnl_prior_mean_log1p.quantile(.75)-
                              audited.pcnl_prior_mean_log1p.quantile(.25)
                              if len(audited) else None),
       "complete_light_range":([float(audited.pcnl_prior_mean_log1p.min()),
                                float(audited.pcnl_prior_mean_log1p.max())]
                                if len(audited) else []),
       "year_min":int(x.year_exp.min()),"year_max":int(x.year_exp.max()),
       "original_PCNLP_years_verified":23,"original_PCNLP_raster_total_bytes":int(
          sum(a["size_bytes"] for row in receipts for a in row["raster_originals_verified"])),
       "output_sha256":sha(outfile),
       "GloPL_PL_outcome_not_used_to_choose_source_cohort":True,
       "historical_light_never_future_experiment":True,
       "this_is_not_effect_of_ALAN_on_MOTH_DEPENDENCE_or_NPMS_plots":True,
       "ALAN_PL_regression_fitted":False}
    (dest/"original_GloPL_PCNLPrior_join_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("ALAN_REAL_GLOPL_EXPERIMENT_YEAR_HISTORIC_PCNL_SUPPORT",
          json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="command",required=True)
    w=sub.add_parser("worker")
    for k in ("source","manifest","raster-dir","out"):
        w.add_argument("--"+k,type=Path,required=True)
    w.add_argument("--shard",type=int,required=True)
    c=sub.add_parser("collect")
    for k in ("source","partials","out"):
        c.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    if a.command=="worker":
        worker(a.source,a.manifest,a.shard,a.raster_dir,a.out)
    else: collect(a.source,a.partials,a.out)

if __name__=="__main__":main()
