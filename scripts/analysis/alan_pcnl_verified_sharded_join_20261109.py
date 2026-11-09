#!/usr/bin/env python3
"""Join actual SHA-pinned 93m iNaturalist island flora source with PCNL 1992-2024.

Six mutually exclusive annual-PCNL shards. Each downloads/MD5-verifies only
its own years, extracts actual cross-sensor harmonized satellite pixels at
the original 0.05-degree photo-cell centers, and exports additive partial sums.
A separate collector requires all 33 years and all six receipts; all prior
years to the observation must be present or exposure is NA. No flower ROI,
no species absence and no biological inference in the join.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.acquisition.alan_pcnl_1992_2024_verified_zenodo_download_20261109 import (
    fetch_one, manifest_read,
)
SOURCE_SHA="e7800476ae6d35acb9ef71525f32521596f884733e2e8ee9f87fda7acb8e9eab"
NROWS=85318
SHARDS=6
YEARS=list(range(1992,2025))
SCHEMA="ALAN_real_full93m_island_PCNL33_year_6shard_exposure_v1"
INPUT_COLUMNS={
 "island","latbin","lonbin","year","lat","lon",
 "all_accepted_recorded_species","moth_visitor_positive_recorded_species",
 "all_accepted_photo_ids","all_accepted_observation_uuids",
 "moth_positive_fraction_OF_RECORDED_WCVP_SPECIES",
}

def require(ok,msg):
    if not ok: raise ValueError("ALAN_SHARDED_PCNL: "+str(msg))

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def source_data(path):
    path=Path(path)
    require(sha256(path)==SOURCE_SHA,"Original full93m source CSV changed")
    x=pd.read_csv(path)
    require(len(x)==NROWS and INPUT_COLUMNS.issubset(x.columns),
            "Original source site-year shape or columns changed")
    require(not x[["island","latbin","lonbin","year"]].duplicated().any(),
            "Duplicate island/site/year")
    require(x.year.isin([2019,2021,2023,2024,2025]).all(),
            "Wrong photographic years")
    require((x.all_accepted_recorded_species>=3).all()
            and (x.all_accepted_observation_uuids>=5).all(),
            "Source photo opportunity prespecification changed")
    require((x.moth_visitor_positive_recorded_species <=
             x.all_accepted_recorded_species).all(),
            "Moth-positive source numerator exceeds total documented species")
    ratio=x.moth_visitor_positive_recorded_species/x.all_accepted_recorded_species
    require(np.allclose(ratio, x.moth_positive_fraction_OF_RECORDED_WCVP_SPECIES,
                        atol=1e-9,rtol=0),
            "Source fraction no longer exact ratio of accepted species")
    return x

def split_years(shard):
    require(shard in range(SHARDS),"Bad shard")
    years=[yr for yr in YEARS if (yr-1992)%SHARDS==shard]
    require(5<=len(years)<=6 and all(y in YEARS for y in years),"Bad years")
    return years

def sample_one_raster(src,x):
    """Raster-native sampling with reprojection; masked/no-data stays absent."""
    from rasterio.warp import transform
    coords_x=x.lon.to_numpy(float)
    coords_y=x.lat.to_numpy(float)
    if src.crs.to_epsg()!=4326:
        coords_x,coords_y=transform("EPSG:4326",src.crs,
                                   coords_x.tolist(),coords_y.tolist())
    length=len(x)
    vals=np.full(length,np.nan,dtype=float)
    for pos in range(0,length,12000):
        subset=list(zip(coords_x[pos:pos+12000],coords_y[pos:pos+12000]))
        raw=src.sample(subset,indexes=1,masked=True)
        batch=np.asarray([np.nan if np.ma.is_masked(v[0]) else float(v[0])
                          for v in raw],dtype=float)
        good=np.isfinite(batch)&(batch>=0)
        vals[pos:pos+len(batch)]=np.where(good,batch,np.nan)
    return vals

def worker(source,manifest,shard,raster_dir,out):
    import rasterio
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite partial exposure")
    out.mkdir(parents=True)
    x=source_data(source)
    files=manifest_read(manifest)["files"]
    observed=x.year.to_numpy(int)
    partial=np.zeros(len(x),dtype=float)
    nvalid=np.zeros(len(x),dtype=np.int16)
    years=split_years(shard)
    audit=[]
    Path(raster_dir).mkdir(parents=True,exist_ok=True)
    for yr in years:
        name=f"PCNL{yr}.tif"
        f=fetch_one(name,files[name],Path(raster_dir)/name)
        with rasterio.open(Path(raster_dir)/name) as src:
            require(src.count==1 and src.crs is not None and src.nodata is not None,
                    "PCNL raster source not spatially valid")
            values=sample_one_raster(src,x)
        eligible=(observed>yr)
        good=eligible & np.isfinite(values)
        partial[good]+=np.log1p(values[good])
        nvalid[good]+=1
        audit.append({"year":yr,**f,
                      "n_before_photo_year":int(eligible.sum()),
                      "n_source_pixels_valid":int(good.sum())})
        print("ALAN_PCNL_YEAR_GENUINE_PIXEL_SOURCE",
              json.dumps(audit[-1],sort_keys=True),flush=True)
    d=pd.DataFrame({"row":np.arange(len(x),dtype=np.int32),
                    "partial_sum":partial,
                    "partial_n":nvalid})
    path=out/f"pcnl_shard_{shard}.csv.gz"
    d.to_csv(path,index=False,compression="gzip")
    receipt={"schema":SCHEMA,"shard":shard,
             "source_site_year_csv_sha256":SOURCE_SHA,
             "n_source_site_year_rows":len(x),"years":years,
             "raster_audit":audit,
             "partial_csv_sha256":sha256(path),
             "no_stored_original_RGB_or_floral_ROI":True}
    (out/f"pcnl_shard_{shard}.json").write_text(json.dumps(receipt,indent=2,
                                                        sort_keys=True)+"\n")
    return receipt

def collect(source,parts,out):
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite completed 33-year exposure")
    x=source_data(source)
    partial_sum=np.zeros(len(x),dtype=float)
    partial_n=np.zeros(len(x),dtype=np.int16)
    seen=[]
    all_rasters=[]
    for shard in range(SHARDS):
        partial_file=Path(parts)/f"pcnl_shard_{shard}.csv.gz"
        receipt_file=Path(parts)/f"pcnl_shard_{shard}.json"
        receipt=json.loads(receipt_file.read_text())
        require(receipt["schema"]==SCHEMA and receipt["shard"]==shard and
                receipt["source_site_year_csv_sha256"]==SOURCE_SHA and
                receipt["n_source_site_year_rows"]==NROWS and
                receipt["years"]==split_years(shard) and
                sha256(partial_file)==receipt["partial_csv_sha256"],
                "One true original PCNL annual shard is unverified")
        d=pd.read_csv(partial_file)
        require(len(d)==NROWS and
                np.array_equal(d.row.to_numpy(),np.arange(NROWS)),
                "Original row key misaligned across yearly PCNL shards")
        require(np.isfinite(d.partial_sum.to_numpy()).all() and
                (d.partial_n.to_numpy()>=0).all(),
                "Invalid source partial exposure sum/count")
        partial_sum+=d.partial_sum.to_numpy(float)
        partial_n+=d.partial_n.to_numpy(np.int16)
        seen.extend(receipt["years"])
        all_rasters.extend(receipt["raster_audit"])
    require(sorted(seen)==YEARS and len(set(seen))==33,
            "Wrong or incomplete 1992–2024 PCNL annual source")
    expected=x.year.to_numpy(int)-1992
    require((partial_n<=expected).all(),
            "Future raster year included or double counted")
    complete=(partial_n==expected)
    x["pcnl_years_expected"]=expected
    x["pcnl_years_observed"]=partial_n
    x["pcnl_1992_to_previous_year_sum_log1p"]=np.where(
        complete,partial_sum,np.nan)
    x["pcnl_1992_to_previous_year_mean_log1p"]=np.where(
        complete,partial_sum/expected,np.nan)
    x["historic_ntl_complete"]=complete
    out.mkdir(parents=True)
    site=out/"ALAN_real_93m_flora_island_2019_2025_PCNLPrior1992_2024.csv.gz"
    x.to_csv(site,index=False,compression="gzip")
    summary={"schema":"ALAN_full93m_real_PCNL_33year_1992_2024_source_join_v1",
      "source_original_photo_IDs":93209351,
      "source_site_year_rows":len(x),
      "source_unique_islands":int(x.island.nunique()),
      "pcnl_complete_prior_year_site_years":int(complete.sum()),
      "pcnl_complete_prior_year_islands":int(x.loc[complete,"island"].nunique()),
      "year_by_year_complete":{str(y):int((complete &
                (x.year.to_numpy(int)==y)).sum()) for y in sorted(x.year.unique())},
      "output_sha256":sha256(site),
      "original_source_sha256":SOURCE_SHA,
      "raster_years_verified":33,
      "raster_bytes_verified":int(sum(d["size_bytes"] for d in all_rasters)),
      "raster_MD5s_from_official_Zenodo_v3":{str(d["year"]):d["md5"]
                                                for d in all_rasters},
      "no_biological_moth_nonvisitor_or_species_absence_inferred":True,
      "observation_year_not_exposed_to_future_nightlight":True}
    (out/"real_pcnl_join_receipt.json").write_text(json.dumps(
        summary,indent=2,sort_keys=True)+"\n")
    print("ALAN_REAL_PCNL_FULL93M_MOTH_SOURCE_JOIN",
          json.dumps(summary,sort_keys=True),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="mode",required=True)
    w=sub.add_parser("worker")
    w.add_argument("--source",required=True,type=Path)
    w.add_argument("--manifest",required=True,type=Path)
    w.add_argument("--shard",required=True,type=int)
    w.add_argument("--raster-dir",required=True,type=Path)
    w.add_argument("--out",required=True,type=Path)
    c=sub.add_parser("collect")
    c.add_argument("--source",required=True,type=Path)
    c.add_argument("--parts",required=True,type=Path)
    c.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    if args.mode=="worker":
        worker(args.source,args.manifest,args.shard,args.raster_dir,args.out)
    else: collect(args.source,args.parts,args.out)

if __name__=="__main__":main()
