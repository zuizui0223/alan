#!/usr/bin/env python3
"""QA-verified wetland vegetation greenness coverage on original 187 ALAN pairs.

Never treat NDVI as wetland area. Never rematch based on botanical outcomes or NDVI.
Stop when both arms of the original matched pairs lack time-series coverage.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from alan_light_development_overlap_20261119 import read_source, pair_sites, SEED

ROOT=Path(__file__).resolve().parents[2]
NDVI=ROOT/"results/all_species_forest/wet_annual/vegetation.csv.gz"
HAB=ROOT/"results/all_species_forest/wet_annual/habitats.csv"
COAST=ROOT/"results/all_species_forest/wet_annual/coastal_screen.csv"
PIN_NDVI="26b18e38149d78d0914e5a89c1f863d766c93e4b"
PIN_HAB="6c464c2c2d2a9ddc71a95cab912175c263a90ff0"
TERMS=["light_difference","d_built","d_built0","d_light0",
       "d_temp_trend","d_rain_trend","d_freq0","d_area0"]
WATER=["d_freq_trend","d_area_trend"]
VEG=["d_ndvi_initial","d_ndvi_trend"]

def git_blob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()

def wetland_vegetation_sites(raw):
    need={"unit_id","y","day","qa","ndvi"}
    if need-set(raw):
        raise ValueError("SOURCE_NDVI_SCHEMA_STOP")
    x=raw.dropna(subset=list(need)).copy()
    x=x[(x.qa==0)&x.y.between(2000,2023)&x.ndvi.between(-.2,1.)]
    x["date"]=pd.to_datetime(x.y.astype(int).astype(str)+"-01-01")+pd.to_timedelta(x.day-1,unit="D")
    x["month"]=x.date.dt.month
    x=x[x.month.between(4,9)].copy()
    x["season"]=((x.month-4)//2).astype(int)
    s=x.groupby(["unit_id","y","season"]).ndvi.median().unstack("season")
    if s.empty:
        return pd.DataFrame(columns=["unit_id","ndvi_initial","ndvi_trend","valid_years"])
    for i in (0,1,2):
        if i not in s:
            s[i]=np.nan
    v=s[[0,1,2]].mean(axis=1).where(s[[0,1,2]].notna().all(axis=1)).rename("ndvi").dropna().reset_index()
    rows=[]
    for uid,g in v.groupby("unit_id"):
        early=g[g.y.between(2000,2004)].ndvi
        late=g[g.y.between(2019,2023)].ndvi
        if len(g)<12 or g.y.max()-g.y.min()<18 or len(early)<2 or len(late)<2:
            continue
        rows.append(dict(unit_id=int(uid),ndvi_initial=float(early.median()),
            ndvi_trend=float(np.polyfit(g.y.to_numpy(float),g.ndvi.to_numpy(float),1)[0]*10),
            valid_years=int(len(g))))
    return pd.DataFrame(rows,columns=["unit_id","ndvi_initial","ndvi_trend","valid_years"])

def fixed_vegetation_overlap(pairs,sites):
    if not {"bright_site_id","comparison_site_id"}<=set(pairs):
        raise ValueError("ORIGINAL_ARM_IDENTITIES_MISSING")
    if sites.unit_id.duplicated().any():
        raise ValueError("DUPLICATE_NDVI_SITE")
    s=sites[["unit_id","ndvi_initial","ndvi_trend"]]
    b=s.rename(columns={"unit_id":"bright_site_id","ndvi_initial":"b_initial","ndvi_trend":"b_trend"})
    c=s.rename(columns={"unit_id":"comparison_site_id","ndvi_initial":"c_initial","ndvi_trend":"c_trend"})
    z=pairs.merge(b,on="bright_site_id",how="left",validate="m:1")
    z=z.merge(c,on="comparison_site_id",how="left",validate="m:1")
    z["d_ndvi_initial"]=z.b_initial-z.c_initial
    z["d_ndvi_trend"]=z.b_trend-z.c_trend
    both=z.dropna(subset=VEG)
    coverage=dict(original_pairs=int(len(pairs)),original_blocks=int(pairs.block.nunique()),
        observed_wetland_vegetation_sites=int(len(sites)),
        bright_arm_ndvi_sites=int(z.b_initial.notna().sum()),
        comparison_arm_ndvi_sites=int(z.c_initial.notna().sum()),
        both_arms_ndvi_original_pairs=int(len(both)),
        both_arms_ndvi_blocks=int(both.block.nunique()))
    return both,coverage

def slope(data,columns):
    x=data[columns].to_numpy(float)
    if len(data)<len(columns)+6 or np.linalg.matrix_rank(x)<x.shape[1]:
        return np.nan
    return 10*float(np.linalg.lstsq(x,data.outcome_difference.to_numpy(float),rcond=None)[0][0])

def nested_on_same_pairs(pairs):
    x=pairs.dropna(subset=TERMS+WATER+VEG+["outcome_difference"]).copy()
    status=dict(same_subset_pairs=int(len(x)),blocks=int(x.block.nunique()))
    if len(x)<30 or x.block.nunique()<10 or len(x)<=len(TERMS+WATER+VEG)+5:
        return {**status,"status":"VEGETATION_SOURCE_SUPPORT_STOP"}
    fits=[slope(x,TERMS),slope(x,TERMS+WATER),slope(x,TERMS+WATER+VEG)]
    if not np.isfinite(fits).all():
        return {**status,"status":"VEGETATION_DESIGN_RANK_STOP"}
    groups=[g for _,g in x.groupby("block")]
    rng=np.random.default_rng(SEED+2)
    draws=[]
    for _ in range(1500):
        bs=pd.concat([groups[i] for i in rng.integers(0,len(groups),len(groups))],ignore_index=True)
        draws.append([slope(bs,TERMS),slope(bs,TERMS+WATER),slope(bs,TERMS+WATER+VEG)])
    draws=np.asarray(draws)
    draws=draws[np.isfinite(draws).all(axis=1)]
    if len(draws)<500:
        return {**status,"status":"VEGETATION_BOOTSTRAP_STOP"}
    ci=np.quantile(draws,[.025,.975],axis=0)
    return {**status,"status":"EXPLORATORY_NESTED_FITS_ONLY",
      "slopes_per_10_DN_pp_per_decade":dict(
         without_water=fits[0],with_water=fits[1],with_water_and_vegetation=fits[2]),
      "block_bootstrap_95pct":dict(without_water=ci[:,0].tolist(),
         with_water=ci[:,1].tolist(),with_water_and_vegetation=ci[:,2].tolist()),
      "valid_block_bootstrap_draws":int(len(draws))}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",default="")
    args=ap.parse_args()
    if git_blob(NDVI)!=PIN_NDVI or git_blob(HAB)!=PIN_HAB:
        raise RuntimeError("PINNED_ORIGINAL_VEGETATION_SOURCE_CHANGED")
    raw,sha=read_source()
    pairs,_,_=pair_sites(raw,False)
    if len(pairs)!=187 or pairs.block.nunique()!=30:
        raise RuntimeError("FROZEN_187_PAIR_SOURCE_CHANGED")
    ids=pd.concat([pairs.bright_site_id,pairs.comparison_site_id])
    if ids.duplicated().any():
        raise RuntimeError("NONREPLACEMENT_MATCH_VIOLATION")
    sites=wetland_vegetation_sites(pd.read_csv(NDVI))
    habitat=pd.read_csv(HAB).merge(pd.read_csv(COAST),on="unit_id",validate="1:1")
    habitat=habitat[(habitat.habitat.astype(str).isin(["411","412"]))&
                    habitat.coastal_classes_in_2km_square.isna()]
    sites=sites.merge(habitat[["unit_id"]],on="unit_id",validate="1:1")
    both,coverage=fixed_vegetation_overlap(pairs,sites)
    result=dict(source_sha256=sha,original_vegetation_git_blob=PIN_NDVI,
        original_corine2000_habitat_git_blob=PIN_HAB,coverage=coverage,
        fit=nested_on_same_pairs(both),
        indicator="QA0 full-three-season MOD13Q1 NDVI, 2000-2023, >=12 valid years, span >=18y, >=2 early and late years",
        scope="2000 CORINE inland marsh/peatland positions only (17 after coastal exclusion)",
        wetland_extent_measured=False,
        original_187_pair_baseline_water_balance_status="SOURCE_SUPPORT_STOP",
        out_of_sample_replication=False,
        causal_mediation=False)
    txt=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)
    print("ALAN_FIXED187_WETLAND_VEGETATION_COVERAGE\n"+txt,flush=True)
    if args.out:
        p=Path(args.out);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(txt+"\n")

if __name__=="__main__":
    main()
