#!/usr/bin/env python3
"""Source-verified *signed* growing-season wetland-position greenness changes.

Compare original MOD13Q1 QA0 annual NDVI with independently mapped 3-epoch
CORINE land-cover change. NDVI measures greenness, NOT wetland area or plant
species abundance. Botanical survey-selection and single brightening site
preclude ALAN effect attribution. This is NOT the 187 matched-pair analysis.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, theilslopes
from alan_wet_vegetation_identical_pairs_20261120 import (
    git_blob, wetland_vegetation_sites
)

ROOT=Path(__file__).resolve().parents[2]
DIR=ROOT/"results/all_species_forest/wet_annual"
CLC=ROOT/"results/all_species_forest/wet_falsification/corine_official_cha_17wetland_20261120.json"
VEG=DIR/"vegetation.csv.gz"
SITE=DIR/"site_trends.csv"
PIN={
  "veg":"26b18e38149d78d0914e5a89c1f863d766c93e4b",
  "site":"6fd31c8bd355104179ddec7394d9cb2744c4b312",
  "corine":"866fa4d7750eb8fae99c56ede68a05bfcbfce6ab"
}
SEED=20261120

def seasonal_annual_ndvi(raw):
    """Three full Apr-May, Jun-Jul, Aug-Sep medians; failed QA is missing."""
    req={"unit_id","y","qa","day","ndvi"}
    if req-set(raw): raise ValueError("SOURCE_SCHEMA_STOP")
    x=raw.dropna(subset=list(req)).copy()
    x=x[(x.qa==0)&x.y.between(2000,2023)&x.ndvi.between(-.2,1.)]
    x["date"]=pd.to_datetime(x.y.astype(int).astype(str)+"-01-01")+pd.to_timedelta(x.day-1,unit="D")
    x=x[x.date.dt.month.between(4,9)].copy()
    x["season"]=((x.date.dt.month-4)//2).astype(int)
    a=x.groupby(["unit_id","y","season"]).ndvi.median().unstack("season")
    for i in range(3):
        if i not in a:a[i]=np.nan
    a=a[[0,1,2]].mean(axis=1).where(a[[0,1,2]].notna().all(axis=1)).rename("ndvi")
    return a.dropna().reset_index().sort_values(["unit_id","y"])

def source_eligibility(annual, valid_sites):
    """Fixed 17 original site IDs plus independent early/late QA time coverage."""
    rows=[]
    for site,g in annual.groupby("unit_id"):
        if int(site) not in set(valid_sites):
            continue
        g=g.sort_values("y")
        early=g[g.y.between(2000,2004)].ndvi
        late=g[g.y.between(2019,2023)].ndvi
        if len(g)<12 or g.y.max()-g.y.min()<18 or len(early)<2 or len(late)<2:
            continue
        slope=10*float(np.polyfit(g.y.to_numpy(float),g.ndvi.to_numpy(float),1)[0])
        ts=10*float(theilslopes(g.ndvi.to_numpy(float),g.y.to_numpy(float))[0])
        rows.append(dict(unit_id=int(site),years=int(len(g)),
            ndvi_initial=float(early.median()),ndvi_final=float(late.median()),
            early_to_late_ndvi=float(late.median()-early.median()),
            ols_ndvi_per_decade=slope, theilsen_ndvi_per_decade=ts))
    return pd.DataFrame(rows)

def join_original_sources(ndvi,site,clc):
    if not (len(site)==17 and site.unit_id.is_unique and len(ndvi)==17 and ndvi.unit_id.is_unique):
        raise RuntimeError("17_ORIGINAL_WETLAND_SITE_SUPPORT_STOP")
    rec=clc["site_rows"]
    if len(rec)!=51 or pd.DataFrame(rec).duplicated(["unit_id","epoch"]).any():
        raise RuntimeError("51_ORIGINAL_CLC_EPOCH_STOP")
    wet=pd.DataFrame(rec).groupby("unit_id")[[
        "changed_area_ha","inland_wetland_loss_ha","inland_wetland_gain_ha"
        ]].sum().reset_index()
    if len(wet)!=17:raise RuntimeError("17_CLC_SITES_NOT_COMPLETE")
    z=site[["unit_id","block","light","difference","candidate_rate","reference_rate"]].merge(
        ndvi,on="unit_id",validate="1:1").merge(wet,on="unit_id",validate="1:1")
    if len(z)!=17:raise RuntimeError("INCOMPLETE_SOURCE_JOIN_STOP")
    if not np.isfinite(z[["ols_ndvi_per_decade","theilsen_ndvi_per_decade",
                         "early_to_late_ndvi","difference"]].to_numpy(float)).all():
        raise RuntimeError("INVALID_SOURCE_VALUES")
    return z.sort_values("unit_id").reset_index(drop=True)

def spatial_block_bootstrap(z):
    groups=[g for _,g in z.groupby("block")]
    if len(groups)<5:
        return {"status":"BLOCK_SUPPORT_STOP","blocks":len(groups)}
    rho=float(spearmanr(z.ols_ndvi_per_decade,z.difference).statistic)
    rng=np.random.default_rng(SEED)
    means=[]; rhos=[]
    for _ in range(2500):
        b=pd.concat([groups[i] for i in rng.integers(0,len(groups),len(groups))],ignore_index=True)
        means.append(float(b.ols_ndvi_per_decade.mean()))
        if b.ols_ndvi_per_decade.nunique()>2 and b.difference.nunique()>2:
            rhos.append(float(spearmanr(b.ols_ndvi_per_decade,b.difference).statistic))
    return {"status":"DESCRIPTIVE_ONLY",
        "independent_blocks":int(z.block.nunique()),
        "site_level_spearman_ndvi_vs_candidate_minus_reference":rho,
        "spatial_block_bootstrap_mean_ndvi95":np.quantile(means,[.025,.975]).tolist(),
        "spatial_block_bootstrap_spearman95":np.quantile(rhos,[.025,.975]).tolist()
            if len(rhos)>=500 else None,
        "valid_spatial_block_draws_for_spearman":len(rhos)}

def summarize(z):
    n=len(z)
    diff=z.ols_ndvi_per_decade
    signed= z.theilsen_ndvi_per_decade
    changes=z.early_to_late_ndvi
    result={"source_n_sites":n,"independent_blocks":int(z.block.nunique()),
      "sites_ndvi_ols_negative":int((diff<0).sum()),
      "sites_ndvi_ols_positive":int((diff>0).sum()),
      "sites_ndvi_theilsen_negative":int((signed<0).sum()),
      "sites_consistent_negative_ols_and_theilsen":int(((diff<0)&(signed<0)).sum()),
      "sites_early_late_ndvi_decreased":int((changes<0).sum()),
      "mean_ols_ndvi_change_per_decade":float(diff.mean()),
      "median_ols_ndvi_change_per_decade":float(diff.median()),
      "mean_early_late_ndvi_delta":float(changes.mean()),
      "median_early_late_ndvi_delta":float(changes.median()),
      "sites_with_corine_wetland_class_loss":int((z.inland_wetland_loss_ha>0).sum()),
      "sites_with_corine_wetland_class_gain":int((z.inland_wetland_gain_ha>0).sum()),
      "sites_with_any_mapped_landcover_change":int((z.changed_area_ha>0).sum()),
      "sites_negative_greenness_without_mapped_wetland_loss":
          int(((diff<0)&(z.inland_wetland_loss_ha==0)).sum()),
      "sites_with_20y_light_brightening":int((z.light>0).sum()),
      "spatial_bootstrap":spatial_block_bootstrap(z),
      "site_level_data":[{k:(int(v) if k in ("unit_id","block","years") else float(v))
          for k,v in row.items()} for row in z.to_dict("records")],
      "caveats":["17 previously selected wetland sites, not 187 original matched pairs",
                 "only one site with increasing long-term nightlights",
                 "NDVI greenness may reflect cutting, phenology or other plants not wetland extent",
                 "MODIS 250-m point and CORINE 1-km circle are different measurement footprints",
                 "CORINE changes 2000–2018; NDVI covers 2000–2023",
                 "vegetation and botanical records not independent of selection",
                 "block bootstraps with eight blocks only are descriptive, not confirmatory"]}
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="")
    args=p.parse_args()
    for key,path in (("veg",VEG),("site",SITE),("corine",CLC)):
        if git_blob(path)!=PIN[key]:
            raise RuntimeError("ORIGINAL_SOURCE_BLOB_PIN_STOP "+key)
    ndvi_source=pd.read_csv(VEG)
    site=pd.read_csv(SITE)
    clc=json.loads(CLC.read_text())
    good=site.unit_id.to_list()
    annual=seasonal_annual_ndvi(ndvi_source)
    strict=source_eligibility(annual,good)
    original=wetland_vegetation_sites(ndvi_source)
    chk=strict.merge(original[["unit_id","ndvi_trend"]],on="unit_id",validate="1:1")
    if len(chk)!=17 or not np.allclose(chk.ols_ndvi_per_decade,chk.ndvi_trend,rtol=1e-9,atol=1e-9):
        raise RuntimeError("NDVI_INDEPENDENT_RECOMPUTATION_MISMATCH")
    joined=join_original_sources(strict,site,clc)
    result={"source_git_blobs":PIN,
      "source_context":"original 17 CORINE-2000 inland marsh/bog sites",
      "site_coverage_policy":"QA0 three-season annual NDVI, >=12 years, >=18 span, early/late >=2 years",
      "matched_187_pair_vegetation_control_estimated":False,
      "associated_with_ALAN_causally":False,
      "results":summarize(joined)}
    txt=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)
    print("ALAN_WETLAND_SIGNED_NDVI_17_SOURCE_RESULT\n"+txt,flush=True)
    if args.out:
        out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(txt+"\n")
if __name__=="__main__":
    main()
