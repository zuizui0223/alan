#!/usr/bin/env python3
"""Real original 2025 island dark-refuge exposure vs *GloBI documented plant* proxy.

Exploratory SOURCE / ASCERTAINMENT diagnostic, NOT moth abundance or pollination.
Pairs matched on original PCNL mean & dark-land fraction, not any outcome.
Independent uncertainty unit is island. Also test 10-km separation of 5-km
overlapping satellite neighbourhoods and recorded-effort negative controls.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.stats import t as student_t, ttest_1samp

EXPOSURE_SHA="32327abb37ee85f3d1d9e7446d9d742fb5adf1457dab60fae1aae791c29f7837"
SOURCE_FIELD="moth_positive_fraction_OF_RECORDED_WCVP_SPECIES"
EXPOSURE_THRESHOLD=(0.,3.)
GEOD=Geod(ellps="WGS84")
SCHEMA="ALAN_2025_island_dark_mosaic_photo_GloBI_positive_source_proxy_matched_v1"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_DARK_MOTH_SOURCE_PROXY: "+str(msg))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for q in iter(lambda:f.read(4<<20),b""):
            h.update(q)
    return h.hexdigest()

def load_real(exposure, exposure_receipt, source, source_receipt):
    e=json.loads(Path(exposure_receipt).read_text())
    s=json.loads(Path(source_receipt).read_text())
    require(e["schema"]=="alan_dark_refugia_original_2025_sites_prior_2024_PCNL_GSHHG_v1" and
            e["source_output_csv_sha256"]==EXPOSURE_SHA and sha(exposure)==EXPOSURE_SHA,
            "Original PCNL land-only connected exposure changed")
    require(s["source_site_year_rows"]==85318 and s["raster_years_verified"]==33 and
            sha(source)==s["output_sha256"],
            "Original authentic photographed plant source drift")
    x=pd.read_csv(exposure,dtype={"island":str})
    y=pd.read_csv(source,dtype={"island":str},low_memory=False)
    keys=["island","latbin","lonbin"]
    require({"threshold_PCNL_native","dark_land_fraction",
             "largest_dark_component_fraction_of_dark","local_mean_PCNL_native",
             "status","lat","lon"}.issubset(x.columns),"Wrong recorded 5km raster schema")
    require({SOURCE_FIELD,"all_accepted_recorded_species",
             "all_accepted_observation_uuids","year"}.issubset(y.columns),
            "Original authentic 93m photo source lacks GloBI documentation support")
    y=y.loc[y.year.eq(2025)].copy()
    require(not y.duplicated(keys).any(),"True source cell repeated after year filter")
    require(not x.duplicated(keys+["threshold_PCNL_native"]).any(),
            "Source site threshold pseudoreplicated")
    x=x.loc[x.status.eq("observed_satellite_geometry_only") &
            x.threshold_PCNL_native.isin(EXPOSURE_THRESHOLD)].copy()
    z=x.merge(y[keys+[SOURCE_FIELD,"all_accepted_recorded_species",
                       "all_accepted_observation_uuids"]],
              how="left",on=keys,validate="many_to_one",indicator=True)
    require(z._merge.eq("both").all() and len(z)==2*3543,
            "Original satellite/photo-site identity source crossmatch changed")
    z["GloBI_positive_source_percent"]=100*z[SOURCE_FIELD]
    z["log_observation_UUIDs"]=np.log1p(z.all_accepted_observation_uuids)
    z["log_recorded_accepted_species"]=np.log1p(z.all_accepted_recorded_species)
    require(z["GloBI_positive_source_percent"].between(0,100).all() and
            z.log_observation_UUIDs.notna().all(),"Invalid real photorecord sources")
    return z

def compare_pairs(x, threshold):
    d=x.loc[x.threshold_PCNL_native.eq(threshold)].copy()
    d=d.loc[d.n_dark_raster_pixels.gt(0) &
            d.largest_dark_component_fraction_of_dark.notna()].copy()
    rows=[]
    for island,g in d.groupby("island",sort=True):
        r=list(g.itertuples(index=False))
        for i in range(len(r)):
            for j in range(i+1,len(r)):
                a,b=r[i],r[j]
                dc=float(a.largest_dark_component_fraction_of_dark-
                         b.largest_dark_component_fraction_of_dark)
                if abs(dc)<.10:continue
                df=abs(float(a.dark_land_fraction-b.dark_land_fraction))
                if df>.10:continue
                dl=abs(np.log1p(float(a.local_mean_PCNL_native))-
                       np.log1p(float(b.local_mean_PCNL_native)))
                if dl>.20:continue
                _,_,m=GEOD.inv(float(a.lon),float(a.lat),float(b.lon),float(b.lat))
                distance_km=float(m)/1000
                higher,lower=(a,b) if dc>0 else (b,a)
                rows.append({
                    "island":island,"distance_km":distance_km,
                    "abs_dark_cover_delta":df,"abs_logmean_light_delta":dl,
                    "abs_patch_connectivity_delta":abs(dc),
                    "diff_photo_GloBI_positive_pp":float(
                        higher.GloBI_positive_source_percent-lower.GloBI_positive_source_percent),
                    "diff_log_photo_UUIDs":float(
                        higher.log_observation_UUIDs-lower.log_observation_UUIDs),
                    "diff_log_recorded_species":float(
                        higher.log_recorded_accepted_species-lower.log_recorded_accepted_species)})
    return pd.DataFrame(rows)

def site_cluster_inference(pairs,column,atleast=8):
    if pairs.empty:return {"status":"no_matched_original_pairs"}
    within=pairs.groupby("island",sort=True)[column].mean()
    n=len(within)
    if n<atleast:return {"status":"too_few_true_islands","islands":n,
                         "pairs":len(pairs)}
    mn=float(within.mean())
    st=float(within.std(ddof=1))
    se=st/np.sqrt(n)
    cv=float(student_t.ppf(.975,n-1))
    p=float(ttest_1samp(within,0).pvalue) if st>0 else None
    return {"status":"observed_source_diagnostic_only",
            "independent_original_islands":n,
            "source_matched_photo_site_pairs":len(pairs),
            "mean_island_equal_higher_connected_minus_lower":mn,
            "nominal_95CI":[mn-cv*se,mn+cv*se],
            "nominal_unadjusted_p":p,
            "unit":"island means, not photo site pairs or repeat thresholds",
            "never_moth_abundance_or_real_pollination":True}

def summarize(frame):
    result={"schema":SCHEMA,
            "original_2025_photo_sites_with_2024_real_raster":3543,
            "documented_moth_visitor_source_is_not_evidence_of_true_effective_pollination":True,
            "land_masked_spectral_PCNL_is_not_lux_or_streetlamp_distance":True,
            "no_unphotographed_species_absence_assumed":True,
            "result_is_post_hoc_ascertainment_diagnostic_not_confirmatory":True,
            "per_threshold":{}}
    for t in EXPOSURE_THRESHOLD:
        pairs=compare_pairs(frame,t)
        out={"all_candidate_matched_pairs":len(pairs),
             "all_candidate_matched_islands":int(pairs.island.nunique()) if len(pairs) else 0,
             "pairs_at_least_10km_separation":
                 int(pairs.distance_km.ge(10).sum()) if len(pairs) else 0,
             "islands_at_least_10km_separation":int(
                pairs.loc[pairs.distance_km.ge(10),"island"].nunique()) if len(pairs) else 0}
        for label,f in [
            ("all_matched_5km_buffer_may_overlap",pairs),
            ("at_least_10km_separation_nonoverlapping_5km_buffers",
             pairs.loc[pairs.distance_km.ge(10)].copy() if len(pairs) else pairs)
        ]:
            out[label]={
                "GloBI_moth_documented_plant_percent":site_cluster_inference(f,"diff_photo_GloBI_positive_pp"),
                "negative_control_photo_observation_effort":site_cluster_inference(f,"diff_log_photo_UUIDs"),
                "negative_control_recorded_species_richness":site_cluster_inference(f,"diff_log_recorded_species")}
        result["per_threshold"][str(t)]=out
    return result

def run(exposure,exposure_receipt,source,source_receipt,out):
    frame=load_real(exposure,exposure_receipt,source,source_receipt)
    result=summarize(frame)
    out=Path(out)
    require(not out.exists(),"No overwrite of source-derived results")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_DARK_PHOTO_SOURCE_PROXY_PAIRED_DIAGNOSTIC",
          json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--exposure",required=True,type=Path)
    p.add_argument("--exposure-receipt",required=True,type=Path)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--source-receipt",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.exposure,a.exposure_receipt,a.source,a.source_receipt,a.out)
