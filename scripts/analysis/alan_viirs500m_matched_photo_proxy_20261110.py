#!/usr/bin/env python3
"""True 2024 500m VIIRS dark configuration source-pair support and GloBI photo bias.

Only original independent 2025 photographed plant site-years and original
GSHHG-masked 2024 satellite pixels; never biological moth abundance or
effective nocturnal pollination. Predefined same-island differences in dark
fraction <=.10, log1p(2024 VIIRS mean) <=.20, connectivity >=.10.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.alan_dark_refugia_photosource_paircheck_20261110 import (
    compare_pairs,site_cluster_inference)

COG_MD5="0b976892d40c2d0d99177c727289eaf2"
FINE_SHA="68ed09378a8bfbcb36572ff0b0578334e4f1795efe22209f437960dbb01d5896"
SOURCE_2025="moth_positive_fraction_OF_RECORDED_WCVP_SPECIES"
SCHEMA="alan_real_VIIRS500m2024_island2025_darkpatch_photo_GloBI_proxy_common_support_v1"

def require(v,m):
    if not bool(v):raise ValueError("ALAN_500M_TRUE_MATCHED_SOURCE: "+str(m))

def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for z in iter(lambda:f.read(4<<20),b""):h.update(z)
    return h.hexdigest()

def prepare(fine_csv,fine_receipt,flora_csv,flora_receipt):
    fr=json.loads(Path(fine_receipt).read_text())
    sr=json.loads(Path(flora_receipt).read_text())
    require(fr["schema"]=="ALAN_original_2025_photo_sites_VIIRS2024_500m_vs_PCNL2024_1km_dark_refugia_v1"
            and fr["output_CSV_SHA256"]==FINE_SHA and sha256(fine_csv)==FINE_SHA
            and fr["source_raster_MD5"]==COG_MD5,
            "Original genuine fine radiance source SHA/MD5 has changed")
    require(sr["source_site_year_rows"]==85318 and sr["raster_years_verified"]==33
            and sha256(flora_csv)==sr["output_sha256"],
            "Authentic original all-93m photographed-plant frame not verified")
    x=pd.read_csv(fine_csv,dtype={"island":str})
    y=pd.read_csv(flora_csv,dtype={"island":str})
    keys=["island","latbin","lonbin"]
    require(x["observation_year"].eq(2025).all() and x.satellite_year.eq(2024).all()
            and not x.duplicated(keys+["VIIRS_dark_threshold_nW_cm2_sr"]).any()
            and {0.,1.,3.}==set(x.VIIRS_dark_threshold_nW_cm2_sr),
            "Fine photo-satellite source has invalid year/site/threshold")
    require({SOURCE_2025,"all_accepted_observation_uuids",
             "all_accepted_recorded_species","year"}.issubset(y.columns),
            "Source flora numerator/effort required for negative controls")
    y=y.loc[y.year.eq(2025)].copy()
    require(not y.duplicated(keys).any(),"Original photo site-year duplicated")
    x=x.loc[x.status.eq("observed_satellite_geometry_only")].copy()
    z=x.merge(y[keys+[SOURCE_2025,"all_accepted_observation_uuids",
                     "all_accepted_recorded_species"]],how="left",on=keys,
              validate="many_to_one",indicator=True)
    require(len(z)>0 and z._merge.eq("both").all() and
            z.loc[z.VIIRS_dark_threshold_nW_cm2_sr.eq(0.),"island"].nunique()==957 and
            z.VIIRS_dark_threshold_nW_cm2_sr.eq(0.).sum()==3513,
            "Genuine original 2025 photo-site VIIRS500m matching differs")
    require(z[SOURCE_2025].between(0,1).all(),"No-source moth visit claim cannot imply moth absence")
    z["GloBI_positive_source_percent"]=100*z[SOURCE_2025]
    z["log_observation_UUIDs"]=np.log1p(z.all_accepted_observation_uuids)
    z["log_recorded_accepted_species"]=np.log1p(z.all_accepted_recorded_species)
    z["local_mean_PCNL_native"]=z.local_mean_viirs_nW_cm2_sr # solely reuse generic matched-pair code
    z["threshold_PCNL_native"]=z.VIIRS_dark_threshold_nW_cm2_sr
    return z,fr

def run(fine_csv,fine_receipt,flora_csv,flora_receipt,out):
    d,fr=prepare(fine_csv,fine_receipt,flora_csv,flora_receipt)
    results={}
    for t in (0.,1.,3.):
        pairs=compare_pairs(d,t)
        matched={}
        for name,p in [("all_matched",pairs),
                       ("nonoverlap_ge_10km",pairs.loc[pairs.distance_km.ge(10.)].copy()
                        if len(pairs) else pairs)]:
            matched[name]={
                "source_unique_independent_islands":int(p.island.nunique()) if len(p) else 0,
                "source_real_matched_site_pairs":len(p),
                "GloBI_documented_moth_visitor_plant_photo_percent":
                    site_cluster_inference(p,"diff_photo_GloBI_positive_pp"),
                "photo_observation_attempts_negative_control":
                    site_cluster_inference(p,"diff_log_photo_UUIDs"),
                "recorded_plant_species_richness_negative_control":
                    site_cluster_inference(p,"diff_log_recorded_species"),
            }
        results[str(t)]=matched
    out=Path(out)
    require(not out.exists(),"Do not overwrite true original source result")
    out.parent.mkdir(parents=True,exist_ok=True)
    result={
        "schema":SCHEMA,
        "source_VIIRS2024_original_500m_MD5":COG_MD5,
        "source_original_2025_site_exposure_CSV_SHA256":FINE_SHA,
        "actual_original_2025_photo_sites_2024_500m_valid":3513,
        "actual_500m_source_islands":957,
        "observed_source_sites_not_random_native_flora":True,
        "conditions":"within same GSHHG island and photo year 2025; dark fraction Δ<=.10, log1p(mean VIIRS radiance) Δ<=.20, conditional 4-neighbor largest dark fraction Δ>=.10; optionally two site centers ≥10km apart",
        "new_VIIRS500m_dark_connectivity_positivity_and_photo_source_proxy":results,
        "PCNL_reference":"2024 ~1km product previously gave 544 matched pairs/235 islands; ≥10km site separation 216 pairs/116 islands; separate products, not direct resolution causal comparison",
        "source_2024_VIIRS_processed_background_mask_used":True,
        "source_moth_GloBI_positive_plant_record_is_NOT_measured_insects":True,
        "not_biological_pollination_or_effective_pollen_deposition":True,
        "no_effect_of_light_on_wild_flora_estimated":True,
        "original_10m_or_40m_scene_not_measured":True,
        "CI_nominal_unadjusted_postexploratory":True,
    }
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_REAL_VIIRS_500M_DARK_PATCH_MATCHED_PHOTOPROXY",
          json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ["fine","fine_receipt","flora","flora_receipt","out"]:
        p.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    run(a.fine,a.fine_receipt,a.flora,a.flora_receipt,a.out)
