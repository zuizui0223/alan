#!/usr/bin/env python3
"""Reconcile *all* source-signed tier07 27,507 true original RGB ROI9 outcomes.

Read the independent first1024 original-photo pilot and 52 remaining warm
512-photo worker original nine-component photographic ROI receipts. Verify the
EXACT 27,507 new observation UUIDs, original photo IDs, GSHHG landmass/source
year, and source palette. Distinguish source photo opportunities from
actually successful original photo nine-palette fractions and ≥3 distinct
successfully photographed accepted taxa. No raw image bytes, no imputation,
no historic-main regression refill and no causal ALAN inference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.acquisition import alan_tier07_remaining26483_original_rgb_roi9 as acq

SOURCE_SHA=acq.SCHEDULE_SHA
FCP_BLOB=acq.FCP_SOURCE_BLOB
SCHEMA="ALAN_verified_27507_independent_tier07_original_RGB_ROI9_v1"

def require(ok,message):
    if not bool(ok):raise ValueError("ALAN_27507_REAL_ROI9_VERIFY: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def source_data(path,receipt):
    return acq.load_signed_original_source(path,receipt)[0]

def read_true_roi(file,success,expected):
    d=pd.read_csv(file,compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==expected and d.photo_id.is_unique and d.observation_id.is_unique,
            "Original photo or observation UUID source count/identity changed")
    require("strict9" in d and set(acq.COMPONENTS)<=set(d.columns),
            "Original FCP nine-colour fractions absent")
    d["measured_roi9"]=d.strict9.str.lower().eq("true")
    require(int(d.measured_roi9.sum())==success,
            "Original-photo nine-palette colour success count mismatch")
    vals=d.loc[d.measured_roi9,list(acq.COMPONENTS)].apply(
         pd.to_numeric,errors="coerce")
    require(vals.notna().all().all() and
            vals.ge(0).all().all() and
            np.allclose(vals.sum(axis=1),1,rtol=0,atol=1e-5),
            "Valid source photo numeric ROI9 component masses changed")
    unknown=d.loc[~d.measured_roi9,list(acq.COMPONENTS)]
    require(unknown.replace("",np.nan).isna().all().all(),
            "A technical original-photo ROI unknown was imputed")
    return d.loc[:,["photo_id","observation_id","measured_roi9"]]

def read_pilot(path):
    root=Path(path)
    file=root/"finished"/"ALAN_tier07_first1024_original_floral_roi9.csv.gz"
    meta=json.loads((root/"finished"/"ALAN_tier07_first1024_original_roi9_receipt.json").read_text())
    n=1024
    require(meta.get("schema")=="ALAN_tier07_1024_real_original_FCP_ROI9_pilot_v1"
           and meta["new_independent_observation_UUIDs"]==n
           and meta["original_image_attempts"]==n
           and meta["source_distinct_GSHHG_islands"]==284
           and meta["source_5photo_3accepted_potential_site_years"]==76
           and meta["per_128_photo_unit_receipts_verified"]==8
           and meta["original_source_schedule_sha256"]==SOURCE_SHA
           and meta["original_frozen_FCP_ROI9_method"]["measurement_git_blob"]==FCP_BLOB
           and meta["real_original_ROI9_photo_outcomes_sha256"]==sha(file),
           "Pilot 1,024 true-photo FCP source numeric receipt not independently valid")
    return read_true_roi(file,meta["original_9palette_numeric_ROI_success"],n),{
        "source_part":"pilot1024","true_photo_attempts":n,
        "measured_original_numeric_ROI9":meta["original_9palette_numeric_ROI_success"],
        "photo_archive_SHA256":sha(file)
    }

def read_remaining(root,task):
    root=Path(root)
    file=root/"finished"/f"alan_tier07_original_photo_ROI9_task{task:02d}.csv.gz"
    meta=json.loads((root/"finished"/"ALAN_tier07_genuine_new_photo_ROI9_task_receipt.json").read_text())
    n=acq.photo_count_for_task(task)
    require(meta.get("schema")=="ALAN_tier07_26483_original_photo_FCP_ROI9_task_v1"
            and meta["task_number"]==task
            and meta["original_photo_image_attempts"]==n
            and meta["independent_new_photo_observation_UUIDs"]==n
            and meta["previous_1024_original_photo_UUIDs_not_reselected"] is True
            and meta["source_schedule_SHA256"]==SOURCE_SHA
            and meta["original_FCP_model_sha_blob"]==FCP_BLOB
            and meta["original_real_photo_colour_file_SHA256"]==sha(file)
            and meta["technical_unknown_not_assigned_an_unnatural_colour"] is True
            and meta["original_RGB_image_bytes_not_persisted"] is True,
            "One 512-photo source worker numeric ROI9 output failed independent SHA check")
    return read_true_roi(file,meta["actual_numerical_9palette_flowerROI_success"],n),{
        "source_task":task,"true_photo_attempts":n,
        "measured_original_numeric_ROI9":meta["actual_numerical_9palette_flowerROI_success"],
        "image_acquisition_failures":meta["original_photo_acquisition_failures"],
        "photo_archive_SHA256":sha(file)
    }

def source_and_actual_sites(frame):
    require(frame.photo_id.is_unique and frame.observation_id.is_unique,
            "The same original observation/photo was counted twice")
    by=frame.groupby("target_id",observed=True).agg(
        original_source_photo_attempts=("observation_id","nunique"),
        source_accepted_species=("wcvp_accepted_plant_name_id","nunique"),
        numeric_original_photoROI9=("measured_roi9","sum"))
    success=frame.loc[frame.measured_roi9]
    measured_species=success.groupby("target_id",observed=True).wcvp_accepted_plant_name_id.nunique()
    by["numeric_measured_accepted_species"]=measured_species.reindex(
        by.index,fill_value=0).astype(int)
    by["opportunity_ge5photo_ge3taxa"]=(
        by.original_source_photo_attempts.ge(5)&by.source_accepted_species.ge(3))
    by["actual_strict_ge5photo_ge3trueROI_ge3taxa"]=(
        by.opportunity_ge5photo_ge3taxa&
        by.numeric_original_photoROI9.ge(3)&
        by.numeric_measured_accepted_species.ge(3))
    return by

def run(schedule,receipt,pilot,remaining,out):
    out=Path(out)
    require(not out.exists(),"Never replace prior signed 27k original photo ROI9 result")
    source=source_data(schedule,receipt)
    require(len(remaining)==acq.TASK_COUNT and
            set(remaining)==set(range(acq.TASK_COUNT)),
            "Partial original source 52 warm photo worker set")
    frames=[];proof=[]
    p,r=read_pilot(pilot)
    frames.append(p);proof.append(r)
    for task,root in sorted(remaining.items()):
        d,meta=read_remaining(root,task)
        frames.append(d);proof.append(meta)
    actual=pd.concat(frames,ignore_index=True)
    require(len(actual)==27507 and actual.photo_id.is_unique and
            actual.observation_id.is_unique,
            "27,507 newly measured source photo/observation IDs not unique")
    joined=source.merge(actual,on=["photo_id","observation_id"],
                        how="inner",validate="one_to_one")
    require(len(joined)==27507 and
            set(joined.photo_id)==set(source.photo_id) and
            set(joined.observation_id)==set(source.observation_id),
            "Measured RGB original photo IDs do not exactly match signed precolour source")
    by_site=source_and_actual_sites(joined)
    require(len(by_site)==9493 and joined.land_id.nunique()==563,
            "Original source GSHHG islands/photo cell-years changed")
    counts=joined.groupby("work_group",observed=True).agg(
        source_attempts=("observation_id","nunique"),
        measured_numeric_roi9=("measured_roi9","sum"),
        distinct_GSHHG_islands=("land_id","nunique"),
        genuine_source_site_years=("target_id","nunique")
    ).reset_index()
    years=joined.groupby("year",observed=True).agg(
        original_photo_attempts=("observation_id","nunique"),
        measured_numeric_roi9=("measured_roi9","sum"),
        source_photo_cell_years=("target_id","nunique"),
        source_GSHHG_islands=("land_id","nunique")
    ).reset_index()
    index=by_site.index.astype(str).str.split("|",regex=False)
    by_site["year"]=[x[0] for x in index]
    src_by_year=by_site.groupby("year",observed=True).agg(
        source_opportunity_sites=("opportunity_ge5photo_ge3taxa","sum"),
        actually_strict_measured_sites=("actual_strict_ge5photo_ge3trueROI_ge3taxa","sum"),
        all_source_photo_cells=("original_source_photo_attempts","size")
    ).reset_index()
    pilot_actual=joined.loc[joined.work_group.eq("0")].copy()
    pilot_site=source_and_actual_sites(pilot_actual)
    require(len(pilot_actual)==1024 and
            int(pilot_site.opportunity_ge5photo_ge3taxa.sum())==76,
            "Original tier07 source pilot precolour 76 sites changed")
    outcome={
       "schema":SCHEMA,
       "signed_original_source_schedule_SHA256":SOURCE_SHA,
       "source_original_observation_UUIDs_and_photo_IDs_verified":27507,
       "original_source_distinct_GSHHG_islands":563,
       "original_source_photo_site_years":9493,
       "original_numeric_9colour_FCP_flowerROI9_successes":int(actual.measured_roi9.sum()),
       "technical_unclassifiable_original_photos":int((~actual.measured_roi9).sum()),
       "source_opportunity_5photo_3accepted_taxon_site_years":int(
          by_site.opportunity_ge5photo_ge3taxa.sum()),
       "actual_strict_5photo_3numericROI9_3measured_accepted_taxon_site_years":int(
          by_site.actual_strict_ge5photo_ge3trueROI_ge3taxa.sum()),
       "pilot_first1024_9colour_successes":int(pilot_actual.measured_roi9.sum()),
       "pilot_first1024_real_strict_site_years":int(
          pilot_site.actual_strict_ge5photo_ge3trueROI_ge3taxa.sum()),
       "first_pilot_and_52_signed_warm_photo_workers_verified":len(proof),
       "all_photo_UUIDs_unique_and_match_source_before_pixel_inference":True,
       "no_source_original_photo_pixels_downloaded_or_reinferred_during_audit":True,
       "unclassified_original_image_ROI_never_imputed":True,
       "observational_source_island_targeted_not_unbiased_global_flora":True,
       "original_historical_main_island_x_ALAN_coefficients_unmodified":True,
       "frozen_source_ROI9_evaluator_git_blob":FCP_BLOB,
       "independent_source_photo_worker_SHA_receipts":proof
    }
    require(outcome["original_numeric_9colour_FCP_flowerROI9_successes"]+
            outcome["technical_unclassifiable_original_photos"]==27507,
            "Measured+unknown genuine new source originals do not sum")
    out.mkdir(parents=True)
    for name,data in [
      ("original_source_photo_group_measured_counts.csv",counts),
      ("original_source_photo_year_measured_counts.csv",years),
      ("original_source_strict_ROI9_site_year_counts.csv",src_by_year)]:
        data.to_csv(out/name,index=False)
        outcome[name.replace(".csv","_SHA256")]=sha(out/name)
    (out/"ALAN_original_tier07_27507_real_photoROI9_receipt.json").write_text(
        json.dumps(outcome,sort_keys=True,indent=2)+"\n")
    print("ALAN_TIER07_27507_ACTUAL_ORIGINAL_RGB_ROI9_SOURCE_RESULT",
          json.dumps({k:outcome[k] for k in (
            "source_original_observation_UUIDs_and_photo_IDs_verified",
            "original_numeric_9colour_FCP_flowerROI9_successes",
            "technical_unclassifiable_original_photos",
            "source_opportunity_5photo_3accepted_taxon_site_years",
            "actual_strict_5photo_3numericROI9_3measured_accepted_taxon_site_years",
            "pilot_first1024_9colour_successes",
            "pilot_first1024_real_strict_site_years")
           },sort_keys=True),flush=True)
    return outcome

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schedule",type=Path,required=True)
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--pilot",type=Path,required=True)
    for task in range(acq.TASK_COUNT):
        p.add_argument(f"--task{task}",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.schedule,args.source_receipt,args.pilot,
        {t:getattr(args,f"task{t}") for t in range(acq.TASK_COUNT)},
        args.out)
