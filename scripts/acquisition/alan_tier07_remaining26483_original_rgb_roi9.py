#!/usr/bin/env python3
"""ALAN scalable original RGB/ROI9 inference, source tier07 work groups 1..26.

Source: exactly 27,507 unique observation UUID/photo ID records frozen and
colour-blind from signed 1,193,829-photo tier07 of the original 93.2m index.
The first 1,024 (source work_group 0) are a SEPARATE existing pilot and
MUST NOT be acquired here. Remaining 26,483 genuine new original-photo UUIDs
are divided into 52 tasks of at most 512 actual source photos:
51 tasks of 512 plus task 51 of 371. One warm FCP YOLO11/EfficientSAM
initialization per task and up to four authenticated 128-photo work units.

No photo is treated as a confirmed flowering native plant, photo not classifiable
remains technical unknown, licensed RGB bytes never enter final artifacts.
The island/site-focused source is not an unbiased global flower-colour sample.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scripts.acquisition import mois_full_index_real_original_photo_measure as worker
from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import (
    valid_nine, COMPONENTS, FRACTIONS,
)

SCHEDULE_SHA="7198455819f1ac260a6cbeacd9f4536c30d766595daa64dc5205c75e9d5a1c62"
SCHEDULE_SCHEMA="ALAN_tier07_27507_new_original_island_1024photo_source_schedule_v1"
FCP_SOURCE_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
FCP_SOURCE_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
SOURCE_ARM="alan_93m_tier07_new_island_observation_ROI_source_v1"
TASK_COUNT=52
TASK_FULL=512
LAST_COUNT=371
FIRST_COMPLETED=1024
REMAINING=26483
SOURCE_TOTAL=27507

def require(test,message):
    if not bool(test):
        raise ValueError("ALAN_TIER07_BULK_SOURCE: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for part in iter(lambda:stream.read(4<<20),b""):
            h.update(part)
    return h.hexdigest()

def task_source_groups(task):
    require(isinstance(task,int) and 0<=task<TASK_COUNT,
            "Task outside 52 source-disjoint jobs")
    return 1+task//2,task%2

def photo_count_for_task(task):
    return LAST_COUNT if task==TASK_COUNT-1 else TASK_FULL

def load_signed_original_source(schedule,receipt):
    require(sha(schedule)==SCHEDULE_SHA,
            "Source tier07 full27507 precolour SHA changed")
    record=json.loads(Path(receipt).read_text())
    require(record.get("schema")==SCHEDULE_SCHEMA and
            record.get("precolour_source_photo_plan_SHA256")==SCHEDULE_SHA and
            record.get("original_source_accepted_new_observation_UUIDs")==SOURCE_TOTAL and
            record.get("disjoint_genuine_source_groups")==27 and
            record.get("source_photo_selection_did_not_use_colour") is True and
            record.get("original_new_image_bytes_downloaded")==0,
            "Source not the signed original colour-blind tier07 study")
    d=pd.read_csv(schedule,compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==SOURCE_TOTAL and d.photo_id.is_unique and
            d.observation_id.is_unique and d.source_arm.eq(SOURCE_ARM).all() and
            d.source_colour_measured.astype(str).str.lower().eq("false").all() and
            d.photo_url_large.str.fullmatch(
                r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
            ).all() and d.photo_id.str.fullmatch(r"\d+").all(),
            "Original source photo IDs/UUIDs or URL formats changed")
    counts=d.groupby("work_group").size()
    require(len(counts)==27 and
            all(counts.get(str(i))==1024 for i in range(26)) and
            counts.get("26")==883 and
            int(d.work_group.eq("0").sum())==FIRST_COMPLETED,
            "Frozen original source 27 groups are not 26x1024 + 883")
    return d,record

def select_task(frame,task):
    group,half=task_source_groups(task)
    q=frame.loc[frame.work_group.eq(str(group))].copy()
    q["within_work_group"]=pd.to_numeric(
        q.within_work_group,errors="raise").astype(int)
    require(q.within_work_group.is_unique and
            q.within_work_group.min()==0 and
            q.within_work_group.max()==len(q)-1,
            "Source precolour original photo group identity/order invalid")
    q=q.loc[q.within_work_group.ge(half*TASK_FULL)&
            q.within_work_group.lt((half+1)*TASK_FULL)].copy()
    q=q.sort_values("within_work_group",kind="stable")
    require(len(q)==photo_count_for_task(task) and
            q.photo_id.is_unique and q.observation_id.is_unique and
            q.photo_id.str.fullmatch(r"\d+").all() and
            q.wcvp_accepted_plant_name_id.ne("").all() and
            q.year.isin(("2019","2021","2023","2024","2025")).all(),
            "Source partial task photo/UUID/year taxonomy drift")
    require(not q.observation_id.isin(set(frame.loc[
       frame.work_group.eq("0"),"observation_id"])).any(),
       "Previously measured original tier07 1024 photos included in new task")
    q["measurement_id"]=q.photo_id.map(
       lambda u:f"alan_tier07_t{task:02d}_source_{u}")
    q["image_filename"]=q.measurement_id+".jpg"
    q["shard"]=(np.arange(len(q))//128).astype(int)
    q["source_snapshot_date"]="2026-09-27"
    units=(len(q)+127)//128
    require(units==(3 if task==TASK_COUNT-1 else 4) and
            q.shard.nunique()==units and
            all(len(part)<=128 for _,part in q.groupby("shard")),
            "Source units do not pack into original FCP 128-photo capacity")
    return q

def prepare(schedule,receipt,task,out):
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite signed original photo source")
    master,_=load_signed_original_source(schedule,receipt)
    photo=select_task(master,task)
    (out/"shards").mkdir(parents=True)
    checks={}
    for idx,part in photo.groupby("shard",sort=True):
        dest=out/"shards"/f"shard-{idx:03d}.csv"
        part.to_csv(dest,index=False,lineterminator="\n")
        checks[str(idx)]=sha(dest)
    record={
      "schema":FCP_SOURCE_SCHEMA,"colour_measured":False,
      "source_group":f"alan_tier07_original_new_26483_task{task:02d}",
      "source_arm":SOURCE_ARM,"work_units":len(checks),
      "worker_ids":list(range(len(checks))),
      "new_photo_attempts_planned":len(photo),
      "photo_shard_sha256":checks,
      "immutable_original_27507_source_SHA256":SCHEDULE_SHA,
      "original_photo_task_number":task,
      "original_work_group":task_source_groups(task)[0],
      "original_work_group_half":task_source_groups(task)[1],
      "source_excludes_original_pilot1024":True,
      "genuine_source_observation_UUIDs_unique":True,
      "source_photo_colour_not_observed_before_selection":True,
      "island_priority_sampling_NOT_worldwide_random":True,
    }
    (out/"manifest.json").write_text(
        json.dumps(record,indent=2,sort_keys=True)+"\n")
    print("ALAN_TIER07_26483_NEW_SOURCE_TASK_FROZEN",
          json.dumps({"task":task,"actual_new_photos":len(photo),
          "GSHHG_islands":int(photo.land_id.nunique()),"units":len(checks)},sort_keys=True),flush=True)
    return record

def collect(source,measured,task,out):
    source,measured,out=map(Path,(source,measured,out))
    require(not out.exists(),"Refuse overwriting authentic original photo outcome")
    manifest=json.loads((source/"manifest.json").read_text())
    require(manifest.get("schema")==FCP_SOURCE_SCHEMA and
            manifest.get("original_photo_task_number")==task and
            manifest.get("immutable_original_27507_source_SHA256")==SCHEDULE_SHA and
            manifest.get("source_excludes_original_pilot1024") is True,
            "Source photo/UUID task contract has drifted")
    frames=[];receipts=[]
    for idx in manifest["worker_ids"]:
        original,locked=worker.input_rows(source,idx)
        folder=measured/f"shard_{idx:03d}"
        receipt=json.loads((folder/"receipt.json").read_text())
        data_file=folder/"photos.csv"
        data=pd.read_csv(data_file,dtype=str,keep_default_na=False)
        require(receipt.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
            and receipt.get("shard")==idx and
            receipt.get("attempted")==len(original) and
            receipt.get("source_sha256")==locked["photo_shard_sha256"][str(idx)] and
            receipt.get("photos_sha256")==sha(data_file) and
            len(data)==len(original) and data.photo_id.is_unique and
            data.observation_id.is_unique and
            set(data.photo_id)==set(original.photo_id) and
            set(data.observation_id)==set(original.observation_id),
            "Authentic FCP worker original photo receipts/identities mismatch")
        require(set(FRACTIONS)<=set(data.columns),
                "Fixed original FCP ROI9 did not record nine measured fractions")
        frames.append(data);receipts.append(receipt)
    all_photo=pd.concat(frames,ignore_index=True)
    methods={json.dumps(r["method"],sort_keys=True) for r in receipts}
    require(len(all_photo)==photo_count_for_task(task) and
            all_photo.photo_id.is_unique and all_photo.observation_id.is_unique
            and len(methods)==1,"Original source photo IDs duplicated or mixed ROI methods")
    method=receipts[0]["method"]
    require(method.get("measurement_git_blob")==FCP_SOURCE_BLOB and
            method.get("location_blind") is True and
            method.get("persisted_image_pixels") is False,
            "Fixed image inference method no longer source-identical")
    good,colour=valid_nine(all_photo)
    all_photo["strict9"]=good
    for i,ch in enumerate(COMPONENTS):
        all_photo[ch]=colour[:,i]
    require(np.allclose(all_photo.loc[good,list(COMPONENTS)].sum(axis=1),1,
                        atol=1.e-6,rtol=0) and
            all_photo.loc[~good,list(COMPONENTS)].isna().all().all(),
            "Technical unclassified ROI9 colour was imputed or original palette invalid")
    out.mkdir(parents=True)
    columns=["photo_id","observation_id","inat_taxon_id",
             "wcvp_accepted_plant_name_id","land_id","target_id",
             "year","observed_on","source_arm","source_priority",
             "work_group","within_work_group","morph","measurement_status",
             "image_sha256","strict9",*COMPONENTS]
    require(set(columns)<=set(all_photo.columns),
            "Actual numeric FCP output not connected to original provenance")
    result_file=out/f"alan_tier07_original_photo_ROI9_task{task:02d}.csv.gz"
    all_photo[columns].to_csv(result_file,index=False,
          compression={"method":"gzip","mtime":0},lineterminator="\n")
    stats={
       "schema":"ALAN_tier07_26483_original_photo_FCP_ROI9_task_v1",
       "task_number":task,
       "source_schedule_SHA256":SCHEDULE_SHA,
       "previous_1024_original_photo_UUIDs_not_reselected":True,
       "original_photo_image_attempts":len(all_photo),
       "independent_new_photo_observation_UUIDs":int(all_photo.observation_id.nunique()),
       "source_distinct_GSHHG_islands":int(all_photo.land_id.nunique()),
       "actual_numerical_9palette_flowerROI_success":int(good.sum()),
       "original_photo_technical_unknown":int((~good).sum()),
       "original_photo_acquisition_failures":int(all_photo.measurement_status.eq(
          "image_acquisition_failed").sum()),
       "source_photo_shard_receipts_verified":len(receipts),
       "original_FCP_model_sha_blob":FCP_SOURCE_BLOB,
       "original_real_photo_colour_file_SHA256":sha(result_file),
       "source_photo_manifest_SHA256":sha(source/"manifest.json"),
       "technical_unknown_not_assigned_an_unnatural_colour":True,
       "original_RGB_image_bytes_not_persisted":True,
       "source_island_and_site_prioritised_not_worldwide_probability":True,
       "original_bloom_native_flora_not_verified":True,
    }
    require(stats["actual_numerical_9palette_flowerROI_success"]+
            stats["original_photo_technical_unknown"]==len(all_photo),
            "Authentic colour successes and technical unknowns do not sum to source photos")
    (out/"ALAN_tier07_genuine_new_photo_ROI9_task_receipt.json").write_text(
        json.dumps(stats,sort_keys=True,indent=2)+"\n")
    print("ALAN_TIER07_26483_TRUE_ORIGINAL_PHOTO_ROI9_TASK",
          json.dumps({k:stats[k] for k in ("task_number","original_photo_image_attempts",
             "independent_new_photo_observation_UUIDs",
             "actual_numerical_9palette_flowerROI_success",
             "original_photo_technical_unknown","source_distinct_GSHHG_islands",
             "original_photo_acquisition_failures")},sort_keys=True),flush=True)
    return stats

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="command",required=True)
    a=sub.add_parser("prepare")
    a.add_argument("--schedule",type=Path,required=True)
    a.add_argument("--receipt",type=Path,required=True)
    a.add_argument("--task",type=int,required=True)
    a.add_argument("--out",type=Path,required=True)
    b=sub.add_parser("collect")
    b.add_argument("--source",type=Path,required=True)
    b.add_argument("--measured",type=Path,required=True)
    b.add_argument("--task",type=int,required=True)
    b.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.command=="prepare":
        prepare(args.schedule,args.receipt,args.task,args.out)
    else:
        collect(args.source,args.measured,args.task,args.out)
