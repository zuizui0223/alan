#!/usr/bin/env python3
"""ALAN source-locked ORIGINAL RGB flower ROI9 for novel-island batches 2–5.

The 8,942 observation UUIDs are pre-colour, disjoint and SHA-verified in the
ALAN progressive source artifact. Each batch is 256 NEW original observations,
one photo per UUID, selected without looking at colour; island-balanced but
NOT a globally representative probability sample. The same original FCP
YOLO11/EfficientSAM model and unmodified source-image fetcher as batch 1 run
in a separate reproducible GitHub Actions worker per batch. No image bytes
persist and technical unknowns stay unknown.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.acquisition import mois_full_index_real_original_photo_measure as worker
from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import (
    valid_nine,HUES,COMPONENTS,FRACTIONS,
)

SCHEDULE_SHA="e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd"
SCHEDULE_SCHEMA="ALAN_8942_new_original_island_UUIDs_progressive_256photo_source_v1"
SOURCE_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
MODEL_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
SOURCE_STUDY_ARM="alan_93m_tier06_new_island_observation_ROI_source_v1"
BATCHES=(2,3,4,5)
SIZE=256

def require(ok,msg):
    if not bool(ok):
        raise ValueError("ALAN_REAL_NEW_ISLAND_ROI9_BATCHES02TO05: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def select_batch(frame,number):
    """Defensive identity-only split; data remain physically unmeasured."""
    require(number in BATCHES,"Invalid source-pinned original colour batch")
    q=frame.copy()
    for col in ("photo_id","observation_id","land_id","target_id","source_batch",
                "photo_url_large","photo_license","source_arm",
                "selection_probability","inat_taxon_id","year","observed_on"):
        require(col in q,"Missing original source observation field "+col)
    require(q.photo_id.is_unique and q.observation_id.is_unique,
            "Full prior source contains reused photo IDs or observation UUIDs")
    require(q.photo_id.str.fullmatch(r"\d+").all() and
            q.selection_probability.eq("").all() and
            q.source_arm.eq(SOURCE_STUDY_ARM).all(),
            "Prioritize-selected original source identity/probability changed")
    require(q.photo_url_large.str.fullmatch(
         r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
    ).all(),"Original source photo URL or extension invalid")
    require(q.colour_measured.astype(str).str.lower().eq("false").all(),
            "Pre-colour scheduled source claims measured photos")
    require(q.source_batch.isin(tuple(str(i) for i in range(1,36))).all(),
            "Original progressive source batch out of range")
    a=q.loc[q.source_batch.eq(str(number))].copy()
    b=q.loc[q.source_batch.eq("1")]
    require(len(a)==SIZE and len(b)==SIZE and
            not set(a.photo_id)&set(b.photo_id) and
            not set(a.observation_id)&set(b.observation_id),
            "Original source batch size or overlap with measured batch1")
    require(a.year.isin(("2019","2021","2023","2024","2025")).all(),
            "Unregistered source photo observation year")
    if number==2:
        require(a.land_id.nunique()==249,
                "Signed batch2 original distinct GSHHG islands changed")
    a["measurement_id"]=a.photo_id.map(lambda s:f"alan93m_b{number:02d}_{s}")
    a["image_filename"]=a.measurement_id+".jpg"
    a["source_snapshot_date"]="2026-09-27"
    a["shard"]=(np.arange(SIZE)//128).astype(int)
    require(a.measurement_id.is_unique and a.shard.value_counts().to_dict()=={0:128,1:128},
            "Original photo 128-shard source packing invalid")
    return a

def verify_full_source(schedule,receipt):
    require(sha(schedule)==SCHEDULE_SHA,"Exact frozen 8,942-photo original source SHA changed")
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")==SCHEDULE_SCHEMA and
            r.get("source_schedule_csv_sha256")==SCHEDULE_SHA and
            r.get("novel_original_observation_UUIDs")==8942 and
            r.get("independent_GSHHG_islands")==369 and
            r.get("first_batch_unchanged") is True and
            r.get("source_batches")==35 and
            r.get("no_photo_images_downloaded_or_hue_ROI_inferred") is True,
            "Not the signed original ALAN 8,942-photo source schedule")
    q=pd.read_csv(schedule,compression="gzip",dtype=str,keep_default_na=False)
    require(len(q)==8942 and q.photo_id.is_unique and
            q.observation_id.is_unique,"Original ALAN source schedule is incomplete")
    return q,r

def freeze(schedule,receipt,batch,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwrite of original source/ROI manifest")
    frame,r=verify_full_source(schedule,receipt)
    q=select_batch(frame,batch)
    (out/"shards").mkdir(parents=True)
    digest={}
    for part in (0,1):
        rows=q.loc[q.shard.eq(part)]
        file=out/"shards"/f"shard-{part:03d}.csv"
        rows.to_csv(file,index=False,lineterminator="\n")
        digest[str(part)]=sha(file)
    signed={
      "schema":SOURCE_SCHEMA,"colour_measured":False,
      "source_group":"alan_93m_tier06_novel_island_UUIDs",
      "source_arm":SOURCE_STUDY_ARM,
      "work_units":2,"worker_ids":[0,1],
      "new_photo_attempts_planned":SIZE,
      "photo_shard_sha256":digest,
      "prior_original_source_schedule_sha256":SCHEDULE_SHA,
      "original_93m_full_tier_registry_sha256":
           r["original_signed_source_8942_sha256"],
      "original_new_photo_batch_number":batch,
      "original_island_count_in_batch":int(q.land_id.nunique()),
      "source_observation_photo_UUIDs_unique":True,
      "selection_probability":None,
      "island_prioritized_not_a_random_global_flora_sample":True,
      "source_photo_pixels_downloaded":0,
      "photo_specific_wild_native_bloom_not_verified":True
    }
    (out/"manifest.json").write_text(json.dumps(signed,sort_keys=True,indent=2)+"\n")
    print("ALAN_NEW_NOVEL_ISLAND_BATCH_SOURCE_FROZEN",json.dumps({
        "batch":batch,"photos":SIZE,"distinct_islands":int(q.land_id.nunique()),
        "original_photo_source_sha256":SCHEDULE_SHA,
        "work_shards_SHA256":digest},sort_keys=True),flush=True)
    return signed

def collect(source,measurements,batch,out):
    source,measurements,out=map(Path,(source,measurements,out))
    require(batch in BATCHES and not out.exists(),"Refuse replacing signed actual new-island pixels")
    m=json.loads((source/"manifest.json").read_text())
    require(m.get("schema")==SOURCE_SCHEMA and m["worker_ids"]==[0,1] and
            m["prior_original_source_schedule_sha256"]==SCHEDULE_SHA and
            m["original_new_photo_batch_number"]==batch,
            "Original source schedule/colour measurement batch mismatch")
    frames=[];receipts=[]
    for part in (0,1):
        frozen,_=worker.input_rows(source,part)
        folder=measurements/f"shard_{part:03d}"
        r=json.loads((folder/"receipt.json").read_text())
        file=folder/"photos.csv"
        d=pd.read_csv(file,dtype=str,keep_default_na=False)
        require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
             and r["shard"]==part and r["attempted"]==128
             and r["photos_sha256"]==sha(file)
             and r["source_sha256"]==m["photo_shard_sha256"][str(part)],
             "Real original-photo ROI worker terminal SHA mismatch")
        require(len(d)==128 and d.photo_id.is_unique and
                d.observation_id.is_unique and
                set(d.photo_id)==set(frozen.photo_id) and
                set(d.observation_id)==set(frozen.observation_id),
                "Original-photo 128-shard IDs or observation UUIDs changed")
        require(set(FRACTIONS)<=set(d.columns),
                "Original FCP evaluator did not retain the nine numeric flower fractions")
        frames.append(d);receipts.append(r)
    full=pd.concat(frames,ignore_index=True)
    require(len(full)==SIZE and full.photo_id.is_unique and
            full.observation_id.is_unique and
            full.source_batch.astype(str).eq(str(batch)).all(),
            "Original source observation collision or wrong frozen batch")
    methods={json.dumps(r["method"],sort_keys=True) for r in receipts}
    require(len(methods)==1,"Original ROI workers changed method mid-batch")
    method=receipts[0]["method"]
    require(method.get("measurement_git_blob")==MODEL_BLOB and
            method.get("location_blind") is True and
            method.get("persisted_image_pixels") is False,
            "Original FCP source blob/location blinding/image ephemerality changed")
    valid,raw=valid_nine(full)
    full["strict9"]=valid
    for i,column in enumerate(COMPONENTS):
        full[column]=raw[:,i]
    full["strict9_white"]=pd.to_numeric(full["strict9_white"],errors="coerce").where(
       full.strict9)
    require(full.loc[valid,"strict9_white"].between(0,1).all() and
            np.allclose(full.loc[valid,list(COMPONENTS)].sum(axis=1),1,atol=1e-6,rtol=0)
            and full.loc[~valid,list(COMPONENTS)].isna().all().all(),
            "Original continuous9 hue unknown was imputed or numeric composition invalid")
    out.mkdir(parents=True)
    fields=["photo_id","observation_id","inat_taxon_id","land_id","target_id",
       "source_batch","year","observed_on","latitude","longitude","source_arm",
       "morph","measurement_status","image_sha256","strict9","strict9_white",
       *COMPONENTS]
    require(set(fields)<=set(full.columns),"Original photo output lacks provenance keys")
    file=out/f"ALAN_original_novel_island_photo_ROI9_batch_{batch:02d}.csv.gz"
    full[fields].to_csv(file,index=False,compression={"method":"gzip","mtime":0})
    status={str(k):int(v) for k,v in full.measurement_status.value_counts().items()}
    stats={
       "schema":"ALAN_novel_source_island_original256_FCP_ROI9_batch_v1",
       "source_batch":batch,
       "original_images_downloaded_and_measured":SIZE,
       "new_distinct_observation_UUIDs":int(full.observation_id.nunique()),
       "actual_distinct_GSHHG_islands":int(full.land_id.nunique()),
       "numeric_original_nine_palette_success":int(valid.sum()),
       "technically_unknown_original_ROI9":int((~valid).sum()),
       "four_hue_source_classified":int(full.morph.isin(HUES).sum()),
       "image_acquisition_failures":int(full.measurement_status.eq("image_acquisition_failed").sum()),
       "technical_photo_statuses":status,
       "source_photo_years":sorted(map(int,full.year.unique())),
       "original_frozen_model_method":method,
       "original_photo_source_schedule_sha256":SCHEDULE_SHA,
       "source_manifest_sha256":sha(source/"manifest.json"),
       "measured_original_colours_sha256":sha(file),
       "source_images_retained_in_artifact":False,
       "technical_unknown_photo_not_assigned_an_ROI9_colour":True,
       "botanical_open_native_flowers_not_independently_verified":True,
       "source_targeted_island_not_probability_complete_global":True,
       "no_original_prior_photo_or_observation_UUIDs_reused":True
    }
    (out/"ALAN_original_novel_island_ROI9_batch_receipt.json").write_text(
       json.dumps(stats,indent=2,sort_keys=True)+"\n")
    print("ALAN_TRUE_ORIGINAL_NEW_NOVEL_ISLAND_ROI9",json.dumps({
       k:stats[k] for k in ("source_batch","original_images_downloaded_and_measured",
          "new_distinct_observation_UUIDs","actual_distinct_GSHHG_islands",
          "numeric_original_nine_palette_success","technically_unknown_original_ROI9",
          "image_acquisition_failures","technical_photo_statuses")
    },sort_keys=True),flush=True)
    return stats

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="command",required=True)
    a=sub.add_parser("prepare")
    a.add_argument("--schedule",type=Path,required=True)
    a.add_argument("--receipt",type=Path,required=True)
    a.add_argument("--batch",type=int,required=True)
    a.add_argument("--out",type=Path,required=True)
    b=sub.add_parser("collect")
    b.add_argument("--source",type=Path,required=True)
    b.add_argument("--measured",type=Path,required=True)
    b.add_argument("--batch",type=int,required=True)
    b.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.command=="prepare":
        freeze(args.schedule,args.receipt,args.batch,args.out)
    else:
        collect(args.source,args.measured,args.batch,args.out)
