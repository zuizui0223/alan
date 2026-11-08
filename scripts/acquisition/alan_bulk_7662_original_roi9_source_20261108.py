#!/usr/bin/env python3
"""High-throughput ALAN 7,662 genuine NEW island original-image ROI9 measurements.

The input is the SHA-authenticated, unmeasured 8,942-photo original source
schedule (35 disjoint 256-photo source batches), prepared blind to colour from
the 93,209,351-photo iNaturalist original index. Original batches 1-5
(1,280 photos) are ALREADY processed/in progress elsewhere and are excluded.

One warm unchanged FCP ROI-v4 YOLO11+EfficientSAM inference process handles
TWO adjacent source batches per job: 15 groups of up to 512 original photos,
each consisting of four original 128-photo work units (last one has 110).
This reduces frozen detector bootstrap from 30 separate jobs to 15.
Unknown colours are retained, no original images persisted, and no site-year
scientific cell inferred until 5 photo attempts, 3 ROI successes, 3 accepted
actual taxa. This arm remains island-prioritized, not globally representative.
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
    valid_nine, HUES, COMPONENTS, FRACTIONS,
)

SOURCE_SHA="e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd"
SOURCE_SCHEMA="ALAN_8942_new_original_island_UUIDs_progressive_256photo_source_v1"
FCP_SOURCE_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
MODEL_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
TAXON_SOURCE_ARM="alan_93m_tier06_new_island_observation_ROI_source_v1"
GROUPS=15
BATCH_FIRST=6
BATCH_LAST=35
N_EACH_SOURCE=256
LAST_BATCH_COUNT=238
EXPECTED_FIVE=1280
EXPECTED_REMAINING=7662
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
ROI_COLS=tuple("strict9_"+p for p in PALETTE)

def require(ok,message):
    if not bool(ok):
        raise ValueError("ALAN_7662_BULK_FROZEN_ROI9: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def group_batches(group):
    require(type(group) is int and 0<=group<GROUPS,
            "Original frozen source group outside 15 disjoint work groups")
    return (BATCH_FIRST+2*group,BATCH_FIRST+2*group+1)

def group_size(group):
    a,b=group_batches(group)
    return (N_EACH_SOURCE+
            (LAST_BATCH_COUNT if b==BATCH_LAST else N_EACH_SOURCE))

def source_read(schedule,receipt):
    require(sha(schedule)==SOURCE_SHA,
            "8,942 original photo source schedule byte hash drift")
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")==SOURCE_SCHEMA and
            r.get("source_schedule_csv_sha256")==SOURCE_SHA and
            r.get("novel_original_observation_UUIDs")==8942 and
            r.get("independent_GSHHG_islands")==369 and
            r.get("source_batches")==35 and
            r.get("first_batch_unchanged") is True and
            r.get("no_photo_images_downloaded_or_hue_ROI_inferred") is True,
            "Original 8,942 observation IDs not precolour source-authenticated")
    s=pd.read_csv(schedule,compression="gzip",dtype=str,keep_default_na=False)
    require(len(s)==8942 and s.photo_id.is_unique and
            s.observation_id.is_unique and
            s.colour_measured.astype(str).str.lower().eq("false").all() and
            s.selection_probability.eq("").all() and
            s.source_arm.eq(TAXON_SOURCE_ARM).all(),
            "Original source-only photo IDs, UUIDs, and acquisition design changed")
    require(s.photo_id.str.fullmatch(r"\d+").all() and
            s.photo_url_large.str.fullmatch(
             r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
            ).all(),"Original iNaturalist open-data source URLs malformed")
    sizes=s.groupby("source_batch").size().to_dict()
    require(len(sizes)==35 and
            all(sizes.get(str(k))==N_EACH_SOURCE for k in range(1,35)) and
            sizes.get("35")==LAST_BATCH_COUNT,
            "Original 8,942 selected source group cardinalities changed")
    require(s.source_batch.astype(int).le(5).sum()==EXPECTED_FIVE and
            s.source_batch.astype(int).gt(5).sum()==EXPECTED_REMAINING,
            "Previously measured 1,280 photos not disjoint from rest")
    return s,r

def select_group(source,group):
    a,b=group_batches(group)
    # Stable original source photo order retained inside both 256-slot waves.
    items=source.loc[source.source_batch.isin((str(a),str(b)))].copy()
    items["photo_order_in_batch"]=items.photo_order_in_batch.astype(int)
    items=items.sort_values(["source_batch","photo_order_in_batch"],
                            key=None,kind="stable")
    require(len(items)==group_size(group) and items.photo_id.is_unique and
            items.observation_id.is_unique,
            "Missing original photo IDs or repeated original UUIDs in group")
    require(not items.observation_id.isin(
       set(source.loc[source.source_batch.astype(int).le(5),"observation_id"])).any(),
       "A previously measured original observation UUID was routed again")
    require(items.wcvp_accepted_plant_name_id.ne("").all() and
            items.photo_license.ne("").all(),
            "Original source lost accepted taxon/photo license")
    items["measurement_id"]=items.photo_id.map(
        lambda x:f"alan_bulk_g{group:02d}_original_{x}")
    items["image_filename"]=items.measurement_id+".jpg"
    items["shard"]=(np.arange(len(items))//128)
    items["source_snapshot_date"]="2026-09-27"
    require(items.shard.nunique()==4 and items.shard.max()==3 and
            items.groupby("shard").size().iloc[:3].eq(128).all(),
            "Original four work units of frozen ROI-v4 group invalid")
    return items

def freeze(schedule,receipt,group,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwritten, previously signed new source pixels")
    source,_=source_read(schedule,receipt)
    d=select_group(source,group)
    (out/"shards").mkdir(parents=True)
    hashes={}
    for unit in range(4):
        file=out/"shards"/f"shard-{unit:03d}.csv"
        d.loc[d.shard.eq(unit)].to_csv(file,index=False,lineterminator="\n")
        hashes[str(unit)]=sha(file)
    contract={
      "schema":FCP_SOURCE_SCHEMA,
      "colour_measured":False,
      "source_arm":TAXON_SOURCE_ARM,
      "work_units":4,"worker_ids":[0,1,2,3],
      "new_photo_attempts_planned":len(d),
      "photo_shard_sha256":hashes,
      "original_source_8942_sha256":SOURCE_SHA,
      "source_group":group,
      "original_source_batches":list(group_batches(group)),
      "source_observation_UUIds_are_all_distinct":True,
      "original_previous_first1280_photo_and_UUIDs_excluded":True,
      "source_precolour_colour_measured":False,
      "sampling_island_prioritized_not_worldwide_random":True
    }
    (out/"manifest.json").write_text(
        json.dumps(contract,indent=2,sort_keys=True)+"\n")
    print("ALAN_SIGNED_BULK512_NEW_ISLAND_SOURCE",json.dumps({
       "group":group,"batches":list(group_batches(group)),
       "new_observations":len(d),
       "true_distinct_islands":int(d.land_id.nunique()),
       "four_work_unit_SHA256":hashes},sort_keys=True),flush=True)
    return contract

def collect(source,measured,group,out):
    source,measured,out=map(Path,(source,measured,out))
    require(group in range(GROUPS) and not out.exists(),
            "Original image source group invalid or output already exists")
    contract=json.loads((source/"manifest.json").read_text())
    require(contract.get("schema")==FCP_SOURCE_SCHEMA and
            contract["source_group"]==group and
            contract["original_source_8942_sha256"]==SOURCE_SHA and
            contract["worker_ids"]==[0,1,2,3] and
            contract["new_photo_attempts_planned"]==group_size(group),
            "Signed ALAN high-throughput original photo source was replaced")
    frames=[];receipts=[]
    for unit in range(4):
        frozen,_=worker.input_rows(source,unit)
        folder=measured/f"shard_{unit:03d}"
        result=json.loads((folder/"receipt.json").read_text())
        file=folder/"photos.csv"
        actual=pd.read_csv(file,dtype=str,keep_default_na=False)
        require(result.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
               and result["shard"]==unit and
               result["attempted"]==len(frozen) and
               result["source_sha256"]==contract["photo_shard_sha256"][str(unit)] and
               result["photos_sha256"]==sha(file) and
               len(actual)==len(frozen) and actual.photo_id.is_unique and
               actual.observation_id.is_unique and
               set(actual.photo_id)==set(frozen.photo_id) and
               set(actual.observation_id)==set(frozen.observation_id),
               "Genuine original photo FCP terminal outcome is not signed/disjoint")
        require(set(FRACTIONS)<=set(actual.columns),
                "Frozen FCP evaluator omitted original nine-flower-fraction fields")
        frames.append(actual);receipts.append(result)
    complete=pd.concat(frames,ignore_index=True)
    methods={json.dumps(r["method"],sort_keys=True) for r in receipts}
    require(len(complete)==group_size(group) and
            complete.photo_id.is_unique and
            complete.observation_id.is_unique and
            len(methods)==1,"Original high-throughput group has duplicated UUIDs/mixed model")
    method=receipts[0]["method"]
    require(method.get("measurement_git_blob")==MODEL_BLOB and
            method.get("location_blind") is True and
            method.get("persisted_image_pixels") is False,
            "Frozen YOLO11/EfficientSAM source/ephemerality changed")
    valid,nine=valid_nine(complete)
    complete["strict9"]=valid
    for k,name in enumerate(ROI_COLS):
        complete[name]=nine[:,k]
    require(np.allclose(
        complete.loc[valid,list(ROI_COLS)].sum(axis=1),1,atol=1e-6,rtol=0)
        and complete.loc[~valid,list(ROI_COLS)].isna().all().all(),
        "Actual original colour palette did not sum to one or unknowns imputed")
    require(complete.source_batch.isin(map(str,group_batches(group))).all(),
            "New ROI photo results attached to wrong source batch")
    require(complete.observed_on.str.match(r"\d{4}-\d{2}-\d{2}").all(),
            "Original source photo observation date missing")
    out.mkdir(parents=True)
    columns=["photo_id","observation_id","inat_taxon_id",
       "wcvp_accepted_plant_name_id","land_id","target_id","source_batch",
       "year","observed_on","latitude","longitude","source_arm",
       "morph","measurement_status","image_sha256","strict9",*ROI_COLS]
    require(set(columns)<=set(complete.columns),
            "Numeric outcome lacks original source metadata or actual ROI9 white")
    dest=out/f"real_new_original_island_ROI9_group{group:02d}.csv.gz"
    complete[columns].to_csv(dest,index=False,
        compression={"method":"gzip","mtime":0},lineterminator="\n")
    actual={"group":group,"original_batches":list(group_batches(group)),
       "original_source_photo_attempts":len(complete),
       "new_original_independent_observation_UUIDs":int(complete.observation_id.nunique()),
       "new_original_GSHHG_islands":int(complete.land_id.nunique()),
       "original_nine_palette_ROI9_true_success":int(valid.sum()),
       "original_photo_ROI9_unknown":int((~valid).sum()),
       "original_photo_acquisition_failed":int(
           complete.measurement_status.eq("image_acquisition_failed").sum()),
       "source_method_original_git_blob":MODEL_BLOB,
       "original_numeric_colour_file_SHA256":sha(dest),
       "source_precolour_manifest_SHA256":sha(source/"manifest.json"),
       "original_same_photo_UUIDs_reprocessed":0,
       "source_image_bytes_persisted":False,
       "unknown_photo_colour_was_not_imputed":True,
       "wild_native_open_flower_unverified":True,
       "source_island_priority_not_worldwide_probability_sample":True,
       "source_photo_technical_statuses":{
           str(k):int(v) for k,v in complete.measurement_status.value_counts().items()},
    }
    require(actual["original_nine_palette_ROI9_true_success"]+
            actual["original_photo_ROI9_unknown"]==len(complete),
            "Original photo unknown-plus-classified count mismatch")
    (out/"real_original_island_bulk512_group_receipt.json").write_text(
        json.dumps(actual,indent=2,sort_keys=True)+"\n")
    print("ALAN_7662_REAL_SOURCE_ROI9_BULK_GROUP",json.dumps(actual,sort_keys=True),flush=True)
    return actual

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="command",required=True)
    a=sub.add_parser("prepare")
    a.add_argument("--schedule",type=Path,required=True)
    a.add_argument("--receipt",type=Path,required=True)
    a.add_argument("--group",type=int,required=True)
    a.add_argument("--out",type=Path,required=True)
    b=sub.add_parser("collect")
    b.add_argument("--source",type=Path,required=True)
    b.add_argument("--measured",type=Path,required=True)
    b.add_argument("--group",type=int,required=True)
    b.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.command=="prepare":
        freeze(args.schedule,args.receipt,args.group,args.out)
    else:
        collect(args.source,args.measured,args.group,args.out)
