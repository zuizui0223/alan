#!/usr/bin/env python3
"""ALAN true-original-photo tier07 1024 new island UUID pilot ROI9 acquisition.

All original photo/observation IDs are pre-colour selected from signed 27,507
source: 256 geographically spread island reference, 768 source site-year
densification. The actual images are acquired with the EXACT unchanged
FCP YOLO11/EfficientSAM original ROI-v4 evaluator in two warm 512-photo
subgroups within ONE GitHub Actions job, retaining 8 signed 128-photo worker
receipts and one source manifest. Technical unknowns remain unknown.
No original RGB image bytes are included in the final output.
"""
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.acquisition import mois_full_index_real_original_photo_measure as worker
from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import (
    valid_nine,PALETTE,COMPONENTS,FRACTIONS,HUES,
)

SOURCE_SHA="7198455819f1ac260a6cbeacd9f4536c30d766595daa64dc5205c75e9d5a1c62"
SOURCE_SCHEMA="ALAN_tier07_27507_new_original_island_1024photo_source_schedule_v1"
WORKER_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
MODEL_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
SOURCE_ARM="alan_93m_tier07_new_island_observation_ROI_source_v1"
EXPECTED=1024

def require(value,message):
    if not bool(value):raise ValueError("ALAN_TIER07_REAL1024_ROI9: "+message)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def original_source(schedule,receipt):
    require(sha(schedule)==SOURCE_SHA,"Pinned full27507 original UUID source photo SHA drift")
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")==SOURCE_SCHEMA and
            r.get("precolour_source_photo_plan_SHA256")==SOURCE_SHA and
            r.get("original_source_accepted_new_observation_UUIDs")==27507 and
            r.get("disjoint_genuine_source_groups")==27 and
            r.get("source_photo_selection_did_not_use_colour") is True and
            r.get("original_new_image_bytes_downloaded")==0,
            "Original independent tier07 source did not pass precolour source gates")
    d=pd.read_csv(schedule,compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==27507 and d.photo_id.is_unique and d.observation_id.is_unique and
            d.work_group.astype(int).between(0,26).all() and
            d.source_arm.eq(SOURCE_ARM).all() and
            d.source_colour_measured.astype(str).str.lower().eq("false").all(),
            "Original 27,507 real source photo IDs/observations changed")
    group=d.loc[d.work_group.eq("0")].copy()
    require(len(group)==EXPECTED and group.photo_id.is_unique and
            group.observation_id.is_unique and group.land_id.nunique()==284 and
            group.target_id.nunique()==324 and
            group.source_priority.value_counts().to_dict()=={
                "coverage_reference":256,"siteyear_densification":768},
            "First1024 colour-blind photo priority pilot source drift")
    require(group.photo_url_large.str.fullmatch(
            r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
            ).all() and group.photo_id.str.fullmatch(r"\d+").all(),
            "Original iNaturalist open-data source URL malformed")
    return group,r

def freeze(schedule,receipt,out):
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite the SHA-locked original 1024 source pilot")
    d,r=original_source(schedule,receipt)
    d["measurement_id"]=d.photo_id.map(lambda x:"alan_tier07_pilot_"+x)
    d["image_filename"]=d.measurement_id+".jpg"
    d["source_snapshot_date"]="2026-09-27"
    d["shard"]=np.arange(len(d))//128
    require(d.shard.nunique()==8 and d.shard.value_counts().eq(128).all(),
            "Original pilot 8×128 work unit packing invalid")
    (out/"shards").mkdir(parents=True)
    checks={}
    for unit in range(8):
        file=out/"shards"/f"shard-{unit:03d}.csv"
        d.loc[d.shard.eq(unit)].to_csv(file,index=False,lineterminator="\n")
        checks[str(unit)]=sha(file)
    m={
      "schema":WORKER_SCHEMA,"colour_measured":False,
      "source_group":"alan_original_tier07_new1024_island_breadth_density_pilot",
      "source_arm":SOURCE_ARM,
      "work_units":8,"worker_ids":list(range(8)),
      "new_photo_attempts_planned":EXPECTED,
      "photo_shard_sha256":checks,
      "original_source_27507_SHA256":SOURCE_SHA,
      "original_source_pilot_group":0,
      "original_uuid_independent":True,
      "photo_original_source_colour_unknown":True,
      "source_shards_selected_blind_to_actual_photo_pixels":True,
    }
    (out/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
    print("ALAN_TIER07_FIRST1024_ORIGINAL_UUID_SOURCE_FROZEN",json.dumps({
        "source_photos":EXPECTED,"islands":int(d.land_id.nunique()),
        "site_years":int(d.target_id.nunique()),"shard_source_SHA256":checks},
        sort_keys=True),flush=True)
    return m

def collect(source,measured,out):
    source,measured,out=map(Path,(source,measured,out))
    require(not out.exists(),"Do not overwrite already measured original-photo ROI9 pilot")
    manifest=json.loads((source/"manifest.json").read_text())
    require(manifest.get("schema")==WORKER_SCHEMA and
            manifest.get("original_source_27507_SHA256")==SOURCE_SHA and
            manifest.get("worker_ids")==list(range(8)),
            "Original pilot source manifest absent or changed")
    source_frames=[];receipts=[]
    for unit in range(8):
        frozen,_=worker.input_rows(source,unit)
        folder=measured/f"shard_{unit:03d}"
        r=json.loads((folder/"receipt.json").read_text())
        file=folder/"photos.csv"
        df=pd.read_csv(file,dtype=str,keep_default_na=False)
        require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
            and r.get("shard")==unit and r.get("attempted")==128 and
            r.get("source_sha256")==manifest["photo_shard_sha256"][str(unit)] and
            r.get("photos_sha256")==sha(file) and
            len(df)==128 and df.photo_id.is_unique and
            df.observation_id.is_unique and set(df.photo_id)==set(frozen.photo_id) and
            set(df.observation_id)==set(frozen.observation_id),
            "Actual original RGB FCP photo/UUID receipt does not match source")
        source_frames.append(df);receipts.append(r)
    actual=pd.concat(source_frames,ignore_index=True)
    methods={json.dumps(r["method"],sort_keys=True) for r in receipts}
    require(len(actual)==EXPECTED and actual.photo_id.is_unique and
            actual.observation_id.is_unique and len(methods)==1,
            "FCP model mixed or original source pilot photo UUIDs are not independent")
    method=receipts[0]["method"]
    require(method.get("measurement_git_blob")==MODEL_BLOB and
            method.get("location_blind") is True and
            method.get("persisted_image_pixels") is False,
            "Frozen source photo ROI FCP model was modified")
    valid,nine=valid_nine(actual)
    actual["strict9"]=valid
    for i,col in enumerate(COMPONENTS):
        actual[col]=nine[:,i]
    require(np.allclose(actual.loc[valid,list(COMPONENTS)].sum(axis=1),1,
                        atol=1.e-6,rtol=0) and
            actual.loc[~valid,list(COMPONENTS)].isna().all().all(),
            "Original nine-hue numeric ROI mass invalid or unknown was imputed")
    original_colour=actual.loc[valid].copy()
    require(original_colour.source_priority.value_counts().index.isin(
      ["coverage_reference","siteyear_densification"]).all(),
      "Source priority secretly changed after looking at image colour")
    out.mkdir(parents=True)
    fields=["photo_id","observation_id","inat_taxon_id",
      "wcvp_accepted_plant_name_id","land_id","target_id","year",
      "source_priority","source_arm","morph","measurement_status",
      "image_sha256","strict9",*COMPONENTS]
    require(set(fields)<=set(actual.columns),"Original photo RGB results lost source identity")
    file=out/"ALAN_tier07_first1024_original_floral_roi9.csv.gz"
    actual[fields].to_csv(file,index=False,
           compression={"method":"gzip","mtime":0},lineterminator="\n")
    actual["real_success"]=valid
    sites=actual.groupby("target_id",observed=True).agg(
        attempts=("observation_id","nunique"),
        photo_success=("real_success","sum"),
        accepted_species=("wcvp_accepted_plant_name_id","nunique"))
    real=actual.loc[valid].groupby("target_id",observed=True).wcvp_accepted_plant_name_id.nunique()
    sites["measured_accepted_species"]=real.reindex(sites.index,fill_value=0)
    strict=sites.attempts.ge(5)&sites.photo_success.ge(3)&sites.measured_accepted_species.ge(3)
    expected_source=sites.attempts.ge(5)&sites.accepted_species.ge(3)
    require(int(expected_source.sum())==76,"Original 76 photo-opportunity-qualified source sites drifted")
    report={
      "schema":"ALAN_tier07_1024_real_original_FCP_ROI9_pilot_v1",
      "original_image_attempts":EXPECTED,
      "new_independent_observation_UUIDs":int(actual.observation_id.nunique()),
      "source_distinct_GSHHG_islands":int(actual.land_id.nunique()),
      "original_9palette_numeric_ROI_success":int(valid.sum()),
      "original_technical_ROI_unknown":int((~valid).sum()),
      "original_RGB_photo_acquisition_failures":int(actual.measurement_status.eq(
          "image_acquisition_failed").sum()),
      "original_site_years":int(len(sites)),
      "source_5photo_3accepted_potential_site_years":int(expected_source.sum()),
      "actual_strict_5photo_3true_ROI9_3accepted_species_sites":int(strict.sum()),
      "source_reference_photo_attempts":256,
      "source_siteyear_targeted_photo_attempts":768,
      "real_original_ROI9_photo_outcomes_sha256":sha(file),
      "original_source_schedule_sha256":SOURCE_SHA,
      "original_source_manifest_sha256":sha(source/"manifest.json"),
      "original_frozen_FCP_ROI9_method":method,
      "per_128_photo_unit_receipts_verified":8,
      "source_original_image_bytes_persisted":False,
      "source_photo_unknown_not_assigned_white":True,
      "new_source_island_and_site_targeted_not_worldwide_random":True
    }
    require(report["original_9palette_numeric_ROI_success"]+
            report["original_technical_ROI_unknown"]==EXPECTED,
            "Actual original source photo outcomes do not sum to selected 1024")
    (out/"ALAN_tier07_first1024_original_roi9_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_TIER07_REAL_FIRST1024_ROI9",json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    actions=p.add_subparsers(dest="command",required=True)
    q=actions.add_parser("prepare")
    q.add_argument("--schedule",type=Path,required=True)
    q.add_argument("--receipt",type=Path,required=True)
    q.add_argument("--out",type=Path,required=True)
    q=actions.add_parser("collect")
    q.add_argument("--source",type=Path,required=True)
    q.add_argument("--measured",type=Path,required=True)
    q.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.command=="prepare":freeze(a.schedule,a.receipt,a.out)
    else:collect(a.source,a.measured,a.out)
