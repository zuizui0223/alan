#!/usr/bin/env python3
"""ALAN independently-new observation 256 genuine RGB ROI-v4 measurements.

Prepare exactly one frozen original-photo work queue from source SHA-pinned
93m tier06 UUID-disjoint observation data. Process on GitHub Actions with the
UNCHANGED YOLO11+EfficientSAM FCP ROI-v4 evaluator. Independently validate
the signed per-photo terminal outputs; use the nine original flower_fraction
components only for technically classified true original ROI source images.

No missing ROI or acquisition failure becomes white/coloured, no licensed
image bytes are saved, and no botanical anthesis/native/pollen claim is made.
Separate source design: do not overwrite the historical-main ALAN panel.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd

from scripts.acquisition import mois_full_index_real_original_photo_measure as worker

SOURCE_QUEUE_SHA="7fb7bbb60c556b3a2a31e10bda60715996e6738a662896c377f5f1569ff3bc56"
TIER_SOURCE_SCHEMA="ALAN_93m_original_tier06_256_new_observation_island_source_v1"
SOURCE_MODEL_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
N=256
PHOTO_UNIT=128
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
FRACTIONS=tuple("flower_fraction_"+x for x in PALETTE)
COMPONENTS=tuple("strict9_"+x for x in PALETTE)
HUES=("white","yellow_orange","red_pink","blue_purple")

def require(x,msg):
    if not bool(x): raise ValueError("ALAN_NEW256_REAL_RGB_ROI9: "+msg)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as file:
        for chunk in iter(lambda:file.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def validate_queue(queue,receipt):
    require(sha(queue)==SOURCE_QUEUE_SHA,"Original ALAN new256 precolour queue SHA drift")
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")==TIER_SOURCE_SCHEMA and
            r["signed_source_tier06"]["tier06_genuinely_indexed_photo_ids"]==597070 and
            r["first_256_blind_new_original_photo_queue"]["original_observations"]==N and
            r["source_photo_per_observation_selected_before_colour"] is True and
            r["source_image_pixels_downloaded"]==0,
            "No original new independent ALAN 93m source certification")
    q=pd.read_csv(queue,compression="gzip",dtype=str,keep_default_na=False)
    require(len(q)==N and q.photo_id.is_unique and q.observation_id.is_unique,
            "New photo selection repeated original image or observation UUID")
    require(q.land_id.nunique()==N and q.target_id.nunique()==N,
            "Pre-colour one-new-photo-per-island independent support changed")
    require(q.source_arm.eq("alan93m_new_island_tier06_conditional_priority").all() and
            q.selection_probability.eq("").all(),
            "Prioritised source must not claim known independent-worldwide inclusion π")
    require(q.photo_url_large.str.fullmatch(
          r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)").all(),
          "Unsafe/unexpected original S3 photographic URL or format")
    require(q.photo_id.str.fullmatch(r"\d+").all() and
            q.observed_on.str.match(r"\d{4}-\d{2}-\d{2}").all(),
            "Original photo ID or date not source authenticated")
    require(q.shard.isin(("0","1")).all() and
            q.groupby("shard").size().to_dict()=={"0":128,"1":128},
            "Original two ROI photo units changed")
    return q,r

def freeze(queue,receipt,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwriting source-only, precolour photo roster")
    table,r=validate_queue(queue,receipt)
    (out/"shards").mkdir(parents=True)
    checks={}
    for i in (0,1):
        frame=table.loc[table.shard.eq(str(i))].copy()
        target=out/"shards"/f"shard-{i:03d}.csv"
        frame.to_csv(target,index=False,lineterminator="\n")
        checks[str(i)]=sha(target)
    manifest={
        "schema":SOURCE_MODEL_SCHEMA,"colour_measured":False,
        "source_group":"alan_93m_tier06_novel_island_UUIDs",
        "source_arm":"alan93m_new_island_tier06_conditional_priority",
        "work_units":2,"worker_ids":[0,1],
        "new_photo_attempts_planned":N,
        "photo_shard_sha256":checks,
        "original_93m_release_manifest_sha256":
            r["underlying_source_SHA256"]["93m_original_tier_registry_SHA256"],
        "source_precolour_queue_sha256":SOURCE_QUEUE_SHA,
        "source_original_selected_observation_UUIDs_are_new_to_canonical":True,
        "selection_probability":None,
        "sampling_island_first_purposive_not_global_random":True,
        "no_original_image_colours_in_source":True,
        "photo_specific_wild_native_bloom_not_verified":True
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print("ALAN_NEW256_SOURCE_PRECOLOUR_FROZEN",json.dumps({
        "photos":N,"observation_uuids":N,
        "islands":int(table.land_id.nunique()),
        "work_unit_SHA256":checks},sort_keys=True),flush=True)
    return manifest

def valid_nine(frame):
    """Return actual existing numeric original ROI9 rather than inferred unknown hue."""
    raw=np.column_stack([
        pd.to_numeric(frame[col],errors="coerce").to_numpy(float)
        for col in FRACTIONS])
    good=np.isfinite(raw).all(axis=1)&(raw>=0).all(axis=1)&
        (raw<=1+1.e-6).all(axis=1)&(np.where(np.isfinite(raw),raw,0).sum(axis=1)>0)
    source4=frame.morph.isin(HUES).to_numpy()
    status=frame.measurement_status.eq("classified_four_state_morph").to_numpy()
    require(np.array_equal(source4,status),
            "Four-hue classification or technical unknown status drift")
    good=good&source4
    sums=np.where(np.isfinite(raw),raw,0).sum(axis=1)
    palette=np.full(raw.shape,np.nan)
    palette[good]=raw[good]/sums[good,None]
    return good,palette

def collect(source,measured,out):
    source=Path(source);measured=Path(measured);out=Path(out)
    require(not out.exists(),"Do not overwrite prior original 9-channel result")
    manifest=json.loads((source/"manifest.json").read_text())
    require(manifest.get("schema")==SOURCE_MODEL_SCHEMA and
            manifest["source_precolour_queue_sha256"]==SOURCE_QUEUE_SHA and
            manifest["worker_ids"]==[0,1],"New256 source manifest changed")
    frames=[];receipts=[]
    for unit in (0,1):
        signed, _=worker.input_rows(source,unit)
        p=measured/f"shard_{unit:03d}"
        rec=json.loads((p/"receipt.json").read_text())
        data=pd.read_csv(p/"photos.csv",dtype=str,keep_default_na=False)
        require(rec["schema"]=="mois_full_original93m_real_RGB_photo_unit_v1" and
                rec["shard"]==unit and rec["attempted"]==PHOTO_UNIT and
                rec["source_sha256"]==manifest["photo_shard_sha256"][str(unit)] and
                rec["photos_sha256"]==sha(p/"photos.csv"),
                "One genuine new original photo source ROI work unit not signed")
        require(len(data)==PHOTO_UNIT and
                data.photo_id.is_unique and data.observation_id.is_unique and
                set(data.photo_id)==set(signed.photo_id),
                "Source and measured original photo ids do not match")
        require(set(FRACTIONS)<=set(data.columns),
                "Fixed FCP worker did not retain nine measured flower ROI fractions")
        frames.append(data)
        receipts.append(rec)
    photo=pd.concat(frames,ignore_index=True)
    require(len(photo)==N and photo.photo_id.is_unique and
            photo.observation_id.is_unique and
            len({json.dumps(r["method"],sort_keys=True) for r in receipts})==1,
            "Duplicate photo/observation or mixed FCP measurement method")
    valid,palette=valid_nine(photo)
    original_four=photo.morph.isin(HUES)
    photo["strict9"]=valid
    for i,col in enumerate(COMPONENTS):
        photo[col]=palette[:,i]
    photo["strict9_white"]=photo["strict9_white"].where(photo.strict9)
    require(photo.loc[valid,"strict9_white"].between(0,1).all() and
            np.allclose(photo.loc[valid,list(COMPONENTS)].sum(axis=1).to_numpy(),
                        1,atol=1e-6,rtol=0),
            "Original signed numeric 9-fraction ROI not compositional")
    cannot=~valid
    require(photo.loc[cannot,list(COMPONENTS)].isna().all().all(),
            "Unknown original image was fabricated as a measured nine-colour hue")
    out.mkdir(parents=True)
    columns=["photo_id","observation_id","inat_taxon_id","source_arm",
         "observed_on","latitude","longitude","land_id","target_id","year",
         "morph","measurement_status","image_sha256","strict9","strict9_white",
         *COMPONENTS]
    for key in columns:
        require(key in photo.columns,"Original photo result missing source field "+key)
    filename=out/"new256_independent_island_actual_original_ROI9.csv.gz"
    photo[columns].to_csv(filename,compression="gzip",index=False)
    by_year=photo.groupby(["year","morph"],observed=True).size().reset_index(name="source_photo_count")
    by_year.to_csv(out/"new256_original_photo_outcomes_by_year.csv",index=False)
    stats={
        "schema":"ALAN_new256_original_island_unique_UUID_real_FCP_ROI9_v1",
        "original_source_photo_attempts":N,
        "independent_new_observation_UUIDs":int(photo.observation_id.nunique()),
        "genuinely_distinct_original_GSHHG_islands":int(photo.land_id.nunique()),
        "photo_years":sorted(map(int,photo.year.unique())),
        "original_four_hue_classified":int(original_four.sum()),
        "original_nine_fraction_strict_colour_success":int(valid.sum()),
        "original_photo_ROI_technically_unknown":int(cannot.sum()),
        "source_original_9colour_unknown_not_cast_to_white_or_coloured":True,
        "raw_original_image_bytes_persisted":False,
        "botanical_open_flower_native_wild_not_independently_verified":True,
        "original_S3_photo_URL_acquisition_and_FCP_model_used":True,
        "source_source_manifest_SHA256":sha(source/"manifest.json"),
        "original_photo_colour_output_SHA256":sha(filename),
        "fixed_FCP_original_method":receipts[0]["method"],
        "technical_photo_statuses":{str(k):int(v) for k,v in photo.measurement_status.value_counts().items()},
        "source_study_is_conditionally_island_prioritized":True,
        "no_preexisting_original_observation_UUID_reused":True,
    }
    (out/"ALAN_new256_original_photo_ROI9_source_receipt.json").write_text(
        json.dumps(stats,indent=2,sort_keys=True)+"\n")
    print("ALAN_NEW_256_TRUE_ORIGINAL_IMAGE_COLOUR",json.dumps(stats,sort_keys=True),flush=True)
    return stats

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    subs=p.add_subparsers(dest="command",required=True)
    a=subs.add_parser("prepare")
    a.add_argument("--queue",type=Path,required=True)
    a.add_argument("--receipt",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    b=subs.add_parser("collect")
    b.add_argument("--source",type=Path,required=True)
    b.add_argument("--measured",type=Path,required=True)
    b.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.command=="prepare":freeze(args.queue,args.receipt,args.out)
    else:collect(args.source,args.measured,args.out)
