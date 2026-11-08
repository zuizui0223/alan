#!/usr/bin/env python3
"""ALAN fixed progressive 8,942-authentic-novel-island-observation photo schedule.

Consume exactly the original SHA-verified 8,942 novel ALAN island observation
UUIDs with accepted WCVP taxon and real same-photo-year VIIRS + WorldPop data
from the executed 93m tier06 metadata source. Preserve the independently frozen
first 256 photos and their observed order byte-for-byte at the identity level.

Schedule the remaining 8,686 with round-robin per GSHHG island, prioritizing
landmasses with fewer previously scheduled original observations and sparse
ALAN original photo support. Exactly one original photo per original UUID.
DO NOT infer colour, flowering state, wild/native status, or obtain image bytes.
Never claim global random sampling after conditional island prioritization.
"""
from __future__ import annotations
import argparse,collections,hashlib,heapq,json
from pathlib import Path
import pandas as pd

EXPECTED_MASTER_SHA="d0b7f7d9f117036d05fbc7ff8b1e62847ef05e7bbc07efc659c9f6cf895d30f9"
EXPECTED_FIRST_SHA="7fb7bbb60c556b3a2a31e10bda60715996e6738a662896c377f5f1569ff3bc56"
EXPECTED_MASTER_N=8942
EXPECTED_FIRST_N=256
EXPECTED_ISLAND_N=369
PHOTO_SALT="ALAN_progressive_island_observation_order_20261108_v1"
SOURCE_MASTER_ARM="alan_93m_tier06_new_island_observation_ROI_source_v1"
SOURCE_FIRST_ARM="alan93m_new_island_tier06_conditional_priority"
OUT_SCHEMA="ALAN_8942_new_original_island_UUIDs_progressive_256photo_source_v1"
BATCH=256

def require(cond,message):
    if not bool(cond):raise ValueError("ALAN_PROGRESSIVE_ISLAND_ROI_SOURCE: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def stable_hash(tag,identity):
    return hashlib.sha256((PHOTO_SALT+"|"+tag+"|"+str(identity)).encode()).hexdigest()

def verify_source(master,first,*,check_sha=True):
    if check_sha:
        require(sha(master)==EXPECTED_MASTER_SHA,
                "Original 8,942 source gzip asset SHA changed")
        require(sha(first)==EXPECTED_FIRST_SHA,
                "Original first 256 genuine-photo queue SHA changed")
    m=pd.read_csv(master,compression="gzip",dtype=str,keep_default_na=False)
    f=pd.read_csv(first,compression="gzip",dtype=str,keep_default_na=False)
    require(len(m)==EXPECTED_MASTER_N and len(f)==EXPECTED_FIRST_N,
            "Actual novel-source photo count changed")
    for obj in (m,f):
        require(obj.photo_id.is_unique and obj.observation_id.is_unique and
                obj.photo_id.str.fullmatch(r"\d+").all() and
                obj.observation_id.ne("").all(),
                "Repeated original photo/observation IDs in novel island source")
        require(obj.land_id.notna().all() and obj.target_id.notna().all() and
                obj.photo_url_large.str.startswith(
                    "https://inaturalist-open-data.s3.amazonaws.com/photos/").all(),
                "Missing GSHHG identity or unsafe frozen source URL")
    require(m.land_id.nunique()==EXPECTED_ISLAND_N and
            f.land_id.nunique()==EXPECTED_FIRST_N,
            "Genuine distinct island support changed")
    require(m.source_arm.eq(SOURCE_MASTER_ARM).all() and
            f.source_arm.eq(SOURCE_FIRST_ARM).all() and
            m.selection_probability.eq("").all() and
            f.selection_probability.eq("").all(),
            "Source design inclusion probabilities were silently invented")
    mset=m.set_index("observation_id",verify_integrity=True)
    require(f.observation_id.isin(mset.index).all(),
            "Frozen first 256 no longer belong to original full 8,942 source")
    q=mset.loc[f.observation_id]
    for attr in ("photo_id","land_id","target_id","inat_taxon_id","year",
                 "photo_url_large"):
        require(q[attr].tolist()==f[attr].tolist(),
                "Same original photo/observation attributes altered: "+attr)
    return m,f

def schedule(master, first):
    """All 8,942 unique real observation UUIDs exactly once, fixed first256."""
    first_ids=set(first.observation_id)
    rest=master.loc[~master.observation_id.isin(first_ids)].copy()
    require(len(rest)==EXPECTED_MASTER_N-EXPECTED_FIRST_N and
            not set(rest.photo_id)&set(first.photo_id),
            "Already-selected original photos leaked into future batches")
    # Put every island's pending observations in reproducible hue-blind order.
    groups={}
    for land,d in rest.groupby("land_id",sort=True):
        group=d.copy()
        group["_rank"]=group.observation_id.map(
             lambda obs:stable_hash("observation",obs))
        group=group.sort_values(["_rank","observation_id"],kind="stable")
        groups[str(land)]=iter(group.drop(columns="_rank").to_dict("records"))
    taken=collections.Counter(first.land_id.astype(str))
    original_counts={str(k):int(v) for k,v in
                     master.groupby("land_id").existing_strict_island_cell_years.first().items()}
    heap=[]
    for land in groups:
        heapq.heappush(heap,(
            taken[land],original_counts[land],
            stable_hash("landmass",land),land))
    order=[]
    while heap:
        used,prior,rank,land=heapq.heappop(heap)
        require(used==taken[land],"Inconsistent original per-island count")
        item=next(groups[land],None)
        if item is None:continue
        order.append(item)
        taken[land]+=1
        heapq.heappush(heap,(
            taken[land],original_counts[land],rank,land))
    require(len(order)==len(rest),"Lost actual novel source observations")
    initial=master.set_index("observation_id").loc[first.observation_id].reset_index()
    suffix=pd.DataFrame(order,columns=master.columns)
    planned=pd.concat([initial,suffix],ignore_index=True)
    require(len(planned)==EXPECTED_MASTER_N and
            planned.observation_id.is_unique and planned.photo_id.is_unique and
            set(planned.photo_id)==set(master.photo_id) and
            set(planned.observation_id)==set(master.observation_id),
            "Progressive schedule dropped/duplicated original UUIDs")
    require(planned.loc[:255,"observation_id"].tolist()==first.observation_id.tolist(),
            "Existing first256 source-photo identities were rescheduled")
    planned["source_batch"]=(planned.index//BATCH)+1
    planned["photo_order_in_batch"]=(planned.index%BATCH)
    planned["colour_measured"]=False
    require(planned.source_batch.max()==35 and
            planned.groupby("source_batch").size().iloc[:-1].eq(BATCH).all(),
            "Source 8,942 must form exactly 34 full and one 238-photo batch")
    return planned

def run(master,first,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwriting signed source photo schedule")
    m,f=verify_source(master,first)
    p=schedule(m,f)
    out.mkdir(parents=True)
    source=out/"progressive_8942_new_island_source_schedule.csv.gz"
    p.to_csv(source,index=False,compression={"method":"gzip","mtime":0},
             lineterminator="\n")
    batches=p.groupby("source_batch",observed=True).agg(
        genuine_observation_UUIDs=("observation_id","nunique"),
        new_GSHHG_islands=("land_id","nunique"),
        distinct_photo_cell_years=("target_id","nunique")).reset_index()
    batches.to_csv(out/"progressive_256_photo_island_batch_support.csv",index=False)
    receipt={
       "schema":OUT_SCHEMA,
       "first_batch_unchanged":True,
       "source_previous_island_queue_photo_sha256":EXPECTED_FIRST_SHA,
       "original_signed_source_8942_sha256":EXPECTED_MASTER_SHA,
       "novel_original_observation_UUIDs":len(p),
       "novel_original_photo_ids":int(p.photo_id.nunique()),
       "independent_GSHHG_islands":int(p.land_id.nunique()),
       "batch_size":BATCH,"source_batches":int(p.source_batch.max()),
       "full_256_photo_batches":int(batches.genuine_observation_UUIDs.eq(BATCH).sum()),
       "final_partial_batch_photos":int(batches.genuine_observation_UUIDs.iloc[-1]),
       "second_batch_genuine_independent_observations":int(
           batches.genuine_observation_UUIDs.iloc[1]),
       "second_batch_independent_islands":int(
           batches.new_GSHHG_islands.iloc[1]),
       "source_schedule_csv_sha256":sha(source),
       "source_batch_support_sha256":sha(
           out/"progressive_256_photo_island_batch_support.csv"),
       "all_photo_ids_are_source_preselected_no_colour":True,
       "no_photo_images_downloaded_or_hue_ROI_inferred":True,
       "no_source_sample_inclusion_probability_invented":True,
       "cannot_stack_island_targeted_batch_as_unbiased_global_flora":True,
    }
    (out/"progressive_island_original_photo_schedule_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("ALAN_SIGNED_PROGRESSIVE_8942_NEW_ISLAND_SOURCE",
          json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master",type=Path,required=True)
    p.add_argument("--first",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.master,a.first,a.out)
