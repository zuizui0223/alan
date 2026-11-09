#!/usr/bin/env python3
"""Precolour SHA-signed island/area/light balanced tier08 original-photo 24,576 wave.

Select 24,576 immutable original iNaturalist S3 photo IDs from the 56,440
tier08 original-year VIIRS/environment island opportunities, AFTER removing
all 910 earlier tier08 source photos, all historical canonical source photo IDs,
and source tier06 / tier07 photo IDs. Historical numeric observation IDs and
new UUIDs are not directly exchangeable: only exact PHOTO ID nonoverlap is
claimed. Colour-blind selection preserves two same-island/year different-VIIRS
sites when available, representation of all eligible island area thirds, and
accepted-taxon variety. Three fixed, nonadaptively sampled 8192-photo waves.
The first 24,576 represent photo attempts, not real numerical ROI success.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

SCHEMA="ALAN_original93m_tier08_additional24576_source_sha_photoid_filtered_v1"
MASTER_SHA="b34a8a564dc9ba3f02934e1e04077b4655cf8a9333ecefcb6bc93c5d9dc757fd"
PREVIOUS910_SHA="9d2e60b9c9bdb2fb56c1c17dd69fca5b6c606d1e64d707020d070dd82e49a552"
TIER06_SHA="e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd"
TIER07_SHA="7198455819f1ac260a6cbeacd9f4536c30d766595daa64dc5205c75e9d5a1c62"
WAVES=(8192,16384,24576)
COUNT=24576
PHOTO_UNIT=128
WORK_GROUP_PHOTOS=512
SHA_SEED="ALAN_TIER08_24576_FROZEN_AREA_LIGHT_SITEYEAR_20261109_v1"
ARM="alan_93m_tier08_new_independent_island_original_photo"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_TIER08_24576_SOURCE: "+msg)
def sha(path):
    return panel.sha(path)
def rank(category,ident):
    return hashlib.sha256((SHA_SEED+"|"+category+"|"+str(ident)).encode()).hexdigest()
def load_verified(master,previous910,canonical,source06,source07):
    for label,path,digest in (
       ("tier08_master",master,MASTER_SHA),("tier08_910",previous910,PREVIOUS910_SHA),
       ("historical_canonical",canonical,panel.PHOTO_SHA),
       ("tier06",source06,TIER06_SHA),("tier07",source07,TIER07_SHA)):
        require(sha(path)==digest,"Changed SHA source original photo roster: "+label)
    m=pd.read_csv(master,compression="gzip",dtype=str,keep_default_na=False)
    prev=pd.read_csv(previous910,compression="gzip",dtype=str,keep_default_na=False,
                     usecols=["photo_id","observation_id"])
    hist=pd.read_csv(canonical,compression="gzip",dtype=str,
                     usecols=["photo_id","observation_id"],low_memory=False)
    s6=pd.read_csv(source06,compression="gzip",dtype=str,
                   usecols=["photo_id","observation_id"])
    s7=pd.read_csv(source07,compression="gzip",dtype=str,
                   usecols=["photo_id","observation_id"])
    require(len(m)==56440 and m.photo_id.is_unique and m.observation_id.is_unique
            and m.target_id.nunique()==13095 and m.land_id.nunique()==695
            and len(prev)==910 and len(s6)==8942 and len(s7)==27507,
            "Old actual 93m source opportunities or prior colour-image cohorts changed")
    require(m.source_arm.eq(ARM).all() and
            m.photo_id.str.fullmatch(r"\d+").all() and
            m.wcvp_accepted_plant_name_id.ne("").all() and
            m.photo_url_large.str.fullmatch(
                r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
            ).all() and m.year.isin(("2019","2021","2023","2024","2025")).all(),
            "Original 93m photo URL, source year, taxon or ID not trusted")
    previously=pd.concat([hist,prev,s6,s7],ignore_index=True)
    known_photos=set(previously.photo_id)
    known_observations=set(prev.observation_id)|set(s6.observation_id)|set(s7.observation_id)
    photo_collision=m.photo_id.isin(known_photos)
    known_uuid=m.observation_id.isin(known_observations)
    require(not known_uuid.any(),"Signed tier08 master repeated newer source observation UUID")
    novel=m.loc[~photo_collision].copy()
    require(novel.photo_id.is_unique and novel.observation_id.is_unique and len(novel)>=COUNT,
            "Insufficient exact ORIGINAL photo IDs after historical and tier06/07/910 removals")
    audit={
       "source_original_tier08_envmatched_island_photo_UUIDs":len(m),
       "already_measured_exact_original_photo_ids_across_source_arms":int(photo_collision.sum()),
       "tier08_eligible_exact_original_photo_ids_after_source_dedup":len(novel),
       "signed_canonical_historical_photo_ids":len(hist),
       "signed_tier06_photo_ids":len(s6),
       "signed_tier07_photo_ids":len(s7),
       "source_previous_tier08_910_photo_ids":len(prev),
       "prior_source_observation_UUID_repeats":int(known_uuid.sum()),
       "historical_observation_ID_formats_not_interchangeable_with_new_UUIDs":True,
    }
    return novel,audit

def frozen_area_light_queue(data, n=COUNT):
    """Island/site round robin, two matched exposure-extreme site-years first.

    All choices are based on precolour photo metadata; no pixel/ROI outcomes.
    Site arrays are grouped by island area terciles, then iterated over islands
    to avoid a few high-activity giant islands dominating the first wave.
    """
    require(len(data)>=n and data.photo_id.is_unique and data.observation_id.is_unique,
            "Cannot choose 24k unique original photos from tier08 source")
    d=data.copy()
    d["log_light"]=pd.to_numeric(d.log_light,errors="raise")
    d["area_km2"]=pd.to_numeric(d.area_km2,errors="raise")
    require(d.area_km2.gt(0).all() and d.log_light.notna().all(),
            "Not source-original same-year VIIRS or island area")
    stats=d.groupby(["land_id","year","target_id"],observed=True,sort=False).agg(
         attempts=("photo_id","size"),
         taxa=("wcvp_accepted_plant_name_id","nunique"),
         log_light=("log_light","first"),
         area_km2=("area_km2","first")).reset_index()
    stats["qualified_target"]=stats.attempts.ge(5)&stats.taxa.ge(3)
    qualifying=stats.loc[stats.qualified_target].copy()
    candidates={}
    for (land,year),g in qualifying.groupby(["land_id","year"],observed=True):
        if len(g)<2:continue
        sorted_light=g.sort_values(["log_light","target_id"],kind="stable")
        low,high=sorted_light.iloc[0],sorted_light.iloc[-1]
        if float(high.log_light-low.log_light)>=0.05 and low.target_id!=high.target_id:
            candidates[str(low.target_id)]="paired_dark_site"
            candidates[str(high.target_id)]="paired_bright_site"
    sites=stats.loc[stats.attempts.ge(3)&stats.taxa.ge(2)].copy()
    sites["paired"]=sites.target_id.astype(str).isin(candidates)
    sites["priority"]=sites.apply(
        lambda x: rank("site",x.target_id),axis=1)
    sites=sites.sort_values(
        ["qualified_target","paired","taxa","attempts","priority"],
        ascending=[False,False,False,False,True],kind="stable")
    areas=sites.groupby("land_id").area_km2.first().sort_values()
    classes=pd.cut(areas.rank(method="first",pct=True),
        bins=[0,1/3,2/3,1],labels=["small","medium","large"],include_lowest=True)
    class_map={str(k):str(v) for k,v in classes.items()}
    by_site={}
    for target,group in d.groupby("target_id",observed=True,sort=False):
        # Source accepted taxa get a first chance before repeated photos.
        group=group.copy()
        group["_h"]=group.photo_id.map(lambda p:rank("photo",p))
        group=group.sort_values("_h",kind="stable")
        first=group.drop_duplicates("wcvp_accepted_plant_name_id",keep="first")
        leftover=group.loc[~group.photo_id.isin(first.photo_id)]
        q=pd.concat([first,leftover],ignore_index=True).head(16)
        q["acquisition_arm"]=candidates.get(str(target),"area_light_site_densification")
        q["area_third"]=class_map.get(str(q.land_id.iloc[0]),"unknown")
        by_site[str(target)]=q.drop(columns="_h")
    # Round robin by area class AND island to maximize transfer, not volume
    # concentrated in a single large island.
    island_targets={}
    for island,group in sites.groupby("land_id",sort=False):
        island_targets[str(island)]=group.target_id.astype(str).tolist()
    classified={name:sorted(
        [land for land in island_targets if class_map.get(land)==name],
        key=lambda v:rank("island",v))
        for name in ("small","medium","large")}
    chosen=[]
    added=set()
    maxdepth=max(map(len,island_targets.values()))
    for depth in range(maxdepth):
        for cat in ("small","medium","large"):
            for land in classified[cat]:
                ids=island_targets[land]
                if depth>=len(ids):continue
                unit=by_site[ids[depth]]
                for row in unit.to_dict("records"):
                    if row["photo_id"] not in added:
                        chosen.append(row);added.add(row["photo_id"])
                        if len(chosen)>=n:break
                if len(chosen)>=n:break
            if len(chosen)>=n:break
        if len(chosen)>=n:break
    if len(chosen)<n:
        # A true original photo with fewer than 3 other opportunities can
        # still complement old source at its exact same location/year.
        remains=d.loc[~d.photo_id.isin(added)].copy()
        remains["_h"]=remains.photo_id.map(lambda p:rank("fallback",p))
        remains=remains.sort_values("_h",kind="stable")
        for row in remains.head(n-len(chosen)).to_dict("records"):
            row["acquisition_arm"]="hash_random_source_opportunity_fallback"
            row["area_third"]=class_map.get(str(row["land_id"]),"unknown")
            row.pop("_h",None)
            chosen.append(row)
    q=pd.DataFrame(chosen)
    require(len(q)==n and q.photo_id.is_unique and q.observation_id.is_unique,
            "Deterministic original photo queue shortage or repeated source UUID")
    q["measurement_id"]=q.photo_id.map(lambda p:"alan_tier08_bulk_"+p)
    q["image_filename"]=q.measurement_id+".jpg"
    q["source_snapshot_date"]="2026-09-27"
    q["shard"]=(np.arange(n)//PHOTO_UNIT).astype(int)
    q["wave"]=(np.arange(n)//WAVES[0]+1).astype(int)
    require(q.shard.nunique()==n//PHOTO_UNIT and
            q.groupby("shard").size().eq(PHOTO_UNIT).all() and
            q.wave.value_counts().to_dict()=={1:8192,2:8192,3:8192},
            "Frozen photo 128-photo shard or three-wave boundaries changed")
    audit={"size":len(q),"selected_distinct_original_GSHHG_islands":int(q.land_id.nunique()),
           "selected_original_source_site_years":int(q.target_id.nunique()),
           "selected_source_years":q.year.value_counts().to_dict(),
           "selected_area_third_photos":q.area_third.value_counts().to_dict(),
           "selected_matched_light_priority_photos":int(q.acquisition_arm.str.startswith("paired_").sum()),
           "selected_taxa_source_WCVP_IDs":int(q.wcvp_accepted_plant_name_id.nunique()),
           "precolour_candidate_island_year_light_extreme_sites":len(candidates),
           "source_selection_did_not_use_flower_ROI_or_photo_pixels":True}
    return q,audit

def prepare(master,previous910,canonical,source06,source07,out):
    out=Path(out)
    require(not out.exists(),"No overwrite of source-locked new 24,576 original photos")
    novel,audit=load_verified(master,previous910,canonical,source06,source07)
    q,selected=frozen_area_light_queue(novel)
    (out/"shards").mkdir(parents=True)
    hashes={}
    for i,frame in q.groupby("shard",sort=True):
        p=out/"shards"/f"shard-{i:03d}.csv"
        frame.to_csv(p,index=False,lineterminator="\n")
        hashes[str(i)]=sha(p)
    receipt={
       "schema":"mois_full_source_observation_sampling_cohort_v1",
       "colour_measured":False,
       "source_group":"alan_tier08_additional24576_balanced_island_area_light",
       "source_arm":ARM,
       "new_photo_attempts_planned":COUNT,
       "work_units":COUNT//PHOTO_UNIT,
       "worker_ids":list(range(COUNT//PHOTO_UNIT)),
       "photo_shard_sha256":hashes,
       "shard_per_warm_group":4,
       "warm_groups":COUNT//WORK_GROUP_PHOTOS,
       "waves":list(WAVES),
       "previous_910_original_photo_input_sha256":PREVIOUS910_SHA,
       "tier08_master_source_SHA256":MASTER_SHA,
       "source_previous_cohort_SHA256":{"canonical":panel.PHOTO_SHA,
                                      "tier06":TIER06_SHA,"tier07":TIER07_SHA},
       "source_overlap_audit":audit,
       "selection_audit":selected,
       "original_photos_not_colours_selected_before_pixels":True,
       "old_historical_observation_ID_vs_new_UUID_crosswalk_NOT_proven":True,
       "no_native_flora_or_causal_ALAN_claim":True,
    }
    (out/"manifest.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    q.drop(columns=["measurement_id","image_filename","source_snapshot_date","shard"]).to_csv(
        out/"precolour_24576_frozen_queue.csv.gz",
        index=False,compression={"method":"gzip","mtime":0},lineterminator="\n")
    print("ALAN_TIER08_24576_PRECOLOUR_ORIGINAL_PHOTO_SOURCE_FROZEN",
          json.dumps({"old_source_overlap":audit,
                      "selected":selected,"source_manifest_sha256":sha(out/"manifest.json")},
                     sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("master","previous910","canonical","source06","source07","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    prepare(a.master,a.previous910,a.canonical,a.source06,a.source07,a.out)
if __name__=="__main__":main()
