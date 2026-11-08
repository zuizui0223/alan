#!/usr/bin/env python3
"""ALAN original-source photo-ROI opportunities: site-year densification, not hues.

From the already signed 8,942 *new* observation UUID photo source, preserve
batches 1-5 (1,280 original photo IDs) and freeze 512 further UNMEASURED source
photo UUIDs with >=5 genuinely available photo observations and >=3 original
WCVP accepted taxa in their corresponding 0.05° GSHHG landmass cell × year.
Prioritize islands with no earlier strict-qualified ALAN photographic cell,
then source-rich site-year units. Do not use the new images' unknown/future ROI
outcomes or invent white-colour labels. This is a *separate, adaptive source
design*, not an unbiased estimate of world/island plant-colour prevalence.
"""
from __future__ import annotations

import argparse,hashlib,json
from collections import Counter
from pathlib import Path

import pandas as pd

SCHEDULE_SHA="e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd"
SOURCE_SCHEMA="ALAN_8942_new_original_island_UUIDs_progressive_256photo_source_v1"
OUT_SCHEMA="ALAN_new_original_island_site_densification_512_source_v1"
PREVIOUS_COMPLETED_OR_EXECUTING_BATCHES=5
SOURCE_N=8942
TARGET_N=512
PER_ISLAND_MAX=64
PER_SITE_MAX=16
SEED="ALAN_SITE_DENSIFICATION_BEFORE_COLOUR_20261108_v1"
SOURCE_ARM="alan93m_site_year_densification_original_uuid_source_v1"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_ISLAND_SITE_SOURCE_512: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def rank(tag,x):
    return hashlib.sha256((SEED+"|"+tag+"|"+str(x)).encode()).hexdigest()

def candidate_cell_support(frame):
    d=frame.copy()
    d["source_batch"]=d.source_batch.astype(int)
    return d.groupby(["target_id","land_id"],observed=True).agg(
       available_original_photo_observations=("observation_id","nunique"),
       accepted_original_taxa=("wcvp_accepted_plant_name_id","nunique"),
       observed_in_original_first_five=("source_batch",lambda s:int((s<=5).sum())),
       historical_strict_island_cells=("existing_strict_island_cell_years","first")
    ).reset_index()

def eligible_source_sites(frame):
    d=candidate_cell_support(frame)
    x=d.loc[(d.available_original_photo_observations>=5)&
            (d.accepted_original_taxa>=3)&
            (d.available_original_photo_observations>d.observed_in_original_first_five)].copy()
    x["zero_strict_prior_island"]=x.historical_strict_island_cells.astype(int).eq(0)
    x["source_rank"]=x.target_id.map(lambda k:rank("site",k))
    return x.sort_values(
       ["zero_strict_prior_island","available_original_photo_observations",
        "accepted_original_taxa","source_rank"],
       ascending=[False,False,False,True],kind="stable")

def only_real_source_row_identity(frame):
    require(len(frame)==SOURCE_N and frame.photo_id.is_unique and
            frame.observation_id.is_unique and
            frame.land_id.nunique()==369 and
            frame.source_batch.astype(int).between(1,35).all(),
            "Source original 8,942 accepted photographed UUIDs changed")
    require(frame.colour_measured.astype(str).str.lower().eq("false").all() and
            frame.selection_probability.fillna("").eq("").all(),
            "Selected source colours or inclusion probabilities were fabricated")
    require(frame.photo_id.astype(str).str.fullmatch(r"\d+").all() and
            frame.photo_url_large.str.fullmatch(
              r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
            ).all(),"Original photo URLs or numeric IDs not valid source")
    require(frame.wcvp_accepted_plant_name_id.notna().all() and
            frame.wcvp_accepted_plant_name_id.ne("").all(),
            "Original source photo has no actually accepted WCVP taxon")
    return True

def select_source(frame,limit=TARGET_N):
    """Differently selected SOURCE rows; never repeat original 1-5 pixels."""
    require(1<=limit<=TARGET_N,"Not an original source-bounded site-wave")
    only_real_source_row_identity(frame)
    batch=frame.source_batch.astype(int)
    first=frame.loc[batch<=PREVIOUS_COMPLETED_OR_EXECUTING_BATCHES].copy()
    rest=frame.loc[batch>PREVIOUS_COMPLETED_OR_EXECUTING_BATCHES].copy()
    require(len(first)==1280 and len(rest)==7662,
            "Previously measured / running first 1280 observation UUIDs changed")
    valid=eligible_source_sites(frame)
    source_count=len(valid)
    picked=[];obs=set();island_counts=Counter();n=0
    # Original same photo-cell-year already shares same genuine original
    # year and precomputed VIIRS/WorldPop; an island code alone is not a site.
    rest_by_cell={k:v for k,v in rest.groupby("target_id",sort=False)}
    for site in valid.itertuples(index=False):
        if n>=limit:break
        land=str(site.land_id)
        if island_counts[land]>=PER_ISLAND_MAX:continue
        group=rest_by_cell.get(site.target_id)
        if group is None:continue
        q=group.copy()
        q["_source_rank"]=q.observation_id.map(lambda v:rank("photo",v))
        q=q.sort_values(["_source_rank","observation_id"],kind="stable")
        one_per_taxon=q.drop_duplicates("wcvp_accepted_plant_name_id",keep="first")
        other=q.loc[~q.observation_id.isin(set(one_per_taxon.observation_id))]
        prioritized=pd.concat([one_per_taxon,other],ignore_index=True)
        take=min(limit-n,PER_SITE_MAX,PER_ISLAND_MAX-island_counts[land],len(prioritized))
        if take<2:continue
        chosen=prioritized.head(take).drop(columns="_source_rank").copy()
        require(chosen.observation_id.is_unique and
                not set(chosen.observation_id)&obs,
                "Repeated original observation UUID in separate source sites")
        picked.append(chosen)
        obs.update(chosen.observation_id)
        island_counts[land]+=take
        n+=take
    require(n==limit,"Not enough verified independent existing island source photos")
    selected=pd.concat(picked,ignore_index=True)
    require(selected.photo_id.is_unique and selected.observation_id.is_unique
            and not set(selected.photo_id)&set(first.photo_id)
            and not set(selected.observation_id)&set(first.observation_id),
            "New site-wave source contains photos already measured in batches 1-5")
    selected["site_wave_photo_order"]=range(len(selected))
    selected["source_arm"]=SOURCE_ARM
    selected["selection_probability"]=""
    selected["site_wave_source_batch"]=(selected.site_wave_photo_order//256)+1
    selected["source_snapshot_date"]="2026-09-27"
    selected["source_photo_pixels_measured"]=False
    return selected,{
       "original_photos_in_first_five_batches_excluded":len(first),
       "source_sites_with_5plus_observations_and_3plus_WCVP_accepted_taxa":source_count,
       "photo_observation_candidates_after_prior_five_source_batches":len(rest),
    }

def number_of_potential_strict_source_sites(source):
    g=source.groupby("target_id",observed=True).agg(
       original_photo_attempts=("observation_id","nunique"),
       original_accepted_taxa=("wcvp_accepted_plant_name_id","nunique"))
    return int((g.original_photo_attempts.ge(5)&
                g.original_accepted_taxa.ge(3)).sum())

def run(schedule,receipt,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite a signed source-photo acquisition roster")
    require(sha(schedule)==SCHEDULE_SHA,"Original 8,942 photo source SHA changed")
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")==SOURCE_SCHEMA and
            r.get("source_schedule_csv_sha256")==SCHEDULE_SHA and
            r.get("novel_original_observation_UUIDs")==SOURCE_N and
            r.get("source_batches")==35 and
            r.get("first_batch_unchanged") is True and
            r.get("no_photo_images_downloaded_or_hue_ROI_inferred") is True,
            "Not the publicly signed source-only, novel original island photo queue")
    data=pd.read_csv(schedule,compression="gzip",dtype=str,keep_default_na=False)
    selected,meta=select_source(data)
    first5=data.loc[data.source_batch.astype(int).le(5)].copy()
    broad7=data.loc[data.source_batch.astype(int).le(7)].copy()
    densified=pd.concat([first5,selected],ignore_index=True)
    candidate={
        "prior_original_five_batches":number_of_potential_strict_source_sites(first5),
        "original_breadth_first_seven_batches":number_of_potential_strict_source_sites(broad7),
        "with_512_new_site_targeted_photographs":number_of_potential_strict_source_sites(densified),
    }
    require(candidate["prior_original_five_batches"]==3 and
            candidate["original_breadth_first_seven_batches"]==6 and
            candidate["with_512_new_site_targeted_photographs"]==60,
            "Previously independently audited source-only 5-photo/3-accepted site potential changed")
    require(int(selected.land_id.nunique())==20 and
            int(selected.target_id.nunique())==61,
            "Original site densification site/island coverage unexpectedly changed")
    out.mkdir(parents=True)
    file=out/"alan_site_year_densification512_original_photo_source.csv.gz"
    selected.to_csv(file,index=False,compression={"method":"gzip","mtime":0},
                   lineterminator="\n")
    cells=selected.groupby(["site_wave_source_batch"],observed=True).agg(
        new_photo_observations=("observation_id","nunique"),
        covered_source_photo_cells=("target_id","nunique"),
        covered_GSHHG_islands=("land_id","nunique")).reset_index()
    cells.to_csv(out/"site_dense512_batch_coverage_aggregate.csv",index=False)
    report={
        "schema":OUT_SCHEMA,
        "precolour_original_full8942_source_SHA256":SCHEDULE_SHA,
        "preselected_source_first5_original_1280_not_reused":True,
        "new_original_photo_observation_UUIDs_in_source":len(selected),
        "new_actual_source_GSHHG_islands":int(selected.land_id.nunique()),
        "new_actual_source_photo_cell_years":int(selected.target_id.nunique()),
        "two_disjoint_256photo_source_work_waves":int(selected.site_wave_source_batch.nunique()),
        "source_candidate_cell_thresholds_not_successful_photo_ROI_results":True,
        "potential_5plus_attempt_3plus_accepted_site_counts":candidate,
        "source_selection_requirements":meta,
        "selected_source_queue_SHA256":sha(file),
        "batch_coverage_csv_SHA256":sha(out/"site_dense512_batch_coverage_aggregate.csv"),
        "no_photo_pixels_obtained_or_colour_relabelled":True,
        "not_confirmed_open_native_flower_or_nightlight_effect":True,
        "cannot_pool_targeted_island_oversample_as_unbiased_global_colour_distribution":True,
        "color_photo_measurements_from_previous_five_batches_were_NOT_used_to_select_targets":True,
    }
    (out/"alan_512_site_densification_original_source_receipt.json").write_text(
         json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_AUTHENTIC_SITE_DENSIFICATION_SOURCE_512",
          json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schedule",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.schedule,a.receipt,a.out)
