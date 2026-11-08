#!/usr/bin/env python3
"""ALAN source-blind original tier07 27,507 photo/observation UUID schedule.

Input: the independently SHA-verified, original-year VIIRS/WorldPop-joined
island photographic source from 1,193,829-photo full93m tier07.
No original photos are downloaded, opened, or inferred by this plan.

Group0 has 256 observation UUIDs on 256 previously strict-unsampled islands
and 768 further distinct observation UUIDs across photo-rich site×years,
where 5 original source photos and 3 genuinely observed WCVP accepted taxa
could support a strict cell after actual ROI9 measurement. Max 12 photos/site
and 48 per island in the densification subset. All remaining 26,483 original
UUIDs are scheduled without replacement into groups 1–26 (1024 photos each,
except 883 in the final group). This is source-conditional, not a globally
representative island-flower sample, and does NOT produce any colour label.
"""
from __future__ import annotations

import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import pandas as pd
import numpy as np

MASTER_SHA="e38161652d3ec85928fa0f4fad91d70c123a7cb2af32f65f199f42504f76e056"
MASTER_N=27507
LANDMASSES=563
SITE_YEARS=9493
GROUP_N=27
GROUP_SIZE=1024
SOURCE_ARM="alan_93m_tier07_new_island_observation_ROI_source_v1"
PLAN_ARM="alan_93m_tier07_source_27k_reference_plus_site_densification_v1"
PLAN_SCHEMA="ALAN_tier07_27507_new_original_island_1024photo_source_schedule_v1"
SALT="ALAN_TIER07_WARM1024_PRECOLOUR_20261108"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_TIER07_27507_SOURCE: "+msg)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def rank(tag,v):
    return hashlib.sha256((SALT+"|"+tag+"|"+str(v)).encode()).hexdigest()

def validated_source(path):
    require(sha(path)==MASTER_SHA,"Authentic tier07 original-photo source SHA changed")
    d=pd.read_csv(path,compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==MASTER_N and d.photo_id.is_unique and
            d.observation_id.is_unique and d.land_id.nunique()==LANDMASSES and
            d.target_id.nunique()==SITE_YEARS,
            "Original 27507 distinct photo/observation/source-site counts drift")
    require(d.source_arm.eq(SOURCE_ARM).all() and
            d.selection_probability.eq("").all() and
            d.photo_id.str.fullmatch(r"\d+").all() and
            d.photo_url_large.str.fullmatch(
              r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
            ).all(),
            "Original source URLs, acceptance or conditional design changed")
    require(d.wcvp_accepted_plant_name_id.ne("").all() and
            d.year.isin(("2019","2021","2023","2024","2025")).all(),
            "Original WCVP accepted taxa/real photo study years absent")
    return d

def site_opportunity_count(frame):
    g=frame.groupby("target_id",observed=True).agg(
       attempts=("observation_id","nunique"),
       species=("wcvp_accepted_plant_name_id","nunique"))
    return int((g.attempts.ge(5)&g.species.ge(3)).sum())

def choose_first1024(original):
    d=original.copy()
    d["_rank"]=d.observation_id.map(lambda x:rank("photo",x))
    d["_oldstrict"]=pd.to_numeric(
        d.existing_strict_island_cell_years,errors="raise").astype(int)
    # Coverage arm: one genuine photo from 256 independent source islands
    # with no previous strict colour cell, selected blind to future pixels.
    across=(d.sort_values(["_oldstrict","_rank"],kind="stable")
              .drop_duplicates("land_id")
              .sort_values(["_oldstrict","_rank"],kind="stable").head(256))
    require(len(across)==256 and across.land_id.nunique()==256 and
            across._oldstrict.eq(0).all(),
            "Original precolour broad source island opportunity has drifted")
    across=across.copy();across["source_priority"]="coverage_reference"
    # Density arm: sites with >=8 existing source photos and >=3 accepted
    # taxa; choose distinct taxa before further observations of same taxon.
    rest=d.loc[~d.observation_id.isin(set(across.observation_id))].copy()
    sites=d.groupby("target_id",observed=True).agg(
         n=("observation_id","nunique"),
         taxa=("wcvp_accepted_plant_name_id","nunique"),
         baseline=("_oldstrict","first"),
         island=("land_id","first")).reset_index()
    sites=sites.loc[sites.n.ge(8)&sites.taxa.ge(3)].copy()
    sites["_site_rank"]=sites.target_id.map(lambda t:rank("site",t))
    sites=sites.sort_values(["baseline","n","taxa","_site_rank"],
                            ascending=[True,False,False,True],kind="stable")
    island_quota=Counter()
    take=[];n=0
    pool={k:v for k,v in rest.groupby("target_id",sort=False)}
    for s in sites.itertuples(index=False):
        if n>=768:break
        land=str(s.island)
        if island_quota[land]>=48:continue
        rows=pool.get(s.target_id)
        if rows is None:continue
        ordered=rows.sort_values(["_rank","observation_id"],kind="stable")
        first=ordered.drop_duplicates("wcvp_accepted_plant_name_id",keep="first")
        more=ordered.loc[~ordered.observation_id.isin(set(first.observation_id))]
        prioritized=pd.concat([first,more],ignore_index=True)
        k=min(12,768-n,48-island_quota[land],len(prioritized))
        if k<5:continue
        chosen=prioritized.head(k).copy()
        take.append(chosen)
        island_quota[land]+=k
        n+=k
    require(n==768,"Original photo-year taxon-supported source no longer admits density arm")
    dense=pd.concat(take,ignore_index=True)
    dense["source_priority"]="siteyear_densification"
    result=pd.concat([across,dense],ignore_index=True)
    require(len(result)==1024 and result.photo_id.is_unique and
            result.observation_id.is_unique,
            "Genuinely new 1024 source photo IDs/UUIDs not unique")
    require(result.land_id.nunique()>=280 and
            site_opportunity_count(result)>=70,
            "Source-blind pilot broad coverage or density opportunity reduced")
    return result.drop(columns=["_rank","_oldstrict"]),{
       "reference_original_photos":256,
       "densification_original_photos":768,
       "reference_independent_islands":int(across.land_id.nunique()),
       "density_independent_islands":int(dense.land_id.nunique()),
       "combined_independent_islands":int(result.land_id.nunique()),
       "combined_distinct_photo_cell_years":int(result.target_id.nunique()),
       "combined_precolour_5photo_3accepted_taxon_sites":site_opportunity_count(result)
    }

def freeze(master,out):
    out=Path(out)
    require(not out.exists(),"Original tier07 photo UUID schedule must never be overwritten")
    d=validated_source(master)
    first,pilot=choose_first1024(d)
    # The remaining photo IDs still enter the final acquisition, regardless
    # of inferred colour or future ROI failure.
    unseen=d.loc[~d.observation_id.isin(set(first.observation_id))].copy()
    unseen["_rank"]=unseen.observation_id.map(lambda v:rank("remaining",v))
    site=d.groupby("target_id",observed=True).agg(
        n=("observation_id","nunique"),
        species=("wcvp_accepted_plant_name_id","nunique")).reset_index()
    site["strict_opportunity"]=site.n.ge(5)&site.species.ge(3)
    unseen=unseen.merge(site[["target_id","strict_opportunity","n","species"]],
                        on="target_id",validate="many_to_one")
    unseen=unseen.sort_values(
         ["strict_opportunity","n","species","_rank"],
         ascending=[False,False,False,True],kind="stable")
    unseen=unseen.drop(columns=["_rank","strict_opportunity","n","species"])
    unseen["source_priority"]="remaining_source_completion"
    scheduled=pd.concat([first,unseen],ignore_index=True)
    require(len(scheduled)==MASTER_N and
            scheduled.photo_id.is_unique and scheduled.observation_id.is_unique and
            set(scheduled.photo_id)==set(d.photo_id) and
            set(scheduled.observation_id)==set(d.observation_id),
            "Source scheduling altered original new independent observation cohort")
    scheduled["work_group"]=np.arange(MASTER_N)//GROUP_SIZE
    scheduled["within_work_group"]=np.arange(MASTER_N)%GROUP_SIZE
    scheduled["source_colour_measured"]=False
    require(scheduled.work_group.nunique()==GROUP_N and
            scheduled.groupby("work_group").size().iloc[:-1].eq(GROUP_SIZE).all()
            and int(scheduled.groupby("work_group").size().iloc[-1])==883,
            "Source 27507 must be 26x1024 plus 883 original-photo groups")
    out.mkdir(parents=True)
    path=out/"alan_tier07_all27507_new_island_original_photo_1024group_source.csv.gz"
    scheduled.to_csv(path,index=False,
                     compression={"method":"gzip","mtime":0},lineterminator="\n")
    g=scheduled.groupby("work_group",observed=True).agg(
       distinct_original_UUIDs=("observation_id","nunique"),
       distinct_islands=("land_id","nunique"),
       distinct_photo_site_years=("target_id","nunique")).reset_index()
    g.to_csv(out/"source_27group_geographic_photo_opportunities.csv",index=False)
    receipt={
      "schema":PLAN_SCHEMA,
      "original_27507_photo_source_SHA256":MASTER_SHA,
      "original_source_accepted_new_observation_UUIDs":MASTER_N,
      "original_source_distinct_GSHHG_islands":LANDMASSES,
      "original_photo_years":[2019,2021,2023,2024,2025],
      "original_source_photo_cell_years":SITE_YEARS,
      "disjoint_genuine_source_groups":GROUP_N,
      "full_groups_of_1024_original_images":26,
      "last_partial_source_group_size":883,
      "first_1024_prespecified_blind_source":pilot,
      "source_photo_selection_did_not_use_colour":True,
      "original_new_image_bytes_downloaded":0,
      "new_original_nine_photo_ROI9_inference_performed":0,
      "no_selection_probability_invented_for_targeted_source":True,
      "different_source_tier07_from_already_completed_8942_tier06":True,
      "precolour_source_photo_plan_SHA256":sha(path),
      "work_group_support_SHA256":sha(out/"source_27group_geographic_photo_opportunities.csv"),
    }
    (out/"ALAN_tier07_27507_original_photo_source_schedule_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("ALAN_TIER07_27507_ORIGINAL_PHOTO_SOURCE_SCHEDULE",
          json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    freeze(a.master,a.out)
