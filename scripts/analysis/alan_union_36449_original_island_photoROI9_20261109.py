#!/usr/bin/env python3
"""Verified cross-tier union of 36,449 authentic independent new island photo ROI9s.

Uses ONLY already measured, SHA-signed original photo results from:
  tier06: 8,942 photos = 5 signed 256-photo outcomes + 15 signed bulk groups
  tier07: 27,507 photos = signed 1,024-photo pilot + 52 signed 371-512-photo tasks

The two original sources were selected WITHOUT photo-colour information from
different immutable full93m original-photo hash tiers. Exact observation UUID
and source photo ID union is required. Strict photo cell×year qualification is
computed from the *union* at the same GSHHG island × 0.05-degree cell × original
observation year, so no overlapping site-year is double-counted, and separately
subthreshold photo sites can become eligible when genuine measured records from
the two disjoint original sources complement one another.

Not a probability-corrected worldwide floral-colour prevalence estimate and
not an ALAN causal/ecological regression. All biological unknowns remain unknown.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import valid_nine
from scripts.acquisition import alan_tier07_remaining26483_original_rgb_roi9 as tier07

SHA06="e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd"
SHA07="7198455819f1ac260a6cbeacd9f4536c30d766595daa64dc5205c75e9d5a1c62"
ROI_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
FIRST06_SHA="7789c5aa1916e5b1ff915d217744a80b990795cb24e362cff20007dda96ee94b"
PILOT07_SHA="2988a8631cfb25bea57690a2a91db1e0a7f00074b8ca3ce7557aba8c09fda113"
EXPECTED={"tier06":(8942,3597,156),"tier07":(27507,10245,769)}
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
COMPONENTS=tuple(f"strict9_{x}" for x in PALETTE)
SCHEMA="ALAN_cross_tier06_tier07_36449_source_real_original_ROI9_union_v1"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_36449_CROSS_TIER_SOURCE_UNION: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def source_rows(path,expected_sha,tier,n):
    require(sha(path)==expected_sha,"Original "+tier+" source photo roster bytes drift")
    d=pd.read_csv(path,compression="gzip",dtype=str,keep_default_na=False)
    needed=("photo_id","observation_id","target_id","land_id",
            "year","wcvp_accepted_plant_name_id")
    require(set(needed)<=set(d.columns),"Source missing original accepted taxon/geo keys")
    require(len(d)==n and d.photo_id.is_unique and d.observation_id.is_unique,
            "Source photo/observation UUIDs not unique")
    require(d.wcvp_accepted_plant_name_id.ne("").all() and
            d.target_id.eq(d.year+"|"+d.target_id.str.split("|",regex=False).str[1]+
                           "|"+d.land_id).all(),
            "Original source WCVP/real GSHHG cell-year construction drift")
    d=d.loc[:,list(needed)].copy()
    d["original_source_tier"]=tier
    return d

def read_numeric_photo(filename,expected_sha=None):
    filename=Path(filename)
    if expected_sha is not None:
        require(sha(filename)==expected_sha,"Measured original photo ROI output SHA changed")
    d=pd.read_csv(filename,compression="gzip" if filename.suffix==".gz" else None,
                  dtype=str,keep_default_na=False)
    require(d.photo_id.is_unique and d.observation_id.is_unique,
            "Duplicated original photo/observation UUID in measured source")
    if "strict9" in d.columns:
        good=d.strict9.str.lower().eq("true")
        require(set(COMPONENTS)<=set(d.columns),"Actual source strict9 numeric components absent")
        nums=d.loc[good,list(COMPONENTS)].apply(pd.to_numeric,errors="coerce")
        require(nums.notna().all().all() and nums.ge(0).all().all() and
                np.allclose(nums.sum(axis=1),1,atol=1e-5,rtol=0),
                "True numeric original nine-colour palette not compositional")
        require(d.loc[~good,list(COMPONENTS)].replace("",np.nan).isna().all().all(),
                "Technical ROI9 unknown was given a colour")
    else:
        require(all("flower_fraction_"+c in d.columns for c in PALETTE),
                "Original source FCP true nine numeric flower-region fractions missing")
        ok,_=valid_nine(d)
        good=pd.Series(ok,index=d.index)
    return d[["photo_id","observation_id"]].assign(real_numeric_ROI9=good.to_numpy(bool))

def read_prior06(directory):
    rows=[];receipts=[]
    for batch in range(1,6):
        if batch==1:
            root=Path(directory)/"first"
            file=root/"finished"/"new256_independent_island_actual_original_ROI9.csv.gz"
            receipt=root/"finished"/"ALAN_new256_original_photo_ROI9_source_receipt.json"
            key="fixed_FCP_original_method";expected=FIRST06_SHA
            success_key="original_nine_fraction_strict_colour_success"
        else:
            root=Path(directory)/f"alan-original-novel-island-ROI9-batch-{batch}"
            file=root/"finished"/f"ALAN_original_novel_island_photo_ROI9_batch_{batch:02d}.csv.gz"
            receipt=root/"finished"/"ALAN_original_novel_island_ROI9_batch_receipt.json"
            key="original_frozen_model_method";expected=None
            success_key="numeric_original_nine_palette_success"
        r=json.loads(receipt.read_text())
        if batch>1: expected=r["measured_original_colours_sha256"]
        require(r[key]["measurement_git_blob"]==ROI_BLOB,
                "Original five-group worker FCP source method changed")
        d=read_numeric_photo(file,expected)
        require(len(d)==256 and int(d.real_numeric_ROI9.sum())==r[success_key],
                "Previously signed first-five original photo outcomes changed")
        rows.append(d)
        receipts.append({"tier":"tier06","batch":batch,"photo_count":len(d),
                         "numeric_roi9":int(d.real_numeric_ROI9.sum()),
                         "photo_file_SHA256":sha(file)})
    return rows,receipts

def read_bulk06(directory):
    rows=[];receipts=[]
    for group in range(15):
        root=Path(directory)/f"alan-real-originalROI9-7662-group-{group}"
        frame=[]
        for shard in range(4):
            folder=root/"measured"/f"shard_{shard:03d}"
            file=folder/"photos.csv"
            r=json.loads((folder/"receipt.json").read_text())
            require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
                    and r["shard"]==shard and
                    r["photos_sha256"]==sha(file) and
                    r["method"]["measurement_git_blob"]==ROI_BLOB and
                    r["method"]["persisted_image_pixels"] is False,
                    "Original tier06 signed 128-photo worker outcome invalid")
            d=read_numeric_photo(file)
            require(len(d)==r["attempted"],"Original source photo shard incomplete")
            frame.append(d)
        q=pd.concat(frame,ignore_index=True)
        require(len(q)==(494 if group==14 else 512) and
                q.photo_id.is_unique and q.observation_id.is_unique,
                "Original completed tier06 bulk group not unique")
        rows.append(q)
        receipts.append({"tier":"tier06","warm_group":group,
                         "photo_count":len(q),
                         "numeric_roi9":int(q.real_numeric_ROI9.sum())})
    return rows,receipts

def read_pilot07(directory):
    root=Path(directory)
    file=root/"finished"/"ALAN_tier07_first1024_original_floral_roi9.csv.gz"
    receipt=root/"finished"/"ALAN_tier07_first1024_original_roi9_receipt.json"
    r=json.loads(receipt.read_text())
    require(r.get("schema")=="ALAN_tier07_1024_real_original_FCP_ROI9_pilot_v1"
            and r["original_frozen_FCP_ROI9_method"]["measurement_git_blob"]==ROI_BLOB
            and r["real_original_ROI9_photo_outcomes_sha256"]==PILOT07_SHA,
            "Pilot tier07 original 1024 photo source method/receipt changed")
    d=read_numeric_photo(file,PILOT07_SHA)
    require(len(d)==1024 and int(d.real_numeric_ROI9.sum())==403,
            "Original tier07 1024 observed numeric ROI source drift")
    return d,{"tier":"tier07","part":"pilot1024","photo_count":len(d),
              "numeric_roi9":403,"original_photo_SHA256":PILOT07_SHA}

def read_bulk07(directory):
    rows=[];receipts=[]
    for task in range(52):
        root=Path(directory)/f"alan-tier07-real-originalROI9-26483-task-{task}"
        file=root/"finished"/f"alan_tier07_original_photo_ROI9_task{task:02d}.csv.gz"
        rec=root/"finished"/"ALAN_tier07_genuine_new_photo_ROI9_task_receipt.json"
        r=json.loads(rec.read_text())
        need=371 if task==51 else 512
        require(r.get("schema")=="ALAN_tier07_26483_original_photo_FCP_ROI9_task_v1" and
                r["task_number"]==task and
                r["original_photo_image_attempts"]==need and
                r["original_FCP_model_sha_blob"]==ROI_BLOB and
                r["source_schedule_SHA256"]==SHA07 and
                r["original_RGB_image_bytes_not_persisted"] is True,
                "Original tier07 FCP photo worker provenance drift")
        d=read_numeric_photo(file,r["original_real_photo_colour_file_SHA256"])
        require(len(d)==need and
                int(d.real_numeric_ROI9.sum())==r["actual_numerical_9palette_flowerROI_success"],
                "Actual original tier07 numeric source photo count changed")
        rows.append(d)
        receipts.append({"tier":"tier07","task":task,"photo_count":len(d),
                         "numeric_roi9":int(d.real_numeric_ROI9.sum()),
                         "photo_file_SHA256":sha(file)})
    return rows,receipts

def site_stats(data):
    require(data.photo_id.is_unique and data.observation_id.is_unique,
            "Source ROI9 observation UUIDs reused")
    group=data.groupby("target_id",observed=True).agg(
       photo_opportunities=("observation_id","nunique"),
       accepted_source_species=("wcvp_accepted_plant_name_id","nunique"),
       actual_numerical_ROI9=("real_numeric_ROI9","sum"),
       island=("land_id","first"),
       year=("year","first"))
    measured=data.loc[data.real_numeric_ROI9]
    species=measured.groupby("target_id",observed=True).wcvp_accepted_plant_name_id.nunique()
    group["successfully_measured_accepted_species"]=species.reindex(group.index,fill_value=0)
    group["source_opportunity_5photo_3taxa"]=(
       group.photo_opportunities.ge(5)&group.accepted_source_species.ge(3))
    group["actual_strict_5photo_3success_3taxa"]=(
       group.source_opportunity_5photo_3taxa&
       group.actual_numerical_ROI9.ge(3)&
       group.successfully_measured_accepted_species.ge(3))
    return group

def summary(data,s):
    return {"independent_original_observations":len(data),
            "original_9palette_numeric_ROI9":int(data.real_numeric_ROI9.sum()),
            "source_islands":int(data.land_id.nunique()),
            "source_cell_years":len(s),
            "source_5attempt_3taxa_opportunity_cells":int(s.source_opportunity_5photo_3taxa.sum()),
            "actual_strict_measured_cells":int(s.actual_strict_5photo_3success_3taxa.sum()),
            "source_islands_with_actual_strict_cells":int(s.loc[
               s.actual_strict_5photo_3success_3taxa,"island"].nunique())}

def run(source06,source07,prior06,bulk06,pilot07,bulk07,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwriting previously signed source cell union result")
    s6=source_rows(source06,SHA06,"tier06",8942)
    s7=source_rows(source07,SHA07,"tier07",27507)
    require(set(s6.observation_id).isdisjoint(set(s7.observation_id)) and
            set(s6.photo_id).isdisjoint(set(s7.photo_id)),
            "Two original full93m source hash tiers share photo/observation IDs")
    first,proof_first=read_prior06(prior06)
    bulk6,proof_bulk6=read_bulk06(bulk06)
    pil,proof_pil=read_pilot07(pilot07)
    bulk7,proof_bulk7=read_bulk07(bulk07)
    measured=pd.concat([*first,*bulk6,pil,*bulk7],ignore_index=True)
    expected=pd.concat([s6,s7],ignore_index=True)
    require(len(measured)==len(expected)==36449 and measured.photo_id.is_unique and
            measured.observation_id.is_unique,
            "Original 36,449 FCP image source measurements not one-to-one")
    d=expected.merge(measured,on=["photo_id","observation_id"],
                     how="inner",validate="one_to_one")
    require(len(d)==36449 and set(d.photo_id)==set(expected.photo_id) and
            set(d.observation_id)==set(expected.observation_id),
            "Measured actual photo identities differ from precolour GSHHG source")
    tier6=d.loc[d.original_source_tier.eq("tier06")].copy()
    tier7=d.loc[d.original_source_tier.eq("tier07")].copy()
    t6=site_stats(tier6);t7=site_stats(tier7);both=site_stats(d)
    require(len(tier6)==EXPECTED["tier06"][0] and
            int(tier6.real_numeric_ROI9.sum())==EXPECTED["tier06"][1] and
            int(t6.actual_strict_5photo_3success_3taxa.sum())==EXPECTED["tier06"][2],
            "Original independently verified tier06 photo/scientific strict-site count drift")
    require(len(tier7)==EXPECTED["tier07"][0] and
            int(tier7.real_numeric_ROI9.sum())==EXPECTED["tier07"][1] and
            int(t7.actual_strict_5photo_3success_3taxa.sum())==EXPECTED["tier07"][2],
            "Original independently verified tier07 photo/scientific strict-site count drift")
    u6=set(t6.index[t6.actual_strict_5photo_3success_3taxa])
    u7=set(t7.index[t7.actual_strict_5photo_3success_3taxa])
    joined=set(both.index[both.actual_strict_5photo_3success_3taxa])
    preexisting=u6|u7
    require(preexisting.issubset(joined) and
            len(joined)>=max(len(u6),len(u7)),
            "True 5attempt/3numericROI/3accepted photo source union lost measured sites")
    byyear=both.groupby("year",observed=True).agg(
       photo_cell_years=("year","size"),
       source_5photo_3taxa_opportunities=("source_opportunity_5photo_3taxa","sum"),
       actual_strict_ROI9_cells=("actual_strict_5photo_3success_3taxa","sum")
    ).reset_index()
    state={
       "schema":SCHEMA,
       "tier06_source_roster_SHA256":SHA06,
       "tier07_source_roster_SHA256":SHA07,
       "original_FCP_ROI9_method_git_blob":ROI_BLOB,
       "tier06_source_measured":summary(tier6,t6),
       "tier07_source_measured":summary(tier7,t7),
       "combined_same_site_year_original_source":summary(d,both),
       "tier06_and_tier07_original_observation_UUID_and_photo_ID_sets_disjoint":True,
       "strict_cells_qualified_in_both_tiers_separately":len(u6&u7),
       "strict_cells_qualified_in_at_least_one_tier_before_union":len(preexisting),
       "additional_strict_sites_created_by_real_photos_across_tiers":len(joined-preexisting),
       "source_sites_photographed_in_both_tiers":int(len(set(t6.index)&set(t7.index))),
       "unique_strict_sites_in_original_new_observation_source_after_union":len(joined),
       "all_photo_source_receipts_joined":len(proof_first)+len(proof_bulk6)+1+len(proof_bulk7),
       "photo_ROI_technical_unknowns_never_assigned_a_hue":True,
       "original_RBG_pixels_not_downloaded_or_reinferred_in_union":True,
       "island_and_site_targeting_not_worldwide_random":True,
       "verified_native_bloom_or_ALAN_causality":False,
       "historical_main_fiveyear_ALAN_regression_unchanged":True
    }
    require(state["all_photo_source_receipts_joined"]==73 and
            state["combined_same_site_year_original_source"]["independent_original_observations"]==36449 and
            state["combined_same_site_year_original_source"]["original_9palette_numeric_ROI9"]==13842,
            "Original source 73 actual-photo worker identity or outcome counts changed")
    out.mkdir(parents=True)
    p=out/"ALAN_36449_original_photoROI9_strict_source_years.csv"
    byyear.to_csv(p,index=False)
    state["year_counts_sha256"]=sha(p)
    (out/"ALAN_36449_tier06_tier07_original_ROI9_union_receipt.json").write_text(
        json.dumps(state,sort_keys=True,indent=2)+"\n")
    print("ALAN_36449_ORIGINAL_ROI9_CROSS_TIER_REAL_SITE_UNION",
          json.dumps(state,sort_keys=True),flush=True)
    return state

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("source06","source07","prior06","bulk06","pilot07","bulk07","out"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.source06,a.source07,a.prior06,a.bulk06,a.pilot07,a.bulk07,a.out)
