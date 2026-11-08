#!/usr/bin/env python3
"""ALAN signed aggregate of 1,280 genuinely new island source photos (5 batches).

Merge the already completed original 256 photographs and four new 256-photo
genuine FCP YOLO11/EfficientSAM ROI9 source workers on EXACT original photo &
observation UUIDs from the immutable 8,942 source. Do not assign technical
unknowns a colour, re-download images, or claim that a source-only 5photo/3spp
site is strict3-success measured. Keep original 2019-25 main model unchanged.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

import numpy as np
import pandas as pd

SCHEDULE_SHA="e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd"
FIRST_ACTUAL_SHA="7789c5aa1916e5b1ff915d217744a80b990795cb24e362cff20007dda96ee94b"
FCP_MODEL_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
SCHEMA="ALAN_real1280_new_island_originalROI9_five_batch_merge_v1"
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
CHANNELS=tuple("strict9_"+x for x in PALETTE)

def require(x,msg):
    if not bool(x):raise ValueError("ALAN_REAL1280_ORIGINAL_ROI9_MERGE: "+msg)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for v in iter(lambda:f.read(4<<20),b""):h.update(v)
    return h.hexdigest()

def source_schedule(schedule,receipt):
    require(sha(schedule)==SCHEDULE_SHA,"Photo-colour blind 8942 source batch schedule SHA changed")
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")=="ALAN_8942_new_original_island_UUIDs_progressive_256photo_source_v1"
            and r.get("source_schedule_csv_sha256")==SCHEDULE_SHA and
            r.get("first_batch_unchanged") is True,
            "Not authentic original ALAN prospective source")
    source=pd.read_csv(schedule,compression="gzip",dtype=str,keep_default_na=False)
    require(len(source)==8942 and source.photo_id.is_unique and
            source.observation_id.is_unique and
            source.source_batch.astype(int).between(1,35).all(),
            "Original source 8942 data are not unique")
    return source

def read_authenticated_batch(root,batch):
    root=Path(root)
    if batch==1:
        file=root/"finished"/"new256_independent_island_actual_original_ROI9.csv.gz"
        receipt=root/"finished"/"ALAN_new256_original_photo_ROI9_source_receipt.json"
        require(sha(file)==FIRST_ACTUAL_SHA,"Original first256 FCP ROI9 source outcome SHA drift")
        info=json.loads(receipt.read_text())
        require(info.get("schema")=="ALAN_new256_original_island_unique_UUID_real_FCP_ROI9_v1"
                and info["original_source_photo_attempts"]==256
                and info["independent_new_observation_UUIDs"]==256
                and info["original_nine_fraction_strict_colour_success"]==92
                and info["original_photo_colour_output_SHA256"]==FIRST_ACTUAL_SHA,
                "Original successfully measured first256 photo receipt changed")
        success=info["original_nine_fraction_strict_colour_success"]
        model=info["fixed_FCP_original_method"]
    else:
        file=root/"finished"/f"ALAN_original_novel_island_photo_ROI9_batch_{batch:02d}.csv.gz"
        receipt=root/"finished"/"ALAN_original_novel_island_ROI9_batch_receipt.json"
        info=json.loads(receipt.read_text())
        require(info.get("schema")=="ALAN_novel_source_island_original256_FCP_ROI9_batch_v1"
                and info["source_batch"]==batch
                and info["original_images_downloaded_and_measured"]==256
                and info["new_distinct_observation_UUIDs"]==256
                and info["measured_original_colours_sha256"]==sha(file)
                and info["technical_unknown_photo_not_assigned_an_ROI9_colour"],
                "New 256 original-photo RGB ROI9 source receipt missing or incorrect")
        success=info["numeric_original_nine_palette_success"]
        model=info["original_frozen_model_method"]
    require(model.get("measurement_git_blob")==FCP_MODEL_BLOB
            and model.get("persisted_image_pixels") is False,
            "The signed model and image-ephemeral procedure changed")
    d=pd.read_csv(file,compression="gzip",dtype=str,keep_default_na=False)
    require(len(d)==256 and d.photo_id.is_unique and
            d.observation_id.is_unique and set(CHANNELS)<=set(d.columns),
            "Real original ROI9 photo IDs or source numeric components invalid")
    d["batch_number"]=batch
    d["real_numeric_ROI9"]=d.strict9.str.lower().eq("true")
    require(int(d.real_numeric_ROI9.sum())==success,
            "Original independent photo nine-palette success totals changed")
    nine=d.loc[d.real_numeric_ROI9,list(CHANNELS)].apply(pd.to_numeric,errors="coerce")
    require(nine.notna().all().all() and
            np.allclose(nine.sum(axis=1),1,atol=1e-5,rtol=0),
            "Real measured photo nine-palette component sum invalid")
    unknown=d.loc[~d.real_numeric_ROI9,list(CHANNELS)]
    require(unknown.replace("",np.nan).isna().all().all(),
            "Technical-unknown source photo was reclassified as a colour")
    return d,{"batch":batch,"measured_photo_source_sha256":sha(file),
              "actual_original_ROI9_success":success,
              "actual_original_ROI9_unknown":256-success,
              "frozen_FCP_method_blob":FCP_MODEL_BLOB}

def actual_colour_site_support(frame):
    q=frame.copy()
    require(q.observation_id.is_unique and q.photo_id.is_unique,
            "Original new observations cannot be pseudo-replicated")
    g=q.groupby("target_id",observed=True).agg(
        original_photo_observations=("observation_id","nunique"),
        actual_numeric_ROI9_photos=("real_numeric_ROI9","sum"))
    x=q.loc[q.real_numeric_ROI9].groupby("target_id",observed=True).wcvp_accepted_plant_name_id.nunique()
    g["real_WCVP_accepted_taxa_with_actual_ROI9"]=x.reindex(g.index,fill_value=0).astype(int)
    return {
      "original_source_photo_cell_years":int(len(g)),
      "source_cell_years_with_5plus_attempts_3plus_real_WCVP_taxa":int(
          ((g.original_photo_observations>=5)&
           (q.groupby("target_id").wcvp_accepted_plant_name_id.nunique().reindex(g.index)>=3)).sum()),
      "actual_strict_5_attempt_3_success_3_accepted_photo_cell_years":int(
          ((g.original_photo_observations>=5)&(g.actual_numeric_ROI9_photos>=3)&
           (g.real_WCVP_accepted_taxa_with_actual_ROI9>=3)).sum()),
      "with_at_least_one_actual_source_ROI9_photo":int(g.actual_numeric_ROI9_photos.ge(1).sum()),
      "independent_source_islands_represented":int(q.land_id.nunique()),
    }

def run(schedule,receipt,first,batches,out):
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite sealed ALAN real photo source aggregate")
    src=source_schedule(schedule,receipt)
    results=[];receipts=[]
    require(len(batches)==4 and set(batches)=={2,3,4,5},
            "Source original pixel batch groups must include exactly 2-5")
    for k,root in [(1,first),*sorted(batches.items())]:
        d,r=read_authenticated_batch(root,k)
        results.append(d);receipts.append(r)
    combined=pd.concat(results,ignore_index=True)
    require(len(combined)==1280 and combined.photo_id.is_unique and
            combined.observation_id.is_unique and
            combined.batch_number.value_counts().to_dict()=={1:256,2:256,3:256,4:256,5:256},
            "Original source ROI9 photo UUIDs duplicated/lost across batches")
    expected=src.loc[src.source_batch.astype(int).le(5)]
    require(len(expected)==1280 and
            set(combined.photo_id)==set(expected.photo_id) and
            set(combined.observation_id)==set(expected.observation_id),
            "ROI9 measured different photo/observation IDs than frozen original source")
    m=combined.merge(expected[["photo_id","observation_id","land_id","target_id",
             "year","wcvp_accepted_plant_name_id","source_batch"]].rename(
             columns={"land_id":"expected_land_id",
                      "target_id":"expected_target_id","year":"expected_year"}),
             on=["photo_id","observation_id"],how="inner",validate="one_to_one")
    require(len(m)==1280 and
            m.land_id.astype(str).eq(m.expected_land_id).all() and
            m.target_id.astype(str).eq(m.expected_target_id).all() and
            m.year.astype(str).eq(m.expected_year).all() and
            m.source_batch.astype(str).eq(m.batch_number.astype(str)).all(),
            "Real photo colour attached to a different island/year original cell")
    support=actual_colour_site_support(m)
    require(support["source_cell_years_with_5plus_attempts_3plus_real_WCVP_taxa"]==3,
            "Original source-only first5 photo opportunity threshold drifted")
    per_batch=m.groupby("batch_number").agg(
       genuine_original_photos=("photo_id","nunique"),
       genuine_new_observation_UUIDs=("observation_id","nunique"),
       actual_numerical_nine_palette_ROIs=("real_numeric_ROI9","sum"),
       source_islands=("land_id","nunique"),
       source_photo_cell_years=("target_id","nunique")).reset_index()
    summary={
      "schema":SCHEMA,
      "source_original_first_five_photo_batches":list(range(1,6)),
      "new_distinct_photo_attempts":len(m),
      "independent_original_observation_UUIDs":int(m.observation_id.nunique()),
      "distinct_source_GSHHG_islands":int(m.land_id.nunique()),
      "genuinely_measured_nine_palette_source_photos":int(m.real_numeric_ROI9.sum()),
      "original_photo_technical_unknowns":int((~m.real_numeric_ROI9).sum()),
      "strict_photo_cell_support_not_added_to_historical_main":True,
      "original_numerical_ROI9_5photo_3success_3accepted_site_support":support,
      "source_photo_mode":"original FCP YOLO11/EfficientSAM ROI-v4 unchanged",
      "original_FCP_method_git_blob":FCP_MODEL_BLOB,
      "per_batch_actual_original_ROI9_receipts":receipts,
      "original_source_full_8942_sha256":SCHEDULE_SHA,
      "no_original_photo_image_bytes_saved":True,
      "technical_unknown_not_imputed_as_white_or_coloured":True,
      "source_island_first_not_random_global_native_flora_sample":True,
    }
    out.mkdir(parents=True)
    file=out/"actual_ROI9_original_first5_batch_aggregate.csv"
    per_batch.to_csv(file,index=False)
    summary["batch_aggregate_SHA256"]=sha(file)
    (out/"ALAN_1280_actual_original_photo_ROI9_source_receipt.json").write_text(
         json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("ALAN_NEW1280_GENUINE_ORIGINAL_ISLAND_ROI9_MERGE",
          json.dumps({key:summary[key] for key in (
              "new_distinct_photo_attempts","independent_original_observation_UUIDs",
              "distinct_source_GSHHG_islands","genuinely_measured_nine_palette_source_photos",
              "original_photo_technical_unknowns",
              "original_numerical_ROI9_5photo_3success_3accepted_site_support")
          },sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schedule",type=Path,required=True)
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--first",type=Path,required=True)
    for batch in (2,3,4,5):
        p.add_argument(f"--batch{batch}",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.schedule,a.source_receipt,a.first,{b:getattr(a,f"batch{b}") for b in (2,3,4,5)},a.out)
