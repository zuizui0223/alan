#!/usr/bin/env python3
"""Collect 24,576 actual original RGB flower ROI9, preserve three fixed waves.

Only original worker CSV plus signed immutable shard manifests are admissible.
Never substitute technical unknown as a flower colour. Recheck all 192
photo-original shard SHA256 receipts and exactly one original photo per ID.
Outputs are three real-photo wave CSVs, no original licensed RGB pixels.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.acquisition import mois_full_index_real_original_photo_measure as worker
from scripts.acquisition.alan_measure_new256_original_photo_roi9_20261108 import (
    COMPONENTS, FRACTIONS, valid_nine,
)
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

SCHEMA="ALAN_real_original_tier08_additional24576_RGB_ROI9_3waves_v1"
SOURCE="alan_tier08_additional24576_balanced_island_area_light"
PHOTO_N=24576
SHARDS=192
WAVE_N=8192
MODEL_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_TIER08_24576_COLLECT: "+msg)
def sha(path):
    return panel.sha(path)

def collect(source,measured,out):
    source,measured,out=map(Path,(source,measured,out))
    require(not out.exists(),"Refuse overwrite of source SHA-authenticated original RGB result")
    m=json.loads((source/"manifest.json").read_text())
    require(m.get("schema")==worker.SCHEMA and
            m.get("source_group")==SOURCE and
            m.get("colour_measured") is False and
            m.get("new_photo_attempts_planned")==PHOTO_N and
            m.get("work_units")==SHARDS and
            m.get("worker_ids")==list(range(SHARDS)) and
            m.get("waves")==[8192,16384,24576] and
            m.get("original_photos_not_colours_selected_before_pixels") is True,
            "Signed large 24,576 original photo precolour source is not genuine")
    batches=[];receipts=[]
    for i in range(SHARDS):
        plan,_=worker.input_rows(source,i)
        names=list(measured.rglob(f"shard_{i:03d}/receipt.json"))
        require(len(names)==1,f"Missing or repeated actual original RGB worker shard {i}")
        folder=names[0].parent
        rec=json.loads(names[0].read_text())
        file=folder/"photos.csv"
        f=pd.read_csv(file,dtype=str,keep_default_na=False)
        require(
            rec.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1" and
            rec.get("shard")==i and
            rec.get("source_sha256")==m["photo_shard_sha256"][str(i)] and
            rec.get("photos_sha256")==sha(file) and
            rec.get("attempted")==len(plan)==len(f)==128 and
            f.photo_id.is_unique and f.observation_id.is_unique and
            set(f.photo_id)==set(plan.photo_id) and
            set(f.observation_id)==set(plan.observation_id) and
            rec.get("raw_licensed_image_bytes_persisted") is False and
            rec["method"].get("measurement_git_blob")==MODEL_BLOB and
            rec["method"].get("location_blind") is True and
            rec["method"].get("persisted_image_pixels") is False,
            f"Source or actual original FCP image receipt / ID SHA mismatch shard {i}")
        require(set(FRACTIONS)<=set(f.columns),
                f"Missing measured flower-region ROI9 fractions {i}")
        batches.append(f)
        receipts.append(rec)
    photo=pd.concat(batches,ignore_index=True)
    require(len(photo)==PHOTO_N and photo.photo_id.is_unique and
            photo.observation_id.is_unique and
            len({json.dumps(r["method"],sort_keys=True) for r in receipts})==1,
            "All 192 true original image photo IDs, observation UUIDs, or method drift")
    good,numeric=valid_nine(photo)
    photo["strict9"]=good
    for i,col in enumerate(COMPONENTS):
        photo[col]=numeric[:,i]
    photo["strict9_white"]=photo["strict9_white"].where(photo.strict9)
    require(photo.loc[photo.strict9,"strict9_white"].between(0,1).all() and
            photo.loc[~photo.strict9,"strict9_white"].isna().all() and
            photo.loc[~photo.strict9,list(COMPONENTS)].isna().all().all() and
            np.allclose(photo.loc[photo.strict9,list(COMPONENTS)].sum(axis=1),
                        1,atol=1e-6,rtol=0),
            "A technical unknown image acquired artificial numeric flower colour")
    require(photo.wave.astype(int).value_counts().to_dict()==
            {1:WAVE_N,2:WAVE_N,3:WAVE_N},
            "Original precolour wave assignments were changed by RGB processing")
    cols=["photo_id","observation_id","inat_taxon_id",
          "wcvp_accepted_plant_name_id","land_id","target_id","year",
          "observed_on","source_arm","acquisition_arm","area_third",
          "area_km2","viirs_median","log_light","wave",
          "morph","measurement_status","image_sha256",
          "strict9",*COMPONENTS]
    require(set(cols)<=set(photo.columns),"Original photo output dropped WCVP/geography/photo source")
    out.mkdir(parents=True)
    wave_summary=[]
    digests={}
    for wave in (1,2,3):
        df=photo.loc[photo.wave.astype(int).eq(wave)].copy()
        fp=out/f"tier08_additional_8192_real_original_photo_ROI9_wave{wave}.csv.gz"
        df[cols].to_csv(fp,index=False,
            compression={"method":"gzip","mtime":0},lineterminator="\n")
        valid=df.loc[df.strict9]
        site=df.groupby("target_id",observed=True).agg(
            original_photo_attempts=("observation_id","nunique"),
            original_accepted_taxa=("wcvp_accepted_plant_name_id","nunique"),
            land=("land_id","first"),year=("year","first"))
        ok=valid.groupby("target_id",observed=True).agg(
            numeric_ROI9=("photo_id","size"),
            measured_accepted_taxa=("wcvp_accepted_plant_name_id","nunique"))
        site=site.join(ok,how="left").fillna({"numeric_ROI9":0,"measured_accepted_taxa":0})
        strict=(site.original_photo_attempts.ge(5)&site.numeric_ROI9.ge(3)&
                site.measured_accepted_taxa.ge(3))
        row={
           "wave":wave,
           "real_original_photo_attempts":len(df),
           "actual_numeric_original_ROI9_photo_success":int(df.strict9.sum()),
           "source_technical_unknown":int((~df.strict9).sum()),
           "original_image_acquisition_failures":int(df.measurement_status.eq(
                "image_acquisition_failed").sum()),
           "source_photo_distinct_islands":int(df.land_id.nunique()),
           "source_photo_site_years":int(df.target_id.nunique()),
           "source_alone_strict_site_years":int(strict.sum()),
           "source_alone_strict_islands":int(site.loc[strict,"land"].nunique()),
           "photo_SHA256":sha(fp),
        }
        digests[str(wave)]=sha(fp)
        wave_summary.append(row)
    receipt={
       "schema":SCHEMA,
       "real_original_image_photo_attempts":PHOTO_N,
       "actual_numeric_original_9palette_flowerROI_success":int(good.sum()),
       "real_original_technical_unknown_photos":int((~good).sum()),
       "real_original_photo_image_retrieval_failures":int(
           photo.measurement_status.eq("image_acquisition_failed").sum()),
       "real_original_measurement_shards":SHARDS,
       "complete_genuine_source_photo_IDs":int(photo.photo_id.nunique()),
       "complete_genuine_source_observation_UUIDs":int(photo.observation_id.nunique()),
       "original_source_no_copyrighted_RGB_pixels_exported":True,
       "original_fixed_FCP_image_measurement_blob":MODEL_BLOB,
       "signed_source_manifest_sha256":sha(source/"manifest.json"),
       "waves":wave_summary,
       "wave_output_SHA256":digests,
       "source_photo_outcomes_unknown_not_imputed":True,
       "neither_causal_ALAN_nor_native_botanical_species_frequencies_verified":True,
       "historical_numeric_observation_ID_vs_new_UUID_crosswalk_unverified":True,
    }
    require(receipt["actual_numeric_original_9palette_flowerROI_success"]+
            receipt["real_original_technical_unknown_photos"]==PHOTO_N,
            "Actual original RGB image successes/unknown do not conserve photo source count")
    (out/"ALAN_tier08_additional24576_actual_original_RGB_ROI9_signed_receipt.json"
    ).write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("ALAN_TIER08_LARGE_24576_TRUE_ORIGINAL_RGB_ROI9_MEASURED",
          json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source","measured","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    collect(a.source,a.measured,a.out)
if __name__=="__main__":main()
