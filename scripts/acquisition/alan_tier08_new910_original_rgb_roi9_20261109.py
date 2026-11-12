#!/usr/bin/env python3
"""Actually measure the 910 SHA-frozen novel tier08 original RGB photographs.

The 256 reference and 654 area/light-paired photos have source-verified UNIQUE
observation UUIDs and photo IDs; each tier08 observation is independent of
historical-main, tier06, tier07. Original FCP YOLO11/EfficientSAM ROI-v4 uses
UNCHANGED Git blob b1d7091... ; unknown results are NEVER filled with white.

'prepare' validates and freezes eight bounded <=128-photo shards.
'collect' verifies all real image worker receipts, numerical ROI9, and
site-year/taxon support, retaining actual image failures as technical unknowns.
No copyrighted original RGB image bytes are distributed in resulting artifacts.
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
    PALETTE, COMPONENTS, FRACTIONS, valid_nine,
)

QUEUE_SHA = "9d2e60b9c9bdb2fb56c1c17dd69fca5b6c606d1e64d707020d070dd82e49a552"
MASTER_SHA = "b34a8a564dc9ba3f02934e1e04077b4655cf8a9333ecefcb6bc93c5d9dc757fd"
SOURCE_SCHEMA = "ALAN_tier08_precolour_independent_uuid_balanced_light_pair_source_v1"
WORKER_SCHEMA = "mois_full_source_observation_sampling_cohort_v1"
MODEL_BLOB = "b1d7091db44682559bbe4b48a9f3b8db648501a6"
SOURCE_ARM = "alan_93m_tier08_new_independent_island_original_photo"
REF_ARM = "hash_random_within_envmatched_tier08_reference"
PAIR_ARM = "area_stratified_sameyear_two_site_light_pair"
COUNT = 910
REF_COUNT = 256
PAIR_COUNT = 654
PHOTO_UNIT = 128
UNITS = 8

def require(ok, message):
    if not bool(ok):
        raise ValueError("ALAN_TIER08_REAL_ORIGINAL_RGB_ROI9: " + message)

def sha(file):
    h = hashlib.sha256()
    with Path(file).open("rb") as f:
        for part in iter(lambda: f.read(4 << 20), b""):
            h.update(part)
    return h.hexdigest()

def load_source(queue, receipt):
    require(sha(queue) == QUEUE_SHA,
            "Frozen original tier08 precolour reference+pair photo queue SHA drift")
    r = json.loads(Path(receipt).read_text())
    cnt = r.get("observation_counts", {})
    require(r.get("schema") == SOURCE_SCHEMA and
            r.get("queue_sha256") == QUEUE_SHA and
            r.get("master_sha256") == MASTER_SHA and
            cnt.get("total_precolour_source_photo_queue") == COUNT and
            cnt.get("reference_actual_original_photos_unmeasured") == REF_COUNT and
            cnt.get("source_paired_actual_original_photos_unmeasured") == PAIR_COUNT and
            cnt.get("new_flower_ROI9_measurements") == 0 and
            r.get("no_unmeasured_colour_inferred") is True,
            "tier08 original source manifest not genuine precolour data")
    q = pd.read_csv(queue, compression="gzip", dtype=str, keep_default_na=False)
    require(len(q) == COUNT and q.photo_id.is_unique and
            q.observation_id.is_unique and
            q.photo_id.str.fullmatch(r"\d+").all() and
            q.source_arm.eq(SOURCE_ARM).all() and
            q.selection_probability.eq("").all() and
            q.photo_url_large.str.fullmatch(
                r"https://inaturalist-open-data\.s3\.amazonaws\.com/photos/\d+/large\.(jpg|jpeg|png)"
            ).all() and
            q.wcvp_accepted_plant_name_id.ne("").all() and
            q.year.isin(("2019", "2021", "2023", "2024", "2025")).all(),
            "Signed photo source uniqueness, taxonomy, years or original URL changed")
    split = q.acquisition_arm.value_counts().to_dict()
    require(split == {REF_ARM: REF_COUNT, PAIR_ARM: PAIR_COUNT},
            "Reference and same-year paired acquisition arms changed")
    pair = q.loc[q.acquisition_arm.eq(PAIR_ARM)].copy()
    require(pair.land_id.nunique() == 28 and pair.target_id.nunique() == 56,
            "Prior selection of 28 independent paired islands / 56 site-years changed")
    require(pair.groupby("land_id").target_id.nunique().eq(2).all() and
            pair.groupby("land_id").year.nunique().eq(1).all(),
            "Paired island source must have exactly two genuine sites IN THE SAME YEAR")
    return q, r

def assign_shards(data):
    """One single immutable 8-shard source roster, last unit only 14 photos."""
    d = data.copy().reset_index(drop=True)
    require(len(d) == COUNT, "Frozen tier08 photo count changed")
    d["measurement_id"] = d.photo_id.map(lambda x: "alan_tier08_910_" + x)
    d["image_filename"] = d.measurement_id + ".jpg"
    d["source_snapshot_date"] = "2026-09-27"
    d["shard"] = (np.arange(len(d)) // PHOTO_UNIT).astype(int)
    sizes = d.groupby("shard").size().to_dict()
    require(sizes == {**{i:128 for i in range(7)},7:14},
            "Expected seven complete original-photo units and last 14-photo unit")
    return d

def prepare(queue, receipt, out):
    out = Path(out)
    require(not out.exists(), "Never overwrite frozen tier08 measurement source")
    q, source_receipt = load_source(queue, receipt)
    q = assign_shards(q)
    (out / "shards").mkdir(parents=True)
    hashes={}
    for i,group in q.groupby("shard",sort=True):
        p = out / "shards" / f"shard-{i:03d}.csv"
        group.to_csv(p,index=False,lineterminator="\n")
        hashes[str(i)] = sha(p)
    m = {
        "schema": WORKER_SCHEMA,
        "colour_measured": False,
        "source_group": "alan_tier08_910_source_reference_256_area_light_paired_654",
        "source_arm": SOURCE_ARM,
        "new_photo_attempts_planned": COUNT,
        "work_units": UNITS,
        "worker_ids": list(range(UNITS)),
        "photo_shard_sha256": hashes,
        "reference_original_photo_UUIDs": REF_COUNT,
        "area_light_pair_original_photo_UUIDs": PAIR_COUNT,
        "source_precolour_queue_SHA256": QUEUE_SHA,
        "source_full_photo_master_SHA256": MASTER_SHA,
        "tier08_source_registry_SHA256":
            source_receipt["source_registry_tier08_audit"]["full_registry_sha256"],
        "original_photos_not_colours_selected_before_pixel_download": True,
        "photo_native_wild_flowering_not_independently_verified": True,
        "sampling_not_global_unbiased_native_flora": True,
    }
    (out / "manifest.json").write_text(
        json.dumps(m,sort_keys=True,indent=2)+"\n")
    print("ALAN_TIER08_NEW910_PRECOLOUR_8_ORIGINAL_PHOTO_WORK_SHARDS_FROZEN",
          json.dumps({"original_photo_UUIDs":COUNT,
            "reference":REF_COUNT,"sameyear_light_pairs":PAIR_COUNT,
            "source_shard_hashes":hashes},sort_keys=True),flush=True)
    return m

def strict_island_sites(measured):
    valid = measured.loc[measured.strict9].copy()
    total = measured.groupby("target_id",observed=True).agg(
        attempts=("observation_id","nunique"),
        area_arm=("acquisition_arm","first"),
        original_taxa=("wcvp_accepted_plant_name_id","nunique"),
        land_id=("land_id","first"),
        year=("year","first"),
    )
    outcomes=valid.groupby("target_id",observed=True).agg(
        true_roi9_success=("observation_id","size"),
        true_measured_taxa=("wcvp_accepted_plant_name_id","nunique"))
    combined=total.join(outcomes,how="left")
    for col in ("true_roi9_success","true_measured_taxa"):
        combined[col]=combined[col].fillna(0).astype(int)
    combined["strict_real"]=combined.attempts.ge(5)&(
        combined.true_roi9_success.ge(3))&combined.true_measured_taxa.ge(3)
    p=combined.loc[combined.area_arm.eq(PAIR_ARM)].reset_index()
    strict=p.loc[p.strict_real]
    both=strict.groupby(["land_id","year"],observed=True).target_id.nunique()
    return {
        "original_source_target_site_years":int(total.shape[0]),
        "real_strict_photo_cell_years_all_arms":int(combined.strict_real.sum()),
        "real_strict_photo_cell_years_targeted_pairs":int(len(strict)),
        "targeted_island_years_with_two_strict_spatial_sites":int(both.ge(2).sum()),
        "reference_actual_strict_photo_cell_years":int(combined.loc[
            combined.area_arm.eq(REF_ARM),"strict_real"].sum()),
    }

def collect(source, measured, out):
    source, measured, out=map(Path,(source, measured, out))
    require(not out.exists(), "Never overwrite genuinely measured tier08 ROI9")
    m=json.loads((source/"manifest.json").read_text())
    require(m.get("schema")==WORKER_SCHEMA and
            m.get("source_precolour_queue_SHA256")==QUEUE_SHA and
            m.get("worker_ids")==list(range(UNITS)),
            "Signed tier08 original-image work manifest drift")
    frames=[];receipts=[]
    for i in range(UNITS):
        fixed,_=worker.input_rows(source,i)
        folder=measured/f"shard_{i:03d}"
        r=json.loads((folder/"receipt.json").read_text())
        csvfile=folder/"photos.csv"
        df=pd.read_csv(csvfile,dtype=str,keep_default_na=False)
        require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1" and
                r.get("shard")==i and r.get("attempted")==len(fixed) and
                r.get("source_sha256")==m["photo_shard_sha256"][str(i)] and
                r.get("photos_sha256")==sha(csvfile) and
                len(df)==len(fixed) and df.photo_id.is_unique and
                df.observation_id.is_unique and
                set(df.photo_id)==set(fixed.photo_id) and
                set(df.observation_id)==set(fixed.observation_id),
                "Measured original RGB worker photo ID/observation/shard SHA mismatch")
        require(set(FRACTIONS)<=set(df.columns),
                "Original FCP did not produce genuine nine component flower ROI data")
        frames.append(df);receipts.append(r)
    raw=pd.concat(frames,ignore_index=True)
    require(len(raw)==COUNT and raw.photo_id.is_unique and
            raw.observation_id.is_unique and
            len({json.dumps(r["method"],sort_keys=True) for r in receipts})==1,
            "Genuine source original photo UUIDs or measurement methods mixed")
    method=receipts[0]["method"]
    require(method.get("measurement_git_blob")==MODEL_BLOB and
            method.get("location_blind") is True and
            method.get("persisted_image_pixels") is False,
            "Changed original YOLO11 EfficientSAM evaluator or persisted licensed RGB")
    good, numeric=valid_nine(raw)
    raw["strict9"]=good
    for i,c in enumerate(COMPONENTS):
        raw[c]=numeric[:,i]
    require(raw.loc[raw.strict9,list(COMPONENTS)].notna().all().all() and
            raw.loc[~raw.strict9,list(COMPONENTS)].isna().all().all() and
            np.allclose(raw.loc[raw.strict9,list(COMPONENTS)].sum(axis=1),1,
                        atol=1e-6,rtol=0),
            "Technical unknown was filled or true nine colour masses not compositional")
    raw["strict9_white"]=raw["strict9_white"].where(raw.strict9)
    require(raw.loc[raw.strict9,"strict9_white"].between(0,1).all() and
            raw.loc[~raw.strict9,"strict9_white"].isna().all(),
            "Original photographic whiteness should exist only for true ROI9 successes")
    status={str(k):int(v) for k,v in raw.measurement_status.value_counts().items()}
    site_support=strict_island_sites(raw)
    out.mkdir(parents=True)
    cols=["photo_id","observation_id","inat_taxon_id",
          "wcvp_accepted_plant_name_id","land_id","target_id","year",
          "observed_on","acquisition_arm","area_third","source_arm",
          "source_light_span","area_km2","viirs_median","log_light",
          "morph","measurement_status","image_sha256",
          "strict9","strict9_white",*COMPONENTS]
    require(set(cols)<=set(raw.columns), "Genuine original worker lost source identifiers")
    outpath=out/"alan_tier08_new910_actual_original_photo_ROI9.csv.gz"
    raw[cols].to_csv(outpath,index=False,
         compression={"method":"gzip","mtime":0},lineterminator="\n")
    report={
        "schema":"ALAN_tier08_910_real_original_photo_FCP_ROI9_v1",
        "original_source_precolour_queue_SHA256":QUEUE_SHA,
        "original_photo_attempts":COUNT,
        "independent_new_original_observation_UUIDs":int(raw.observation_id.nunique()),
        "real_numeric_original_flower_ROI9_success":int(good.sum()),
        "real_original_technical_unknown_photos":int((~good).sum()),
        "real_original_image_retrieval_failures":int(raw.measurement_status.eq(
            "image_acquisition_failed").sum()),
        "source_reference_photo_attempts":REF_COUNT,
        "source_sameyear_pair_photo_attempts":PAIR_COUNT,
        "source_distinct_islands":int(raw.land_id.nunique()),
        "source_photo_measurement_statuses":status,
        **site_support,
        "genuine_source_original_photo_CSV_sha256":sha(outpath),
        "original_source_manifest_sha256":sha(source/"manifest.json"),
        "model_source_blob":MODEL_BLOB,
        "frozen_image_model_method":method,
        "signed_original_128photo_worker_receipts":UNITS,
        "raw_original_rgb_image_bytes_persisted":False,
        "source_colour_unknown_never_imputed_as_white_or_nonwhite":True,
        "source_island_targeted_not_global_native_flora_prevalence":True,
        "no_assumed_night_pollination_or_causal_flower_evolution":True,
    }
    require(report["real_numeric_original_flower_ROI9_success"]+
            report["real_original_technical_unknown_photos"]==COUNT,
            "Source photo RGB outcomes do not add to actual source attempt count")
    (out/"ALAN_tier08_new910_real_original_ROI9_receipt.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_TIER08_NEW910_REAL_ORIGINAL_RGB_FCP_ROI9_FULL_SOURCE_VERIFIED",
          json.dumps(report,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    cmd=p.add_subparsers(dest="command",required=True)
    q=cmd.add_parser("prepare")
    for x in ("queue","receipt","out"):
        q.add_argument("--"+x,type=Path,required=True)
    q=cmd.add_parser("collect")
    for x in ("source","measured","out"):
        q.add_argument("--"+x,type=Path,required=True)
    a=p.parse_args()
    if a.command=="prepare":
        prepare(a.queue,a.receipt,a.out)
    else:
        collect(a.source,a.measured,a.out)

if __name__=="__main__":
    main()
