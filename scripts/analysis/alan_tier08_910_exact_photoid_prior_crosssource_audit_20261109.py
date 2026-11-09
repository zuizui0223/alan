#!/usr/bin/env python3
"""Exact original photo-ID cross-source audit of finished tier08 910 real ROI9.

Read only signed, already-measured source CSVs; never redownload RGB, invent
flower colours, or equate legacy numeric observation IDs with original-source
UUIDs. tier08 can only increase unique PHOTO ID denominator after deleting
identical original photo IDs present in historical canonical, tier06, tier07.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

SCHEMA="ALAN_tier08_910_original_photoID_crosssource_reuse_audit_v1"
SHA={
  "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
  "tier06":"e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd",
  "tier07":"7198455819f1ac260a6cbeacd9f4536c30d766595daa64dc5205c75e9d5a1c62",
  "tier08":"37a329e54b0bbced86e40ed16d45ca095e530aa4faa64b3aa674c29023f8c1f9"
}
N={"canonical":1032882,"tier06":8942,"tier07":27507,"tier08":910}
def require(ok,message):
    if not bool(ok):raise ValueError("ALAN_TIER08_ORIGINAL_PHOTOID_REUSE: "+message)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for p in iter(lambda:f.read(2<<20),b""):h.update(p)
    return h.hexdigest()

def restore_four(canonical,tier06,tier07,tier08):
    paths={"canonical":canonical,"tier06":tier06,"tier07":tier07,"tier08":tier08}
    d={}
    for arm,p in paths.items():
        require(sha(p)==SHA[arm],arm+" signed photo ID bytes drift")
        cols=(["photo_id","observation_id","source_subarm","observed_on"]
              if arm=="canonical" else
              ["photo_id","observation_id","target_id",
               "wcvp_accepted_plant_name_id","year","strict9","acquisition_arm"]
              if arm=="tier08" else
              ["photo_id","observation_id"])
        q=pd.read_csv(p,compression="gzip",dtype=str,keep_default_na=False,
                      usecols=cols)
        require(len(q)==N[arm] and q.photo_id.is_unique and
                q.observation_id.is_unique and
                q.photo_id.str.fullmatch(r"\d+").all(),
                arm+" authentic source photo/observation IDs incomplete")
        d[arm]=q
    require(not set(d["tier06"].photo_id)&set(d["tier07"].photo_id),
            "Old source tier06 and tier07 photographed IDs overlap")
    return d

def strict_sites(frame):
    if frame.empty:
        return {"strict_cell_years":0,
                "strict_pair_island_years":0}
    g=frame.groupby("target_id",observed=True).agg(
        attempts=("observation_id","nunique"),
        area_arm=("acquisition_arm","first"))
    numeric=frame.loc[frame.strict9.astype(str).str.lower().eq("true")]
    success=numeric.groupby("target_id",observed=True).agg(
        roi9=("observation_id","size"),
        measured_accepted_species=("wcvp_accepted_plant_name_id","nunique"))
    g=g.join(success,how="left").fillna(0)
    qualified=(g.attempts.ge(5)&g.roi9.ge(3)&
               g.measured_accepted_species.ge(3))
    # Match the actual source collector's *cell-level* arm definition:
    # first original photo fixes the site acquisition arm. Reference and
    # targeted frames may have different UUIDs in the SAME source site×year,
    # but a site is never counted in both arm-specific strict outcomes.
    valid_targets=set(g.index[qualified & g.area_arm.eq(
        "area_stratified_sameyear_two_site_light_pair")])
    good_pair=frame.loc[frame.target_id.isin(valid_targets)].copy()
    qualifying_pair_cell_years=int(len(valid_targets))
    # Exact island/year key: location landmass from original same site×year.
    good_pair["land_id"]=good_pair.target_id.str.rsplit("|",n=1).str[-1]
    good_pair["year"]=good_pair.target_id.str.split("|",regex=False).str[0]
    n_sites=good_pair.groupby(["land_id","year"]).target_id.nunique()
    return {
       "strict_cell_years":int(qualified.sum()),
       "strict_pair_arm_site_years":qualifying_pair_cell_years,
       "strict_pair_island_years":int(n_sites.ge(2).sum())
    }

def audit(d):
    hist=d["canonical"]
    source06,source07,new=d["tier06"],d["tier07"],d["tier08"]
    hids=set(hist.photo_id)
    hmain=set(hist.loc[hist.source_subarm.eq("historical_main"),"photo_id"])
    prior06=set(source06.photo_id)
    prior07=set(source07.photo_id)
    require(not prior06&prior07,"Source tiers 06 and 07 must not duplicate photos")
    current=set(new.photo_id)
    overlap={
       "exact_photo_in_all_canonical":len(current&hids),
       "exact_photo_in_historical_main":len(current&hmain),
       "exact_photo_in_tier06":len(current&prior06),
       "exact_photo_in_tier07":len(current&prior07),
       "exact_photo_in_at_least_one_previous_arm":len(current&(hids|prior06|prior07)),
    }
    reused=current&(hids|prior06|prior07)
    cleaned=new.loc[~new.photo_id.isin(reused)].copy()
    require(cleaned.photo_id.is_unique and
            not set(cleaned.photo_id)&(hids|prior06|prior07) and
            len(cleaned)+len(reused)==910,
            "Conservative exact original-photo ID union must be disjoint")
    before= strict_sites(new)
    after= strict_sites(cleaned)
    labels={
        "new_original_tier08_real_photo_attempts_before_exact_photo_dedup":len(new),
        "new_original_tier08_real_numeric_roi9_before_dedup":int(new.strict9.str.lower().eq("true").sum()),
        "new_original_tier08_photo_attempts_after_exact_photo_dedup":len(cleaned),
        "new_original_tier08_real_numeric_roi9_after_dedup":int(
             cleaned.strict9.str.lower().eq("true").sum()),
        "cross_source_photo_ID_overlap":overlap,
        "strict_site_support_before":before,
        "strict_site_support_after":after,
        "historical_canonical_observation_id_text_matches_tier08":int(
            len(set(hist.observation_id)&set(new.observation_id))),
        "old_observation_ID_format_to_original_UUID_crosswalk_not_verified":True,
        "remaining_tier08_different_photo_ids_not_proven_independent_biological_observations":True,
        "original_RGB_pixels_downloaded_for_this_audit":0,
        "new_image_ROI9_inferences_for_this_audit":0,
        "original_numerical_photo_unknowns_not_cast_to_white":True,
        "all_photo_IDs_deduplicated_without_colour_favourable_replacement":True,
        "full_sameyear_island_environment_not_automatically_refitted":True,
    }
    require(labels["new_original_tier08_real_numeric_roi9_before_dedup"]==330 and
            before["strict_cell_years"]==43 and
            before["strict_pair_island_years"]==15,
            "Previously authentic tier08 910 original-photo ROI9 outcome drift")
    return labels

def run(canonical,tier06,tier07,tier08,out):
    out=Path(out)
    require(not out.exists(),"Cannot overwrite SHA-authenticated source-only crossarm audit")
    d=restore_four(canonical,tier06,tier07,tier08)
    a=audit(d)
    out.mkdir(parents=True)
    r={
       "schema":SCHEMA,
       "source_photo_counts":N,
       "source_sha256":SHA,
       "audit":a,
       "source_selection_limited_to_2019_2021_2023_2024_2025":True,
       "tiers_06_07_08_separately_colour_blind_preselected":True,
       "not_a_new_positive_ALAN_ecology_effect":True,
    }
    (out/"ALAN_tier08_910_exact_photoid_prior_source_audit.json").write_text(
        json.dumps(r,sort_keys=True,indent=2)+"\n")
    print("ALAN_TIER08_910_AFTER_ALL_HISTORICAL_TIER06_TIER07_ORIGINAL_PHOTOID_DEDUP",
          json.dumps(a,sort_keys=True),flush=True)
    return r

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("canonical","tier06","tier07","tier08","out"):
        p.add_argument("--"+name,type=Path,required=True)
    x=p.parse_args()
    run(x.canonical,x.tier06,x.tier07,x.tier08,x.out)
if __name__=="__main__":
    main()
