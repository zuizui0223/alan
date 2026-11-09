#!/usr/bin/env python3
"""ALAN: genuinely ADD tier08 910 measured original RGB ROI9 to old+new pooled source.

Import and preserve the 2026-11-09 rigorously reconstructed historical-main
and tier06+tier07 true original photo ROI9 result. Use the identical original
GSHHG exact photo cell/year, WCVP accepted species, original same-year VIIRS/
WorldPop, climate, and original GSHHG mainland coast vertex distances. Nothing
is inferred from source-only photo opportunities or unclassified images.

The 910 actually measured tier08 original photos form a source with a different
acquisition design (256 reference, 654 area-stratified same-year light pairs).
Deduplicate by exact PHOTO IDs against all original photos; legacy observation
IDs and new UUIDs are not proven joinable without a signed crosswalk. Never
claim 910 globally independent biological observations without that crosswalk.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
import numpy as np
from scripts.analysis import alan_pooled_historical_new36449_photoROI9_20261109 as base
from scripts.analysis import alan_new36449_photoROI9_island_light_replication_20261109 as fresh
from scripts.analysis import alan_isolation_area_light_20261108 as isolation
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

SCHEMA="ALAN_historical_main_tier06_tier07_tier08_910_all_original_photo_pooled_v1"
TIER08_SOURCE_SHA="37a329e54b0bbced86e40ed16d45ca095e530aa4faa64b3aa674c29023f8c1f9"
TIER08_SOURCE_QUEUE_SHA="9d2e60b9c9bdb2fb56c1c17dd69fca5b6c606d1e64d707020d070dd82e49a552"
TIER08_SCHEMA="ALAN_tier08_910_real_original_photo_FCP_ROI9_v1"
FCP_MODEL_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
TIER08_PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
TIER08_COMPONENTS=tuple("strict9_"+x for x in TIER08_PALETTE)
TIER08_ARM="tier08_910_balanced_reference_and_light_pairs"
AREA_TERM="log_light:log10_area_c"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_TIER08_NEW_PHOTO_POOLED: "+msg)

def load_real_tier08(photo_csv, receipt):
    photo_csv,receipt=Path(photo_csv),Path(receipt)
    require(panel.sha(photo_csv)==TIER08_SOURCE_SHA,
            "Genuine original tier08 910-photo numerical source SHA changed")
    r=json.loads(receipt.read_text())
    require(r.get("schema")==TIER08_SCHEMA and
            r.get("genuine_source_original_photo_CSV_sha256")==TIER08_SOURCE_SHA and
            r.get("original_source_precolour_queue_SHA256")==TIER08_SOURCE_QUEUE_SHA and
            r.get("original_photo_attempts")==910 and
            r.get("independent_new_original_observation_UUIDs")==910 and
            r.get("real_numeric_original_flower_ROI9_success")==330 and
            r.get("real_original_technical_unknown_photos")==580 and
            r.get("signed_original_128photo_worker_receipts")==8 and
            r["frozen_image_model_method"]["measurement_git_blob"]==FCP_MODEL_BLOB and
            r["raw_original_rgb_image_bytes_persisted"] is False,
            "Tier08 original photos not image measured and source signed")
    d=pd.read_csv(photo_csv,compression="gzip",dtype=str,keep_default_na=False)
    needed={"photo_id","observation_id","land_id","target_id","year",
            "wcvp_accepted_plant_name_id","strict9","strict9_white",
            "acquisition_arm",*TIER08_COMPONENTS}
    require(needed<=set(d.columns) and len(d)==910 and d.photo_id.is_unique
            and d.observation_id.is_unique,
            "Tier08 original image photo IDs, accepted taxa, years or measured ROI missing")
    good=d.strict9.str.lower().eq("true")
    white=pd.to_numeric(d.strict9_white,errors="coerce")
    components=d.loc[good,list(TIER08_COMPONENTS)].apply(pd.to_numeric,errors="coerce")
    require(int(good.sum())==330 and
            white.loc[good].between(0,1).all() and
            white.loc[~good].isna().all() and
            components.notna().all().all() and
            np.allclose(components.sum(axis=1),1,atol=1e-6,rtol=0) and
            d.wcvp_accepted_plant_name_id.ne("").all(),
            "Measured flower ROI numeric eligibility or technical missingness changed")
    # Collector exports strict9_white both as named component and original
    # white-fraction field; the two must be identical when present.
    if "strict9_white.1" in d:
        duplicate=pd.to_numeric(d["strict9_white.1"],errors="coerce")
        require(np.allclose(white.loc[good],duplicate.loc[good],
                            atol=1e-10,rtol=0),
                "Two recorded real white ROI columns disagree")
    require(d.target_id.eq(d.year+"|"+
       d.target_id.str.split("|",regex=False).str[1]+"|"+d.land_id).all(),
       "Source tier08 original site/year GSHHG identity changed")
    q=d[["photo_id","observation_id","target_id","land_id","year",
         "wcvp_accepted_plant_name_id"]].copy()
    q["year"]=q.year.astype(int)
    q["real_numeric_ROI9"]=good.astype(bool)
    q["source_photo_white_fraction"]=white
    q["photo_source"]=TIER08_ARM
    return q,{
       "actual_tier08_source_photos":len(q),
       "actual_tier08_numeric_9ROI_success":int(good.sum()),
       "source_reference_photo_attempts":int(
           d.acquisition_arm.eq("hash_random_within_envmatched_tier08_reference").sum()),
       "source_light_pair_photo_attempts":int(
           d.acquisition_arm.eq("area_stratified_sameyear_two_site_light_pair").sum()),
       "photo_source_original_ROI9_SHA256":TIER08_SOURCE_SHA,
       "frozen_fcp_original_evaluator_git_blob":FCP_MODEL_BLOB,
    }

def same_photo_reconciliation(prior,tier08):
    common=set(prior.photo_id)&set(tier08.photo_id)
    collisions=prior.loc[prior.photo_id.isin(common),[
       "photo_id","target_id","land_id","year","wcvp_accepted_plant_name_id"
    ]].merge(
       tier08.loc[tier08.photo_id.isin(common),[
          "photo_id","target_id","land_id","year","wcvp_accepted_plant_name_id"]],
       on="photo_id",how="inner",validate="one_to_one",suffixes=("_prior","_tier08"))
    for key in ("target_id","land_id","year","wcvp_accepted_plant_name_id"):
        require(collisions[key+"_prior"].equals(collisions[key+"_tier08"]),
                "Same original photo ID disagrees on source site, year or WCVP taxon")
    out=tier08.loc[~tier08.photo_id.isin(common)].copy()
    require(len(out)==len(tier08)-len(common) and
            not set(prior.photo_id).intersection(set(out.photo_id)),
            "Tier08 historical/new original photo collisions were counted twice")
    return out,{"tier08_original_photo_ids_already_in_prior":len(common),
                "tier08_original_photo_ids_truly_additional":len(out),
                "tier08_observation_ID_format_independence_from_historical_unproven":True}

def three_arm_equal_source_photo_white(arm1,arm2,arm3,eligible):
    """Every present photo-source arm gets equal weight WITHIN each cell-year."""
    sources=[("old",arm1),("tier06_07",arm2),("tier08",arm3)]
    q=eligible.copy()
    ncols=[]
    for key,d in sources:
        allcells,_=base.photo_cell_stats(d)
        name="source_"+key+"_white"
        q=q.merge(allcells[["target_id","white_equal"]].rename(
            columns={"white_equal":name}),on="target_id",how="left",
            validate="one_to_one")
        ncols.append(name)
    support=q[ncols].notna().sum(axis=1)
    require(support.ge(1).all(),"No authentic measured source for eligible site")
    q["white_pp"]=100*q[ncols].mean(axis=1)
    return q,{"cells_with_3_measured_sources":int(support.eq(3).sum()),
              "cells_with_2_measured_sources":int(support.eq(2).sum()),
              "cells_with_1_measured_source":int(support.eq(1).sum())}

def run(paths,coast,source06,source07,prior06,bulk06,pilot07,bulk07,tier08,tier08_receipt,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwrite of signed real expanded photo pool")
    canonical,wcvp=panel.load_frozen_photos(paths)
    old,audit=base.old_original_photos(canonical,wcvp,paths)
    new,_=base.new_original_photos(source06,source07,prior06,bulk06,pilot07,bulk07)
    require(len(old)==50829 and len(new)==36449,
            "Historical-main or prior tier06+07 genuine photographic source changed")
    shared=set(old.photo_id)&set(new.photo_id)
    require(len(shared)==329,"Old/new source observed photo-ID collision count drift")
    new_clean=new.loc[~new.photo_id.isin(shared)].copy()
    prior=pd.concat([old,new_clean],ignore_index=True)
    require(len(prior)==86949 and prior.photo_id.is_unique,
            "Old+new original photo ID-deduplicated denominator changed")
    tier08all,source_report=load_real_tier08(tier08,tier08_receipt)
    tier08clean,cross_report=same_photo_reconciliation(prior,tier08all)
    require(tier08clean.observation_id.is_unique and
            tier08clean.photo_id.is_unique,
            "New genuine tier08 original photographic identifiers duplicated")
    combined=pd.concat([prior,tier08clean],ignore_index=True)
    require(combined.photo_id.is_unique,"Repeated exact photographic IDs after all arms")
    environment=fresh.original_sameyear_environment(paths)
    _,before_strict=base.photo_cell_stats(prior)
    before, before_audit=base.environment_qualified(before_strict,environment)
    _,after_strict=base.photo_cell_stats(combined)
    after,after_audit=base.environment_qualified(after_strict,environment)
    require(len(before)==2503 and before.land_id.nunique()==193,
            "Prior successful original photo pooled 2,503-site result changed")
    require(len(after)>=len(before) and
            set(before.target_id).issubset(set(after.target_id)),
            "Adding true original photos lost previously strict site-years")
    recent=new_clean if tier08clean.empty else pd.concat([new_clean,tier08clean],
        ignore_index=True)
    _,recentstrict=base.photo_cell_stats(recent)
    recentenv,recentaudit=base.environment_qualified(recentstrict,environment)
    alternate,support=three_arm_equal_source_photo_white(
        old,new_clean,tier08clean,after)
    land=isolation.load_signed_land(coast)
    meta=land[["land_id","land_class","area_km2"]]
    distances=isolation.isolated_mainland_distance_km(
        land,set(before.land_id)|set(after.land_id))
    fits={}
    for name,sites in (
       ("previous_true_photo_pool",before),
       ("tier08_incremental_original_photos_in_combined_pool",after),
       ("new_source_tier06_07_plus_tier08_only",recentenv),
       ("pooled_three_source_equal_contribution_per_site",alternate)):
        fits[name]=base.fit_photo_arm(sites,meta,distances,name)
    reference=fits["previous_true_photo_pool"]["models"]["joint_island_equal"][
        "terms"][AREA_TERM]["beta"]
    require(abs(reference-0.5777987743303025)<1e-5,
            "Prior successful 86,949-photo regression not numerically reproduced")
    overlap_counts={
        "old_historical_main_original_photo_attempts":len(old),
        "earlier_new_tier06_07_photo_records":len(new),
        "earlier_old_new_exact_original_photo_ID_collision":len(shared),
        "previous_source_deduplicated_original_photos":len(prior),
        "tier08_real_original_photo_attempts":len(tier08all),
        "tier08_genuine_additional_photo_ids":len(tier08clean),
        "all_source_photoID_unique_original_photos":len(combined),
        "original_numeric_ROI9_in_previous_photo_pool":int(prior.real_numeric_ROI9.sum()),
        "original_numeric_ROI9_new_tier08_photoID_clean":int(tier08clean.real_numeric_ROI9.sum()),
        "total_numeric_original_ROI9_photo_source":int(combined.real_numeric_ROI9.sum()),
        "strict_island_site_years_before_tier08":len(before),
        "strict_island_site_years_after_tier08":len(after),
        "strict_island_site_years_gained_from_tier08":int(len(after)-len(before)),
        "distinct_islands_with_strict_after_tier08":int(after.land_id.nunique()),
        "independent_spatial_islands_in_pooled_tier08_model":int(
            fits["tier08_incremental_original_photos_in_combined_pool"].get(
            "n_spatial_multisite_islands",0)),
        "within_island_spatial_cells_after_tier08":int(
            fits["tier08_incremental_original_photos_in_combined_pool"].get(
            "n_spatial_multisite_cells",0)),
    }
    report={
      "schema":SCHEMA,
      "source_years":[2019,2021,2023,2024,2025],
      "counts":overlap_counts,
      "tier08_910_actual_original_RGB_measured_source_report":source_report,
      "tier08_photo_id_collision_vs_prior_source_report":cross_report,
      "same_cell_source_equal_weight_sensitivity":support,
      "same_year_VIIRS_WorldPop_source_quality_before":before_audit,
      "same_year_VIIRS_WorldPop_source_quality_after":after_audit,
      "same_year_VIIRS_WorldPop_source_quality_new_only":recentaudit,
      "models":fits,
      "baseline_old_new_original_photo_id_pool_numerically_reproduced":True,
      "old_historical_observation_ID_vs_new_UUID_crosswalk_still_unresolved":True,
      "source_photo_not_unbiased_native_flora_or_night_pollination":True,
      "no_colours_imputed_and_no_images_re_downloaded_in_this_join":True,
      "source_output_genuine_original_tier08_photo_sha256":TIER08_SOURCE_SHA
    }
    out.mkdir(parents=True)
    path=out/"ALAN_actual_historical_36449_tier08_910_photoID_unique_pooled.json"
    path.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    short={"counts":overlap_counts,
           "area_interactions":{
             label:({modelname:r["models"][modelname]["terms"][AREA_TERM]
                     for modelname in ("joint_island_equal","joint_cell_equal")}
                    if r["estimable"] else {"estimable":False})
             for label,r in fits.items()}}
    print("ALAN_HISTORICAL_NEW36449_PLUS_TIER08_910_ACTUAL_ORIGINAL_PHOTO_POOL",
          json.dumps(short,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("canonical","wcvp","geo2019","viirs2019","pop2019",
                "geo2021","viirs2021","pop2021","geo2023","viirs2023","pop2023",
                "geo2024","cn2024","uc2024","uc2024receipt",
                "geo2025","geo2025receipt","viirs2025"):
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    for name in ("coast","source06","source07","prior06","bulk06",
                 "pilot07","bulk07","tier08","tier08-receipt","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    source_keys=("canonical","wcvp","geo2019","viirs2019","pop2019",
                 "geo2021","viirs2021","pop2021","geo2023","viirs2023","pop2023",
                 "geo2024","cn2024","uc2024","uc2024receipt",
                 "geo2025","geo2025receipt","viirs2025")
    run({k:getattr(a,k) for k in source_keys},a.coast,a.source06,
        a.source07,a.prior06,a.bulk06,a.pilot07,a.bulk07,
        a.tier08,a.tier08_receipt,a.out)
