#!/usr/bin/env python3
"""Source-backed ALAN intake audit of measured original 1,181,146-photo ROI9 union.

Find genuine EXTRA independent observations, not repeated photo views of one
previous observation, and determine which REAL island photo cell-years have
existing same-year VIIRS + WorldPop coverage versus require new exposure joins.

This is NOT a missing-colour sensitivity model or a new flora prevalence
estimate. Distinguish source design from the frozen historical-main study.
Never publish individual photo/observation IDs or coordinates in results.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import pandas as pd
import numpy as np
import geopandas as gpd
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import prepare_land,assign_points
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import cell_id
from scripts.analysis.mois_global2025_original_source_20261107 import COAST_SHA

RESTORED_SHA="faf6951505f1002a1327f64a2ced13431b562cf13893707125fafff918160f62"
YEARS=(2019,2021,2023,2024,2025)
STATES=("island","mainland")
EXPECTED_ORIGINAL=1181146
EXPECTED_STRICT=451643

def photo_source(archive, canonical):
    panel.require(panel.sha(archive)==RESTORED_SHA,
                  "The original georestored v2 ROI9 archive SHA drifted")
    wanted=["photo_id","observation_id","inat_taxon_id",
            "observed_on","latitude","longitude","source_design",
            "strict9","strict9_white"]
    df=pd.read_parquet(archive,columns=wanted)
    panel.require(len(df)==EXPECTED_ORIGINAL,
                  "Original 1,181,146 archive is incomplete")
    for key in ("photo_id","observation_id","inat_taxon_id"):
        df[key]=df[key].fillna("").astype(str)
    panel.require(df.photo_id.is_unique and int(df.strict9.sum())==EXPECTED_STRICT,
                  "Duplicated original photo IDs or strict9 source count drift")
    photo_ids=set(canonical.photo_id.astype(str))
    observation_ids=set(canonical.observation_id.astype(str))
    source_overlap=int(df.photo_id.isin(photo_ids).sum())
    panel.require(source_overlap>=900000,
                  "Expected original full-source and canonical photo overlap absent")
    extra=df.loc[~df.photo_id.isin(photo_ids)].copy()
    extra["new_independent_observation"]=~extra.observation_id.isin(observation_ids)
    extra["photo_year"]=pd.to_numeric(
        extra.observed_on.fillna("").astype(str).str.slice(0,4),
        errors="coerce")
    extra["latitude"]=pd.to_numeric(extra.latitude,errors="coerce")
    extra["longitude"]=pd.to_numeric(extra.longitude,errors="coerce")
    extra["usable_date"]=extra.photo_year.isin(YEARS)
    extra["usable_coords"]=extra.latitude.between(-90,90) & extra.longitude.between(-180,180)
    extra["real_strict9"]=extra.strict9.astype(bool) & extra.strict9_white.notna()
    extra["strict9_white"]=pd.to_numeric(extra.strict9_white,errors="coerce")
    panel.require(extra.loc[extra.real_strict9,"strict9_white"].between(0,1).all(),
                  "Original numerical strict9 white fraction is invalid")
    source_counts={
        "signed_original_attempted_photos":len(df),
        "signed_original_strict9_photos":int(df.strict9.sum()),
        "already_used_canonical_photo_ids":source_overlap,
        "extra_original_photo_ids_not_in_canonical":len(extra),
        "extra_photo_ids_with_new_observation_uuid":int(extra.new_independent_observation.sum()),
        "extra_true_strict9_photos":int(extra.real_strict9.sum()),
        "extra_strict9_with_new_observation_uuid":int((extra.real_strict9&
            extra.new_independent_observation).sum()),
        "extra_new_observation_photos_with_supported_date_and_coordinates":int((
            extra.new_independent_observation & extra.usable_date & extra.usable_coords).sum()),
        "extra_new_observation_real9_with_supported_date_and_coordinates":int((
            extra.new_independent_observation & extra.usable_date & extra.usable_coords &
            extra.real_strict9).sum()),
        "other_photo_views_of_original_canonical_observations_not_independent":int((
            ~extra.new_independent_observation).sum())
    }
    return extra,source_counts

def source_wcvp_accepted(extra,wcvp):
    ref=wcvp.loc[wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
                 ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    panel.require(ref.inat_taxon_id.is_unique,"Original WCVP accepted source duplicates")
    ref["inat_taxon_id"]=ref.inat_taxon_id.astype(str)
    return extra.merge(ref,on="inat_taxon_id",how="left",validate="many_to_one")

def geography_with_actual_landmasses(extra,coast):
    source=Path(coast)
    receipt=json.loads((source/"source_receipt.json").read_text())
    panel.require(receipt.get("original_archive_sha256")==COAST_SHA,
                  "Original signed GSHHG archive SHA drift")
    for rel,checksum in receipt["extracted_sha256"].items():
        panel.require(panel.sha(source/rel)==checksum,
                      "Original GSHHG source fragment drift")
    land=prepare_land(gpd.read_file(source/"GSHHS_h_L1.shp"))
    lakes=gpd.read_file(source/"GSHHS_h_L2.shp")
    q=extra.loc[extra.new_independent_observation & extra.usable_date &
                extra.usable_coords,["photo_id","observation_id","inat_taxon_id",
                 "latitude","longitude","source_design","photo_year",
                 "real_strict9","strict9_white","wcvp_accepted_plant_name_id"]].copy()
    q=q.rename(columns={"photo_id":"unit_id"})
    panel.require(q.unit_id.is_unique and q.observation_id.is_unique,
                  "Source extras still duplicate observation UUIDs")
    placed=assign_points(q,land,lakes)
    q=q.merge(placed[["unit_id","land_id","land_class"]],on="unit_id",
              validate="one_to_one")
    q["cell_id"]=cell_id(q.longitude,q.latitude)
    good=q.land_class.isin(STATES)&q.land_id.notna()
    q=q.loc[good].copy()
    q["year"]=q.photo_year.astype(int)
    q["land_id"]=q.land_id.astype(str)
    q["target_id"]=(q.year.astype(str)+"|"+q.cell_id+"|"+q.land_id)
    panel.require(q.target_id.str.count("\\|").eq(2).all(),
                  "Original same-year GSHHG photo cell target malformed")
    return q

def summarize(source_photo,existing_year_frames):
    available=[]
    for year, (full,strict) in existing_year_frames.items():
        a=full.copy()
        a["year"]=int(year)
        a["year_and_cell"]=a.target_id.astype(str).str.rsplit("|",n=1).str[0]
        available.append(a)
    orig=pd.concat(available,ignore_index=True)
    panel.require(orig.target_id.is_unique,
                  "Duplicated original same-year VIIRS/WorldPop geocell key")
    statuses=orig.population_status if "population_status" in orig else orig.population_source_status
    density=orig.population_density_per_km2 if "population_density_per_km2" in orig else orig.population_density
    orig["env_ok"]=orig.viirs_valid_months.ge(3)&orig.viirs_median.notna()&statuses.eq("OK")&density.notna()
    orig["climate_ok"]=orig[["bio5","bio12","elevation"]].notna().all(axis=1)
    env=orig[["target_id","year","land_class","env_ok","climate_ok",
              "photo_attempts","strict9_success","accepted_species_n"]].copy()
    env["accepted_species_n"]=env.accepted_species_n.fillna(0).astype(int)
    base_strict=set(pd.concat([s.target_id for _,s in existing_year_frames.values()]))
    b=source_photo.merge(env.rename(columns={"land_class":"existing_land_class"}),
                         on=["target_id","year"],how="left",validate="many_to_one")
    panel.require(b.loc[b.existing_land_class.notna(),"land_class"].eq(
        b.loc[b.existing_land_class.notna(),"existing_land_class"]).all(),
        "New photo coordinate assigned to different landmass at existing source site")
    accepted=b.wcvp_accepted_plant_name_id.notna()
    b["has_existing_exposure"]=(b.env_ok.fillna(False)&b.climate_ok.fillna(False))
    b["real_strict_accepted"]=(b.real_strict9&accepted)
    original_classified=set(b.loc[b.has_existing_exposure,"target_id"])
    island=b.land_class.eq("island")
    by_year_design=b.groupby(["year","land_class","source_design"],dropna=False).agg(
       independent_original_new_observation_photos=("unit_id","size"),
       true_strict9_original_photos=("real_strict9","sum"),
       strict9_wcvp_accepted_photos=("real_strict_accepted","sum"),
       real_photos_with_existing_exposure=("has_existing_exposure","sum")).reset_index()
    # Source-only environment opportunity audit: count unique REAL photo cells,
    # not unmeasured photo IDs as independent ecological outcomes.
    q=b.loc[b.real_strict_accepted].copy()
    existing=q.loc[q.has_existing_exposure]
    missing=q.loc[~q.has_existing_exposure]
    def group_stats(frame):
        return {
            "photos":int(len(frame)),
            "distinct_observations":int(frame.observation_id.nunique()),
            "distinct_original_photo_cell_years":int(frame.target_id.nunique()),
            "distinct_GSHHG_islands":int(frame.loc[frame.land_class.eq("island"),"land_id"].nunique()),
            "island_photo_records":int(frame.land_class.eq("island").sum()),
            "island_photo_cell_years":int(frame.loc[frame.land_class.eq("island"),"target_id"].nunique()),
        }
    names=set(orig.target_id)
    report={
       "signed_baseline_strict_5photo_3species_photo_cell_years":7982,
       "signed_baseline_strict_island_cell_years":710,
       "signed_baseline_broad_3photo_2species_photo_cell_years":20543,
       "signed_baseline_broad_island_cell_years":2026,
       "new_independent_observation_photos_with_landmass_and_selected_year":int(len(b)),
       "new_independent_observation_photos_true_strict9":int(b.real_strict9.sum()),
       "new_independent_observation_photos_true_strict9_and_WCVP":int(q.shape[0]),
       "new_observation_source_all_supplementary_arms_not_unbiased_random_sample":True,
       "strict9_WCVP_new_observations_with_existing_source_ALAN_exposure":group_stats(existing),
       "strict9_WCVP_new_observations_requiring_new_ALAN_exposure_or_other_env":group_stats(missing),
       "new_strict9_accepted_photos_on_existing_frozen_strict_sites":int(
           existing.target_id.isin(base_strict).sum()),
       "new_strict9_accepted_photo_cells_on_sites_below_frozen_strict_gate":int(
           existing.loc[~existing.target_id.isin(base_strict),"target_id"].nunique()),
       "new_true_strict9_WCVP_observation_ids_are_unique":bool(b.observation_id.is_unique),
       "photo_geography_independent_point_join_to_original_GSHHG":True,
       "no_supplementary_photo_retreated_as_original_historical_main":True,
       "no_new_pixels_downloaded_or_roi_inferred":True,
       "no_photo_colour_imputation":True,
    }
    return report,by_year_design

def run(paths,source_file,coast,out):
    out=Path(out);panel.require(not out.exists(),"Refuse overwrite of prior ALAN intake")
    original,wcvp=panel.load_frozen_photos(paths)
    source,whole=photo_source(source_file,original)
    with_wcvp=source_wcvp_accepted(source,wcvp)
    geocoded=geography_with_actual_landmasses(with_wcvp,coast)
    years={}
    for year in (2019,2021,2023):
        years[year]=panel.year_historic(year,paths,original,wcvp,return_all=True)
    years[2024]=panel.year2024(paths,original,wcvp,return_all=True)
    years[2025]=panel.year2025(paths,original,wcvp,return_all=True)
    supported,by_design=summarize(geocoded,years)
    report={
      "schema":"ALAN_extra_original1181146_georestored_v2_sameyear_island_intake_audit_v1",
      "extra_source_photo_asset_sha256":RESTORED_SHA,
      "canonical_photo_source_sha256":panel.PHOTO_SHA,
      "wcvp_source_sha256":panel.WCVP_SHA,
      "source_photo_years":list(YEARS),
      "original_world":whole,
      "ALAN_additional_true_photo_opportunities":supported,
      "caution":"Only real original measured pixels/accepted-taxonomy and georeferenced independent photos; new source arms had different inclusion probabilities, no main model refit or wild-native flower prevalence claim.",
      "all_source_input_SHA256":{**{k:panel.sha(v) for k,v in paths.items()},
                                 "extra_original_photo_archive":panel.sha(source_file),
                                 "coast_original_archive_sha256":COAST_SHA},
    }
    out.mkdir(parents=True)
    by_design.to_csv(out/"extra_photo_source_arm_year_landclass_aggregate.csv",index=False)
    report["aggregate_csv_sha256"]=panel.sha(out/"extra_photo_source_arm_year_landclass_aggregate.csv")
    (out/"additional_measured_photo_opportunities_summary.json").write_text(
         json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_EXTRA_REAL_1181146_SOURCE_PHOTO_INTAKE",
          json.dumps({"whole":whole,"support":supported},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp",
          "geo2019","viirs2019","pop2019",
          "geo2021","viirs2021","pop2021",
          "geo2023","viirs2023","pop2023",
          "geo2024","cn2024","uc2024","uc2024receipt",
          "geo2025","geo2025receipt","viirs2025")
    for key in keys:
        p.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    p.add_argument("--source-parquet",type=Path,required=True)
    p.add_argument("--coast",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in keys},a.source_parquet,a.coast,a.out)
