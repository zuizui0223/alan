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
    # Year loaders have two official population column conventions. The
    # concatenated panel contains BOTH names with missing values by year;
    # coalesce at the row level, never select one global column.
    statuses=(orig.get("population_status",pd.Series(index=orig.index,dtype=object))
              .combine_first(orig.get("population_source_status",
                         pd.Series(index=orig.index,dtype=object))))
    density=(orig.get("population_density_per_km2",
                      pd.Series(index=orig.index,dtype=float))
             .combine_first(orig.get("population_density",
                         pd.Series(index=orig.index,dtype=float))))
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


def candidate_observation_rescues(extra,canonical):
    """Actual extra photo success of a canonical *unknown* observation.

    Each known observation UUID is one ecological event, even when it has
    multiple genuine photos. Preserve first-photo colour and photo identities;
    the extra colour is a separately measured sibling replacement candidate.
    """
    original=canonical[["observation_id","photo_id","inat_taxon_id","source_subarm",
                        "observed_on","latitude","longitude","strict9","strict9_white"]].copy()
    original=original.rename(columns={
        "photo_id":"canonical_photo_id",
        "inat_taxon_id":"canonical_inat_taxon_id",
        "strict9":"canonical_strict9",
        "strict9_white":"canonical_strict9_white",
        "observed_on":"canonical_observed_on",
        "latitude":"canonical_latitude","longitude":"canonical_longitude"
    })
    orig_ids=set(original.observation_id)
    siblings=extra.loc[extra.real_strict9 &
                       extra.observation_id.isin(orig_ids)].copy()
    siblings=siblings.merge(original,on="observation_id",
                            how="inner",validate="many_to_one")
    panel.require(siblings.photo_id.ne(siblings.canonical_photo_id).all(),
                  "Same photo counted again as a new sibling colour")
    matched=siblings.inat_taxon_id.eq(siblings.canonical_inat_taxon_id)
    unknown=~siblings.canonical_strict9.astype(bool) | (
        siblings.canonical_strict9_white.isna())
    report={
       "extra_genuine_strict9_sibling_photos":len(siblings),
       "sibling_photo_records_with_conflicting_source_taxon_id":int((~matched).sum()),
       "strict9_extra_photos_for_canonical_already_classifiable_observations":int(
           (matched & ~unknown).sum()),
       "strict9_extra_photos_for_canonical_technical_unknowns":int(
           (matched & unknown).sum()),
       "technical_unknown_observation_UUIDs_with_1plus_measured_sibling_photos":int(
           siblings.loc[matched & unknown,"observation_id"].nunique()),
       "same_observation_UUID_siblings_are_not_new_independent_observations":True,
    }
    subset=siblings.loc[matched & unknown].copy()
    if subset.empty:
        return pd.DataFrame(columns=[
            "observation_id","source_subarm","latitude","longitude",
            "observed_on","wcvp_taxon_id","rescued_white_fraction",
            "extra_strict9_photos_on_same_observation"]),report
    group=subset.groupby("observation_id",observed=True).agg(
        source_subarm=("source_subarm","first"),
        latitude=("canonical_latitude","first"),
        longitude=("canonical_longitude","first"),
        observed_on=("canonical_observed_on","first"),
        wcvp_taxon_id=("canonical_inat_taxon_id","first"),
        rescued_white_fraction=("strict9_white","mean"),
        extra_strict9_photos_on_same_observation=("photo_id","nunique")
    ).reset_index()
    panel.require(group.observation_id.is_unique and
                  group.rescued_white_fraction.between(0,1).all(),
                  "Original photo rescue must be observation-unique with measured white")
    report["unique_canonical_observations_rescuable_original_roi9"]=len(group)
    report["historical_main_observations_rescuable_original_roi9"]=int(
        group.source_subarm.eq("historical_main").sum())
    return group,report

def map_rescues_to_existing_ALAN_year_exposure(rescues,wcvp,years,coast):
    """Only original canonical dated & located observations; no new species/geo claims."""
    if rescues.empty:
        return {"source_years":list(YEARS),"rescuable_with_year_and_valid_source_coordinates":0,
                "rescuable_historical_main_island_accepted_with_existing_ALAN_exposure":0,
                "rescuable_historical_main_island_photocell_years_with_exposure":0,
                "rescuable_historical_main_distinct_islands_with_exposure":0,
                "all_missing_exposure_remain_unlabelled":True},pd.DataFrame(),pd.DataFrame()
    ref=wcvp.loc[wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
        ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    ref.inat_taxon_id=ref.inat_taxon_id.astype(str)
    q=rescues.copy()
    q["year"]=pd.to_numeric(q.observed_on.fillna("").astype(str).str[:4],
                             errors="coerce")
    q["latitude"]=pd.to_numeric(q.latitude,errors="coerce")
    q["longitude"]=pd.to_numeric(q.longitude,errors="coerce")
    q=q.loc[q.year.isin(YEARS)&q.latitude.between(-90,90)&
            q.longitude.between(-180,180)].copy()
    q=q.rename(columns={"wcvp_taxon_id":"inat_taxon_id"})
    q.inat_taxon_id=q.inat_taxon_id.astype(str)
    q=q.merge(ref,on="inat_taxon_id",how="left",validate="many_to_one")
    src=Path(coast)
    land=prepare_land(gpd.read_file(src/"GSHHS_h_L1.shp"))
    lakes=gpd.read_file(src/"GSHHS_h_L2.shp")
    points=q.rename(columns={"observation_id":"unit_id"})
    placed=assign_points(points,land,lakes)
    q=points.merge(placed[["unit_id","land_id","land_class"]],
                    on="unit_id",validate="one_to_one")
    q=q.loc[q.land_id.notna()&q.land_class.isin(STATES)].copy()
    q["cell_id"]=cell_id(q.longitude,q.latitude)
    q["land_id"]=q.land_id.astype(str)
    q["target_id"]=q.year.astype(int).astype(str)+"|"+q.cell_id+"|"+q.land_id
    source=[]
    for year,(full,_) in years.items():
        e=full.copy()
        status=e["population_status"] if "population_status" in e.columns else e["population_source_status"]
        pop=e["population_density_per_km2"] if "population_density_per_km2" in e.columns else e["population_density"]
        e["env_ok"]=(e.viirs_valid_months.ge(3)&e.viirs_median.notna()&
                     status.eq("OK")&pop.notna()&
                     e[["bio5","bio12","elevation"]].notna().all(axis=1))
        e["year"]=int(year)
        source.append(e[["year","target_id","env_ok"]])
    known=pd.concat(source,ignore_index=True)
    panel.require(known.target_id.is_unique,
                  "Original five-year ALAN environment has duplicate site IDs")
    q=q.merge(known,on=["year","target_id"],how="left",
              validate="many_to_one")
    qualified=q.loc[q.source_subarm.eq("historical_main")&
         q.wcvp_accepted_plant_name_id.notna()&q.env_ok.fillna(False)&
         q.land_class.eq("island")]
    by_source=q.groupby(["year","land_class","source_subarm"],observed=True).agg(
        rescued_obs=("unit_id","nunique"),
        restored_hue_photos=("extra_strict9_photos_on_same_observation","sum"),
        wcvp_accepted_rescues=("wcvp_accepted_plant_name_id",
                                lambda s:int(s.notna().sum())),
        source_photo_cells=("target_id","nunique")).reset_index()
    report={
        "source_years":list(YEARS),
        "rescuable_with_year_and_valid_source_coordinates":len(q),
        "rescuable_with_actual_WCVP_accepted_taxon":int(
             q.wcvp_accepted_plant_name_id.notna().sum()),
        "rescuable_historical_main_island_accepted_with_existing_ALAN_exposure":len(qualified),
        "rescuable_historical_main_island_photocell_years_with_exposure":int(
             qualified.target_id.nunique()),
        "rescuable_historical_main_distinct_islands_with_exposure":int(
             qualified.land_id.nunique()),
        "all_missing_exposure_remain_unlabelled":True,
        "observation_UUIDs_not_photo_IDs_are_counted_as_newly_classifiable":True
    }
    return report,by_source,q


def admit_real_sibling_roi_to_photo_site_thresholds(rescues,canonical,wcvp,years):
    """Exact historical-main site threshold gains from genuinely measured sibling ROI.

    The source observation already existed. Only previously unclassifiable
    canonical observation-level photo ROI labels are replaced by true sibling
    original-photo ROI measurements. Do not increase photo/observation attempts.
    Changes are source-supported eligibility counts, NOT refitted ALAN slopes.
    """
    whole=[]
    for year,(full,_) in years.items():
        f=full.copy()
        f["year"]=int(year)
        f["cell_id"]=f.target_id.astype(str).str.split("|",regex=False).str[1]
        p=(f["population_status"] if "population_status" in f.columns
           else f["population_source_status"])
        density=(f["population_density_per_km2"] if "population_density_per_km2" in f.columns
                 else f["population_density"])
        f["eligibility_environment"]=(f.viirs_valid_months.ge(3)&
            f.viirs_median.notna()&p.eq("OK")&density.notna()&
            f[["bio5","bio12","elevation"]].notna().all(axis=1))
        whole.append(f[["target_id","year","cell_id","land_id","land_class",
                        "photo_attempts","strict9_success","accepted_species_n",
                        "eligibility_environment"]])
    old=pd.concat(whole,ignore_index=True)
    panel.require(old.target_id.is_unique,"Repeated signed source cell year")
    old["accepted_species_n"]=old.accepted_species_n.fillna(0).astype(int)
    # All true original photo acceptance conditions must reproduce existing
    # strict and broad source denominators before admitting any sibling photo.
    def eligible(frame,a,s,t):
        return (frame.eligibility_environment&
                frame.photo_attempts.ge(a)&
                frame.strict9_success.ge(s)&
                frame.accepted_species_n.ge(t))
    base={}
    for name,a,s,t in (("strict_5photos_3accepted",5,3,3),
                       ("broad_3photos_2accepted",3,2,2)):
        flags=eligible(old,a,s,t)
        baseline=old.loc[flags]
        base[name]={"cell_years":len(baseline),
                    "island_cell_years":int(baseline.land_class.eq("island").sum())}
    panel.require(base["strict_5photos_3accepted"]=={
                    "cell_years":7982,"island_cell_years":710} and
                  base["broad_3photos_2accepted"]=={
                    "cell_years":20543,"island_cell_years":2026},
                  "Original strict/broad photo-cell-year eligibility drift")

    accepted_lookup=wcvp.loc[
       wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
       ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    accepted_lookup.inat_taxon_id=accepted_lookup.inat_taxon_id.astype(str)
    # Reconstruct already observed accepted taxa on original historical-main
    # strict ROI photos for EXACT source cells; do not use unknown source taxa.
    c=canonical.loc[
        canonical.source_subarm.eq("historical_main")&
        canonical.strict9.astype(bool)&canonical.strict9_white.notna()&
        canonical.latitude.between(-90,90)&canonical.longitude.between(-180,180),
        ["observation_id","inat_taxon_id","observed_on","latitude","longitude"]].copy()
    c["year"]=pd.to_numeric(c.observed_on.astype(str).str[:4],errors="coerce")
    c=c.loc[c.year.isin(YEARS)].copy()
    c["cell_id"]=cell_id(c.longitude,c.latitude)
    # Ambiguous geographic cells are excluded from assignment rather than
    # silently reassigned to a nearby island or the continent.
    ambiguity=old.duplicated(["year","cell_id"],keep=False)
    unique_geo=old.loc[~ambiguity,["year","cell_id","target_id"]]
    c=c.merge(accepted_lookup,on="inat_taxon_id",how="left",validate="many_to_one")
    c=c.loc[c.wcvp_accepted_plant_name_id.notna()]
    c=c.merge(unique_geo,on=["year","cell_id"],how="inner",validate="many_to_one")
    source_pairs=c[["target_id","wcvp_accepted_plant_name_id"]].drop_duplicates()
    pre_n=source_pairs.groupby("target_id").size()
    baseline_old_n=old.set_index("target_id").loc[pre_n.index,"accepted_species_n"]
    mismatched=(pre_n.astype(int).ne(baseline_old_n.astype(int)))
    panel.require(not mismatched.any(),
          "Original accepted-taxon exact site reconstruction does not match frozen counts")
    # All source rescued observation IDs still existed, so original opportunity
    # count never changes. Only authentic photo-colour successes change.
    q=rescues.loc[rescues.source_subarm.eq("historical_main")&
                  rescues.target_id.isin(set(old.target_id))].copy()
    panel.require(not q.unit_id.duplicated().any(),
                 "One original observation counted twice as an ROI rescue")
    added=q.groupby("target_id").unit_id.nunique()
    new_pairs=q.loc[q.wcvp_accepted_plant_name_id.notna(),
         ["target_id","wcvp_accepted_plant_name_id"]].drop_duplicates()
    if len(new_pairs):
        tagged=new_pairs.merge(source_pairs.assign(_existing_species=1),
                  on=["target_id","wcvp_accepted_plant_name_id"],how="left",
                  validate="one_to_one")
        incremental=tagged.loc[tagged._existing_species.isna()]
        new_species_count=incremental.groupby("target_id").size()
    else:
        new_species_count=pd.Series(dtype=int)
    post=old.copy()
    post["new_genuine_roi9_observation_measurements"]=post.target_id.map(added).fillna(0).astype(int)
    post["new_distinct_accepted_species"]=post.target_id.map(new_species_count).fillna(0).astype(int)
    post["strict9_success"]=post.strict9_success+post.new_genuine_roi9_observation_measurements
    post["accepted_species_n"]=post.accepted_species_n+post.new_distinct_accepted_species
    result={}
    for name,a,s,t in (("strict_5photos_3accepted",5,3,3),
                       ("broad_3photos_2accepted",3,2,2)):
        b=eligible(old,a,s,t)
        e=eligible(post,a,s,t)
        panel.require(not (b&~e).any(),
           "Authentic sibling ROI colour rescue removed a previously eligible cell")
        newcomer=post.loc[e&~b]
        result[name]={
            "previous_cell_years":int(b.sum()),
            "source_rescued_eligible_cell_years":int(e.sum()),
            "additional_actual_photo_colour_eligible_cell_years":len(newcomer),
            "previous_island_cell_years":int((b&old.land_class.eq("island")).sum()),
            "source_rescued_island_cell_years":int((e&post.land_class.eq("island")).sum()),
            "additional_island_cell_years":int(newcomer.land_class.eq("island").sum()),
            "different_island_landmasses_adding_eligible_cells":int(
                 newcomer.loc[newcomer.land_class.eq("island"),"land_id"].nunique()),
        }
    return {
       "schema":"ALAN_signed_real_sibling_photo_ROI_observation_rescue_eligibility_v1",
       "new_genuine_strict9_observation_measurements_with_original_source_site":int(
            post.new_genuine_roi9_observation_measurements.sum()),
       "additional_accepted_species_per_photo_cell_sum":int(
            post.new_distinct_accepted_species.sum()),
       "no_new_distinct_observation_uuids_or_photo_site_claim":True,
       "no_new_photo_opportunity_included_in_eligibility_denominator":True,
       "same_original_historical_main_and_sameyear_environment_source":True,
       "technical_unknown_first_photo_not_relabelled_itself":True,
       "source_sibling_original_photo_strict9_is_actual_measured_colour":True,
       "new_eligible_uses_same_success_and_species_thresholds":True,
       "results":result
    }

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
    rescues,rescue_counts=candidate_observation_rescues(source,original)
    rescue_geography,rescue_by_design,rescue_photo_cells=map_rescues_to_existing_ALAN_year_exposure(
        rescues,wcvp,years,coast)
    rescued_eligibility=admit_real_sibling_roi_to_photo_site_thresholds(
        rescue_photo_cells,original,wcvp,years)
    report={
      "schema":"ALAN_extra_original1181146_georestored_v2_sameyear_island_intake_audit_v1",
      "extra_source_photo_asset_sha256":RESTORED_SHA,
      "canonical_photo_source_sha256":panel.PHOTO_SHA,
      "wcvp_source_sha256":panel.WCVP_SHA,
      "source_photo_years":list(YEARS),
      "original_world":whole,
      "ALAN_additional_true_photo_opportunities":supported,
      "original_observation_multi_photo_colour_rescue":rescue_counts,
      "rescuable_original_observation_geographic_ALAN_exposure":rescue_geography,
      "source_rescued_original_photo_eligibility":rescued_eligibility,
      "caution":"Only real original measured pixels/accepted-taxonomy and georeferenced independent photos; new source arms had different inclusion probabilities, no main model refit or wild-native flower prevalence claim.",
      "all_source_input_SHA256":{**{k:panel.sha(v) for k,v in paths.items()},
                                 "extra_original_photo_archive":panel.sha(source_file),
                                 "coast_original_archive_sha256":COAST_SHA},
    }
    out.mkdir(parents=True)
    by_design.to_csv(out/"extra_photo_source_arm_year_landclass_aggregate.csv",index=False)
    rescue_by_design.to_csv(out/"rescuable_sibling_photo_original_observations_aggregate.csv",index=False)
    report["aggregate_csv_sha256"]=panel.sha(out/"extra_photo_source_arm_year_landclass_aggregate.csv")
    (out/"additional_measured_photo_opportunities_summary.json").write_text(
         json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_EXTRA_REAL_1181146_SOURCE_PHOTO_INTAKE",
          json.dumps({"whole":whole,"support":supported,"rescue":rescue_counts,"rescue_environment":rescue_geography,"new_eligibility":rescued_eligibility},sort_keys=True),flush=True)
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
