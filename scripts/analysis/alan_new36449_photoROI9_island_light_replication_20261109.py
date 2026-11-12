#!/usr/bin/env python3
"""Independent ORIGINAL-photo island area × nighttime light source replication.

This is a NEW, conditionally island/site-targeted sample, entirely disjoint in
observation UUID and photo ID from the original historical-main 2019-25 sample.
It does NOT have independent photo selection probability or validate native
flowering status. Never silently combine it as an unweighted global census.

Input: all 73 SHA-verified source groups already measured in original RGB:
 tier06 8942 new original photo UUIDs; tier07 27507 new original photo UUIDs.
Accepted WCVP plant source taxa, exact signed GSHHG cell×year, and actual FCP
ROI9 white mass are connected to the original matched same-year VIIRS, WorldPop,
bio5, rain, elevation. Compare the new-source strict quality/within-island
sampling support (>=5 original attempts, >=3 genuine ROI9 successes,
>=3 distinct accepted WCVP species actually measured) across
tier06 / tier07 / the exact source union.

Fit original ALAN-style island fixed effects, year fixed effects, source-island
cluster robust errors, same covariates and separate area/isolation moderators.
This tests photographic association transfer across DIFFERENT acquisition
designs, NOT causal effects of nightlights on floral pigment, species loss or
pollinator health. Never hide failed numeric or source joins.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import alan_union_36449_original_island_photoROI9_20261109 as union
from scripts.analysis import alan_island_area_light_20261108 as area
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as historic_tier
from scripts.analysis import alan_isolation_area_light_20261108 as isolation
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import cell_id

SCHEMA="ALAN_36449_new_source_real_white_ROI_island_area_isolation_ALAN_replication_v1"
SOURCE_EXPECTED={"tier06":(8942,3597,156),"tier07":(27507,10245,769),
                 "union":(36449,13842,1137)}
YEARS=(2019,2021,2023,2024,2025)
PALETTE=union.PALETTE
WHITE_COL="strict9_white"
FRACTIONS=tuple("flower_fraction_"+p for p in PALETTE)

def require(ok,message):
    if not bool(ok):raise ValueError("ALAN_NEW36449_COLOUR_LIGHT_REPLICATION: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda:f.read(4<<20),b""):h.update(buf)
    return h.hexdigest()

def original_colour_from_file(path):
    """Read original measured white *fraction*; all true nine palette gates intact."""
    path=Path(path)
    base=union.read_numeric_photo(path)
    raw=pd.read_csv(path,compression="gzip" if path.suffix==".gz" else None,
                    dtype=str,keep_default_na=False)
    require(base.photo_id.equals(raw.photo_id) and
            base.observation_id.equals(raw.observation_id),
            "Original numeric ROI validation and exact photographed source not aligned")
    if WHITE_COL in raw.columns and "strict9" in raw.columns:
        frac=pd.to_numeric(raw[WHITE_COL],errors="coerce")
    else:
        require(set(FRACTIONS)<=set(raw.columns),
                "No genuine original flower-region nine-colour numeric masses")
        fractions=raw.loc[:,list(FRACTIONS)].apply(
            pd.to_numeric,errors="coerce").to_numpy(dtype=float)
        valid=base.real_numeric_ROI9.to_numpy(bool)
        sums=np.where(np.isfinite(fractions),fractions,0).sum(axis=1)
        frac=pd.Series(np.where(valid,fractions[:,0]/np.maximum(sums,1e-100),
                                np.nan),index=raw.index)
    valid=base.real_numeric_ROI9.to_numpy(bool)
    require(frac.loc[valid].between(0,1).all() and
            frac.loc[~valid].isna().all(),
            "Technical unknown original image was invented as white or nonwhite")
    base["source_photo_white_fraction"]=frac.to_numpy()
    return base

def true_source_colour(prior06,bulk06,pilot07,bulk07):
    """Use pre-existing signed original ROI9 measurements, never run FCP again."""
    # This also enforces SHA, FCP method blob, 73 independent source receipts.
    p6,c1=union.read_prior06(prior06)
    b6,c2=union.read_bulk06(bulk06)
    p7,c3=union.read_pilot07(pilot07)
    b7,c4=union.read_bulk07(bulk07)
    authorised=pd.concat([*p6,*b6,p7,*b7],ignore_index=True)
    require(len(authorised)==36449 and authorised.photo_id.is_unique and
            authorised.observation_id.is_unique and
            int(authorised.real_numeric_ROI9.sum())==13842 and
            len(c1)+len(c2)+1+len(c4)==73,
            "Original signed 73-image-result provenance changed")
    paths=[]
    for batch in range(1,6):
        if batch==1:
            p=Path(prior06)/"first"/"finished"/"new256_independent_island_actual_original_ROI9.csv.gz"
        else:
            p=Path(prior06)/f"alan-original-novel-island-ROI9-batch-{batch}"/"finished"/f"ALAN_original_novel_island_photo_ROI9_batch_{batch:02d}.csv.gz"
        paths.append(p)
    for group in range(15):
        for shard in range(4):
            paths.append(Path(bulk06)/f"alan-real-originalROI9-7662-group-{group}"/"measured"/f"shard_{shard:03d}"/"photos.csv")
    paths.append(Path(pilot07)/"finished"/"ALAN_tier07_first1024_original_floral_roi9.csv.gz")
    for task in range(52):
        paths.append(Path(bulk07)/f"alan-tier07-real-originalROI9-26483-task-{task}"/"finished"/f"alan_tier07_original_photo_ROI9_task{task:02d}.csv.gz")
    actual=pd.concat([original_colour_from_file(p) for p in paths],
                     ignore_index=True)
    require(len(actual)==len(authorised) and actual.photo_id.is_unique and
            actual.observation_id.is_unique and
            set(actual.photo_id)==set(authorised.photo_id),
            "True original 9colour photo fractions do not match complete signed image IDs")
    check=actual.merge(authorised,on=["photo_id","observation_id"],validate="one_to_one",
             suffixes=("","_signed"))
    require(len(check)==36449 and
            check.real_numeric_ROI9.eq(check.real_numeric_ROI9_signed).all(),
            "True original photographic colour validation differs from SHA-signed source")
    return actual,{"source_fcp_model_git_blob":union.ROI_BLOB,
                  "actual_original_measurement_files":len(paths),
                  "original_source_signed_work_results":73,
                  "source_photo_UUIDs":len(actual),
                  "true_numeric_ROI9_success":int(actual.real_numeric_ROI9.sum()),
                  "actual_missing_technical_ROI9":int((~actual.real_numeric_ROI9).sum())}

def photo_cell_species_equal(frame):
    """Count original photo opportunities, then measure means by accepted taxa."""
    require(frame.photo_id.is_unique and frame.observation_id.is_unique,
            "Original photos/UUIDs duplicated within source-only target")
    trial=frame.groupby("target_id",observed=True).agg(
        photo_attempts=("observation_id","nunique"),
        source_species_n=("wcvp_accepted_plant_name_id","nunique"),
        strict9_success=("real_numeric_ROI9","sum"),
        land_id=("land_id","first"),
        year=("year","first")).reset_index()
    valid=frame.loc[frame.real_numeric_ROI9].copy()
    measured_taxon=valid.groupby(["target_id","wcvp_accepted_plant_name_id"],
        observed=True).source_photo_white_fraction.mean().reset_index()
    across=measured_taxon.groupby("target_id",observed=True).agg(
        accepted_species_n=("wcvp_accepted_plant_name_id","nunique"),
        white_equal=("source_photo_white_fraction","mean")
        ).reset_index()
    trial=trial.merge(across,on="target_id",how="left",validate="one_to_one")
    trial["accepted_species_n"]=trial.accepted_species_n.fillna(0).astype(int)
    trial["real_strict_eligible"]=(trial.photo_attempts.ge(5)&
        trial.strict9_success.ge(3)&trial.accepted_species_n.ge(3))
    eligible=trial.loc[trial.real_strict_eligible].copy()
    eligible["white_pp"]=100*eligible.white_equal
    require(eligible.white_pp.between(0,100).all() and
            eligible.strict9_success.ge(3).all() and
            eligible.accepted_species_n.ge(3).all(),
            "Strict photo site has invalid FCP ROI9 white or real measured species")
    return trial,eligible

def original_sameyear_environment(paths):
    """Original signed env matched by exact year × cell × island GSHHG identity."""
    canonical,wcvp=panel.load_frozen_photos(paths)
    tables={}
    for year in YEARS:
        if year in (2019,2021,2023):
            full,_=panel.year_historic(year,paths,canonical,wcvp,return_all=True)
        else:
            full,_=getattr(panel,f"year{year}")(
                paths,canonical,wcvp,return_all=True)
        require(full.target_id.is_unique and
                full.land_class.isin(("island","mainland")).all(),
                "Original matched same-year photographic VIIRS/WorldPop cells changed")
        q=full.copy()
        status=q["population_status"] if "population_status" in q.columns else q["population_source_status"]
        density=q["population_density_per_km2"] if "population_density_per_km2" in q.columns else q["population_density"]
        q["environment_complete"]=(q.viirs_valid_months.ge(3)&q.viirs_median.notna()&
            status.eq("OK")&density.notna()&
            q[["bio5","bio12","elevation"]].notna().all(axis=1))
        q["log_light"]=np.log1p(q.viirs_median.clip(lower=0))
        q["log_pop"]=np.log1p(density.clip(lower=0))
        q["log_rain"]=np.log1p(q.bio12.clip(lower=0))
        q["year"]=int(year)
        tables[year]=q[["target_id","land_id","land_class","year",
                        "environment_complete","log_light","log_pop",
                        "bio5","log_rain","elevation"]].copy()
    source=pd.concat(tables.values(),ignore_index=True)
    require(source.target_id.is_unique and
            set(source.year)==set(YEARS) and
            source.log_light.notna().sum()>10000,
            "Original time-indexed exposure not independently recoverable")
    source.land_id=source.land_id.astype(str)
    return source

def model_supported(cohort,environment,land,distances):
    """Fit joint area/isolation on genuinely repeated *different* sites."""
    dsource=cohort.copy()
    dsource["year"]=pd.to_numeric(dsource.year,errors="raise").astype(int)
    d= dsource.merge(environment,on=["target_id","land_id","year"],
                   how="left",validate="one_to_one",suffixes=("","_env"))
    missing=int(d.environment_complete.isna().sum())
    invalid=int(d.environment_complete.fillna(False).eq(False).sum())
    require(not missing and not invalid and d.land_class.eq("island").all(),
            "Genuine strict source islands lack ORIGINAL-year VIIRS/WorldPop/climate exposure")
    require(d[["white_pp","log_light","log_pop","bio5","log_rain","elevation"]].notna().all().all(),
            "Original adjusted regression contains nonfinite source covariates")
    d["is_island"]=d.land_class.eq("island").astype(int)
    d["log_attempts"]=np.log1p(d.photo_attempts.astype(float))
    islands,selected,g=area.join_area_and_select(d,land)
    support={"strict_photo_cell_years":len(d),
             "strict_distinct_GSHHG_islands":int(islands.land_id.nunique()),
             "spatially_repeated_islands":int(selected.land_id.nunique()),
             "spatially_repeated_source_cell_years":len(selected),
             "largest_island_share_of_spatial_sample":(
                  float(selected.land_id.value_counts().iloc[0]/len(selected))
                  if len(selected) else None),
             "islands_have_genuine_different_geographic_0p05deg_sites":True,
             "source_pop_and_VIIRS_same_original_photo_year":True}
    if selected.land_id.nunique()<10 or len(selected)<80:
        return {"estimable":False,"reason":"Insufficient independent multi-site source GSHHG islands",
                "source_support":support}
    df=isolation.complete_island_frame(selected,distances)
    center_area=float(np.log10(df.groupby("land_id").area_km2.first()).median())
    center_dist=float(np.log10(1+df.groupby("land_id")
                     .distance_mainland_vertex_km.first()).median())
    fit={}
    for kind in ("joint",):
        for weight,eq in (("cell_equal",False),("island_equal",True)):
            fit[f"{kind}_{weight}"]=isolation.fit_light_moderation(
                df,center_area,center_dist,kind,eq)
    # Remove the one most photographed island to evaluate concentration,
    # retaining the other real source GSHHG cell×year rows.
    top=str(df.land_id.value_counts().index[0])
    trim=df.loc[df.land_id.ne(top)]
    if trim.land_id.nunique()>=10 and len(trim)>=80:
        fit["leave_top_island_out_equal_island"]=isolation.fit_light_moderation(
             trim,center_area,center_dist,"joint",True)
    support["source_centre_area_log10"]=center_area
    support["source_centre_mainland_distance_log10_1pkm"]=center_dist
    support["largest_island_omitted_for_diagnostic"]=True
    return {"estimable":True,"source_support":support,
            "models":fit,
            "post_selection_observational_exploratory_not_causal":True}

def same_site_taxon_photo_repeatability(canonical,wcvp,old_strict,new_strict,new_photos):
    """Actual old/new photographic source agreement at the SAME site and year.

    Different original photo/UUID samples can represent different photographed
    accepted taxa even at the same geographic site and year. No assumption
    that absence from volunteer photos means a plant is absent in the flora.
    """
    old=old_strict.copy()
    new=new_strict.copy()
    overlap=old[["target_id","land_id","white_pp","accepted_species_n"]].merge(
        new[["target_id","white_pp","accepted_species_n"]],
        on="target_id",how="inner",validate="one_to_one",
        suffixes=("_old","_new"))
    # Sharing a landmass in both source samples is NOT identical to sharing
    # an exact original site×year. The union source inventory independently
    # established 44 shared *landmass IDs*, but only a subset can have
    # exact shared strict source photo cell-years.
    require(len(overlap)==301 and
            overlap.land_id.nunique()>=1 and
            overlap.land_id.nunique()<=44,
            "Original historical/new source exact site-year overlap count drifted")
    id_parts=overlap.target_id.str.split("|",regex=False,expand=True)
    require(id_parts.shape[1]==3,"Original source cell year target malformed")
    targets=pd.DataFrame({
       "target_id":overlap.target_id,
       "year":id_parts[0].astype(int),
       "cell_id":id_parts[1],
    })
    require(not targets.duplicated(["year","cell_id"]).any(),
            "Old geographic site/year lands are mixed and cannot be assigned")
    p=canonical.loc[
       canonical.source_subarm.eq("historical_main")&
       canonical.strict9.astype(bool)&canonical.strict9_white.notna()&
       canonical.latitude.between(-90,90)&
       canonical.longitude.between(-180,180),
       ["photo_id","observation_id","inat_taxon_id","observed_on",
        "latitude","longitude","strict9_white"]].copy()
    p["year"]=pd.to_numeric(p.observed_on.astype(str).str[:4],errors="coerce")
    p=p.loc[p.year.isin(YEARS)].copy()
    p["cell_id"]=cell_id(p.longitude,p.latitude)
    p=p.merge(targets,on=["year","cell_id"],how="inner",validate="many_to_one")
    lookup=wcvp.loc[wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
                     ["inat_taxon_id","wcvp_accepted_plant_name_id"]]
    require(lookup.inat_taxon_id.is_unique,"Old original WCVP accepted taxon source changed")
    p=p.merge(lookup,on="inat_taxon_id",how="inner",validate="many_to_one")
    old_species=p.groupby(["target_id","wcvp_accepted_plant_name_id"],
          observed=True).strict9_white.mean().reset_index(name="old_photo_fraction")
    require(old_species.old_photo_fraction.between(0,1).all(),
            "Old photographed strict source white fraction invalid")
    expected=overlap.set_index("target_id").accepted_species_n_old
    observed=old_species.groupby("target_id").size().reindex(expected.index,fill_value=0)
    require(np.array_equal(observed.to_numpy(int),expected.to_numpy(int)),
            "Old strict site accepted taxa not reconstructed exactly from photographs")
    n=new_photos.loc[new_photos.real_numeric_ROI9 &
                      new_photos.target_id.isin(set(overlap.target_id))].copy()
    new_species=n.groupby(["target_id","wcvp_accepted_plant_name_id"],
          observed=True).source_photo_white_fraction.mean().reset_index(
             name="new_photo_fraction")
    observed_new=new_species.groupby("target_id").size().reindex(expected.index,fill_value=0)
    require(np.array_equal(observed_new.to_numpy(int),
                           overlap.set_index("target_id").accepted_species_n_new.to_numpy(int)),
            "New source successfully measured WCVP taxa do not match strict site")
    matched=old_species.merge(new_species,
         on=["target_id","wcvp_accepted_plant_name_id"],how="inner",
         validate="one_to_one")
    taxon_old={k:set(g.wcvp_accepted_plant_name_id) for k,g in old_species.groupby("target_id")}
    taxon_new={k:set(g.wcvp_accepted_plant_name_id) for k,g in new_species.groupby("target_id")}
    jj=[len(taxon_old[k]&taxon_new[k])/len(taxon_old[k]|taxon_new[k])
        for k in overlap.target_id]
    delta=overlap.white_pp_new-overlap.white_pp_old
    result={
      "shared_exact_original_site_and_observation_year_cells":len(overlap),
      "shared_original_GSHHG_islands":int(overlap.land_id.nunique()),
      "old_strict_original_photo_wcvp_taxon_site_pairs":len(old_species),
      "new_strict_original_photo_wcvp_taxon_site_pairs":len(new_species),
      "same_original_site_year_taxon_pairs_photographed_in_both_sources":len(matched),
      "same_site_years_with_at_least_one_accepted_taxon_photographed_in_both_sources":int(
           matched.target_id.nunique()),
      "same_source_same_site_taxon_jaccard_median":float(np.median(jj)),
      "same_source_same_site_taxon_jaccard_mean":float(np.mean(jj)),
      "same_site_years_with_zero_overlap_in_actually_photographed_taxa":int(
           sum(v==0 for v in jj)),
      "photographic_site_whiteness_old_new_pearson_r":float(
           overlap.white_pp_old.corr(overlap.white_pp_new)),
      "mean_site_new_minus_old_white_ROI_pp":float(delta.mean()),
      "median_abs_site_old_new_white_ROI_difference_pp":float(delta.abs().median()),
      "p90_abs_site_old_new_white_ROI_difference_pp":float(delta.abs().quantile(.9)),
      "matched_same_taxon_photo_fraction_median_abs_difference_pp":(
          float((100*(matched.new_photo_fraction-matched.old_photo_fraction)).abs().median())
          if len(matched) else None),
      "actual_new_and_old_observation_UUIDs_disjoint_but_same_island_environment":True,
      "photos_different_species_presence_not_a_field_absence_census":True,
      "matched_taxon_same_year_photo_does_not_prove_same_genotype_or_pigment_adaptation":True,
      "same_cell_photo_source_comparison_is_descriptive_not_randomized":True
    }
    require(np.isfinite([
       result["same_source_same_site_taxon_jaccard_mean"],
       result["photographic_site_whiteness_old_new_pearson_r"],
       result["median_abs_site_old_new_white_ROI_difference_pp"]]).all(),
       "Independent exact site-year source photographic comparison is unstable")
    return result

def run(paths,coast,source06,source07,prior06,bulk06,pilot07,bulk07,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwrite of new source ALAN ecological result")
    first=union.source_rows(source06,union.SHA06,"tier06",8942)
    last=union.source_rows(source07,union.SHA07,"tier07",27507)
    require(not set(first.photo_id)&set(last.photo_id) and
            not set(first.observation_id)&set(last.observation_id),
            "New source 93m photo tiers overlap observation UUID")
    numeric,proof=true_source_colour(prior06,bulk06,pilot07,bulk07)
    ref=pd.concat([first,last],ignore_index=True)
    allphotos=ref.merge(numeric,on=["photo_id","observation_id"],
                        how="inner",validate="one_to_one")
    require(len(allphotos)==36449 and allphotos.photo_id.is_unique and
            allphotos.observation_id.is_unique,
            "True real numeric photo mass not one-to-one with precolour original source")
    groups={}
    for key,n in SOURCE_EXPECTED.items():
        sub=allphotos if key=="union" else allphotos.loc[
            allphotos.original_source_tier.eq(key)]
        total,selected=photo_cell_species_equal(sub)
        require(len(sub)==n[0] and int(sub.real_numeric_ROI9.sum())==n[1] and
                len(selected)==n[2],
                "Previously independently signed original source photo/strict cell numbers drift: "+key)
        groups[key]=(total,selected)
    environment=original_sameyear_environment(paths)
    original_land=isolation.load_signed_land(coast)
    ids=set().union(*(
       set(d.land_id.astype(str)) for _,d in groups.values()))
    # Exact same GSHHG source island polygon area and closest MAINLAND coast
    # VERTEX distances for every island in the comparison (not centroids).
    distances=isolation.isolated_mainland_distance_km(original_land,ids)
    landmeta=original_land[["land_id","land_class","area_km2"]]
    # Measure true target-population overlap instead of asserting that new
    # original-photo IDs imply independent GSHHG island environments.
    older=historic_tier.by_year_tier(
        historic_tier.build_frozen_five_years(paths))["strict_5p3r3s"]
    older=older.loc[older.is_island.eq(1)].copy()
    require(len(older)==710 and older.target_id.is_unique and
            int(older.land_id.nunique())==101,
            "Original historical-main 710-cell 101-island photo support altered")
    fresh=groups["union"][1]
    old_sites=set(older.target_id)
    new_sites=set(fresh.target_id)
    canonical,wcvp=panel.load_frozen_photos(paths)
    paired_photography=same_site_taxon_photo_repeatability(
        canonical,wcvp,older,fresh,allphotos)
    old_lands=set(older.land_id.astype(str))
    new_lands=set(fresh.land_id.astype(str))
    shared_lands=old_lands&new_lands
    old_multi=area.join_area_and_select(older,landmeta)[1]
    overlap={
       "old_historical_main_strict_island_cells":710,
       "old_historical_main_island_landmasses":101,
       "old_historical_main_true_multisite_islands":int(old_multi.land_id.nunique()),
       "new_acquisition_source_strict_island_cells":len(fresh),
       "new_acquisition_source_strict_islands":len(new_lands),
       "strict_cell_years_present_in_both_acquisition_and_historical_main":len(old_sites&new_sites),
       "strict_distinct_GSHHG_islands_present_in_both_source_cohorts":len(shared_lands),
       "old_original_strict_cell_years_on_shared_islands":int(
            older.land_id.astype(str).isin(shared_lands).sum()),
       "new_original_strict_cell_years_on_shared_islands":int(
            fresh.land_id.astype(str).isin(shared_lands).sum()),
       "old_and_new_source_original_photo_observation_UUIDs_are_disjoint":True,
       "overlapping_sites_are_not_independent_environmental_replicates":True,
    }
    stats={}
    for key,(allcells,strict) in groups.items():
        stats[key]={
          "real_original_photo_observation_UUIDs":SOURCE_EXPECTED[key][0],
          "numeric_original_nine_palette_ROI9":SOURCE_EXPECTED[key][1],
          "actual_5attempt_3ROI9_3accepted_photo_cell_years":len(strict),
          "strict_photo_islands":int(strict.land_id.nunique()),
          "sameyear_exposure_join_and_area_isolation_model":
              model_supported(strict,environment,landmeta,distances)
        }
    unionModel=stats["union"]["sameyear_exposure_join_and_area_isolation_model"]
    # Explicitly track whether the new source can be interpreted as an
    # independent directional replication; do not force positive slope.
    if unionModel.get("estimable"):
        terms=unionModel["models"]["joint_island_equal"]["terms"]
        stats["union"]["independent_new_source_area_light_sign"]=(
            "positive" if terms["log_light:log10_area_c"]["beta"]>0 else "negative_or_zero")
        stats["union"]["independent_new_source_area_moderator_95pct_ci_excludes_zero"]=(
            terms["log_light:log10_area_c"]["ci_low"]>0 or
            terms["log_light:log10_area_c"]["ci_high"]<0)
    result={
      "schema":SCHEMA,
      "photo_source_years":list(YEARS),
      "independent_source_observation_UUIDs":36449,
      "source_original_numeric_9palette_ROI9_successes":13842,
      "source_photo_cell_years_strict_after_tier06_tier07_union":1137,
      "source_last_verified_73_SHA_signed_original_image_results":proof,
      "source_sameyear_environment_cells":len(environment),
      "source_model_results":stats,
      "original_historical_strict_photo_source_geographic_overlap":overlap,
      "shared_same_island_site_year_old_new_source_taxon_and_white_ROI_repeatability":paired_photography,
      "source_photo_white_ROI_not_verified_native_plant_population_white_fraction":True,
      "island_and_site_selected_acquisition_not_unbiased_worldwide_flora_sampling":True,
      "independent_original_photo_IDs_not_independent_island_environment_realizations":True,
      "no_original_photo_image_bytes_downloaded_or_reinferred":True,
      "no_missing_photo_colour_imputation":True,
      "historical_main_fiveyear_ALAN_coefficients_unchanged":True,
      "source_SHA256":{**{k:panel.sha(v) for k,v in paths.items()},
          "original_coast_GSHHG_sha256":isolation.COAST_SHA,
          "source_tier06_photo_schedule":union.SHA06,
          "source_tier07_photo_schedule":union.SHA07}
    }
    out.mkdir(parents=True)
    file=out/"ALAN_36449_new_photo_source_island_nightlight_replication.json"
    file.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    display={k:{"n_strict":v["actual_5attempt_3ROI9_3accepted_photo_cell_years"],
                "n_islands":v["strict_photo_islands"],
                "model":v["sameyear_exposure_join_and_area_isolation_model"]}
             for k,v in stats.items()}
    print("ALAN_NEW36449_ORIGINAL_PHOTO_ISLAND_AREA_ALAN_REPLICATION",
          json.dumps(display,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp","geo2019","viirs2019","pop2019",
          "geo2021","viirs2021","pop2021","geo2023","viirs2023","pop2023",
          "geo2024","cn2024","uc2024","uc2024receipt",
          "geo2025","geo2025receipt","viirs2025")
    for key in keys:
        parser.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    for item in ("coast","source06","source07","prior06","bulk06",
                 "pilot07","bulk07","out"):
        parser.add_argument("--"+item,type=Path,required=True)
    args=parser.parse_args()
    run({k:getattr(args,k) for k in keys},args.coast,args.source06,args.source07,
        args.prior06,args.bulk06,args.pilot07,args.bulk07,args.out)
