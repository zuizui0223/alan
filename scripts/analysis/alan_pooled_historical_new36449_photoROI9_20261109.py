#!/usr/bin/env python3
"""ALAN observational ORIGINAL-photo pooled historical + source-disjoint tier06/07.

The historical 1,032,882-observation canonical registry is NOT entirely the
historical-main ecological model sample. First reconstruct only original
'historical_main' photographs at SHA-authenticated, unambiguous GSHHG
0.05-degree islands in the 5 comparable years. Then combine actual DISTINCT
photo+observation UUIDs with tier06 8,942 and tier07 27,507 measured images.
Recompute accepted-taxon-equal white ROI at each EXACT island/site/year, including
previously individually ineligible sites that acquire genuine source support.
Do not count overlapping site-years twice or treat old/new photo arms as
independent island/environment realizations.

This is a deliberately descriptive, *non-probability-weighted* pooled-photo
model; keep the old-only and new-only estimands and an equal-SOURCE-arm
within-cell sensitivity separately. No imputation, night-pollination claims or
redistribution of copyrighted image pixels.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import alan_new36449_photoROI9_island_light_replication_20261109 as fresh
from scripts.analysis import alan_union_36449_original_island_photoROI9_20261109 as union
from scripts.analysis import alan_isolation_area_light_20261108 as isolation
from scripts.analysis import alan_island_area_light_20261108 as area
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as historic
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis.mois_2024_sameyear_global_island_photo_white_20261107 import cell_id

SCHEMA = "ALAN_actual_historical_new36449_original_photo_pooled_island_ROI9_v1"
YEARS = (2019, 2021, 2023, 2024, 2025)
PHOTO_N_NEW = 36449
NEW_ROI_N = 13842
STRICT_OLD = 710
STRICT_NEW = 1137
SITES_SHARED = 301
ISLANDS_SHARED = 44
AREA_TERM = "log_light:log10_area_c"
SOURCE_OLD = "historical_main"
SOURCE_NEW = "new_tier06_tier07"

def require(ok, message):
    if not bool(ok):
        raise ValueError("ALAN_OLD_NEW_REAL_PHOTO_POOLED: " + message)

def unique_original_gev(geo_path, year):
    """Use globally unambiguous GSHHG 0.05-degree cells, not guessed island IDs."""
    g = pd.read_csv(geo_path, dtype={"cell_id":str, "target_id":str,
                                    "land_id":str}, low_memory=False)
    require({"cell_id","target_id","land_id","land_class"} <= set(g.columns),
            "Original GSHHG geography has no actual cell / island membership")
    require(g.target_id.is_unique, "GSHHG target-year source contains duplicates")
    require(g.target_id.str.startswith(str(year)+"|").all(),
            "GSHHG source year does not match original photo observation")
    ambiguous = g.cell_id.duplicated(keep=False)
    u = g.loc[~ambiguous & g.land_class.eq("island"),
              ["cell_id","target_id","land_id"]].copy()
    require(u.cell_id.is_unique and u.target_id.is_unique,
            "A mixed-island grid cell was assigned to an arbitrary island")
    return u, {"year":year,"ambiguous_original_geographic_rows":int(ambiguous.sum()),
               "safe_original_island_cells":len(u)}

def old_original_photos(canonical, wcvp, paths):
    """Reconstruct exact old-photo source, including quality failures as attempts."""
    require(len(canonical) == 1032882 and canonical.photo_id.is_unique and
            canonical.observation_id.is_unique,
            "Frozen historical canonical 1,032,882 original photo IDs drift")
    base = canonical.loc[
       canonical.source_subarm.eq(SOURCE_OLD) &
       canonical.latitude.between(-90,90) &
       canonical.longitude.between(-180,180)
    ].copy()
    base["year"] = pd.to_numeric(
        base.observed_on.astype(str).str[:4], errors="coerce")
    base = base.loc[base.year.isin(YEARS)].copy()
    base["year"] = base.year.astype(int)
    base["cell_id"] = cell_id(base.longitude,base.latitude)
    lookup = wcvp.loc[
       wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
       ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    require(lookup.inat_taxon_id.is_unique, "WCVP accepted species join duplicated")
    frames, audits = [], []
    for year in YEARS:
        geo, audit = unique_original_gev(paths[f"geo{year}"],year)
        q=base.loc[base.year.eq(year)].merge(
            geo,on="cell_id",how="inner",validate="many_to_one")
        q=q.merge(lookup,on="inat_taxon_id",how="left",validate="many_to_one")
        require(q.photo_id.is_unique and q.observation_id.is_unique,
                "A historical photo was duplicated by island or accepted-taxon join")
        require(q.strict9.notna().all(), "Original historical strict9 label missing")
        q["real_numeric_ROI9"]=q.strict9.astype(bool)
        q["source_photo_white_fraction"]=pd.to_numeric(
            q.strict9_white,errors="coerce").where(q.real_numeric_ROI9)
        require(q.loc[q.real_numeric_ROI9 &
                      q.wcvp_accepted_plant_name_id.notna(),
                      "source_photo_white_fraction"].between(0,1).all(),
                "Historic original strict ROI9 has missing/invalid white mass")
        q["photo_source"]=SOURCE_OLD
        fields=["photo_id","observation_id","target_id","land_id","year",
                "wcvp_accepted_plant_name_id","real_numeric_ROI9",
                "source_photo_white_fraction","photo_source"]
        frames.append(q[fields].copy())
        audit["old_original_photos_with_unambiguous_GSHHG_island"]=len(q)
        audit["real_numeric_photo_ROI9_original"]=int(q.real_numeric_ROI9.sum())
        audits.append(audit)
    q=pd.concat(frames,ignore_index=True)
    require(q.photo_id.is_unique and q.observation_id.is_unique,
            "Old independent observation IDs repeated across years")
    return q, audits

def new_original_photos(source06,source07,prior06,bulk06,pilot07,bulk07):
    a=union.source_rows(source06,union.SHA06,"tier06",8942)
    b=union.source_rows(source07,union.SHA07,"tier07",27507)
    source=pd.concat([a,b],ignore_index=True)
    measured, signed = fresh.true_source_colour(
        prior06,bulk06,pilot07,bulk07)
    joined=source.merge(measured,on=["photo_id","observation_id"],
                        validate="one_to_one")
    require(len(joined)==PHOTO_N_NEW and joined.photo_id.is_unique and
            joined.observation_id.is_unique and
            int(joined.real_numeric_ROI9.sum())==NEW_ROI_N,
            "New independent actual ROI9 source photos not complete")
    joined["year"]=joined.year.astype(int)
    joined["photo_source"]=SOURCE_NEW
    keep=["photo_id","observation_id","target_id","land_id","year",
          "wcvp_accepted_plant_name_id","real_numeric_ROI9",
          "source_photo_white_fraction","photo_source"]
    return joined[keep].copy(),signed

def photo_cell_stats(photos):
    """Equal accepted taxa within each geographical/observation-year cell."""
    require(photos.photo_id.is_unique and photos.observation_id.is_unique,
            "Repeated independent photo or observation UUID in pooled sample")
    grouped=photos.groupby("target_id",observed=True).agg(
        photo_attempts=("observation_id","nunique"),
        strict9_success=("real_numeric_ROI9","sum"),
        year=("year","first"),land_id=("land_id","first")
    ).reset_index()
    numeric=photos.loc[
        photos.real_numeric_ROI9 &
        photos.wcvp_accepted_plant_name_id.notna() &
        photos.source_photo_white_fraction.notna()].copy()
    require(numeric.source_photo_white_fraction.between(0,1).all(),
            "Actual flower ROI9 fractions must be between zero and one")
    species=numeric.groupby(
        ["target_id","wcvp_accepted_plant_name_id"],observed=True
    ).source_photo_white_fraction.mean().reset_index()
    cells=species.groupby("target_id",observed=True).agg(
        white_equal=("source_photo_white_fraction","mean"),
        accepted_species_n=("wcvp_accepted_plant_name_id","nunique")
    ).reset_index()
    grouped=grouped.merge(cells,on="target_id",how="left",validate="one_to_one")
    grouped["accepted_species_n"]=grouped.accepted_species_n.fillna(0).astype(int)
    grouped["strict_photo_gate"]=(grouped.photo_attempts.ge(5)&
            grouped.strict9_success.ge(3)&
            grouped.accepted_species_n.ge(3))
    good=grouped.loc[grouped.strict_photo_gate].copy()
    good["white_pp"]=100*good.white_equal
    require(good.white_pp.between(0,100).all(), "Unmeasured photographic white value imputed")
    return grouped,good

def environment_qualified(cells, environment):
    e=cells.merge(environment,on=["target_id","land_id","year"],
                  how="left",validate="one_to_one")
    missing=int(e.environment_complete.isna().sum())
    invalid=int((~e.environment_complete.fillna(False)).sum())
    e=e.loc[e.environment_complete.eq(True) &
            e.land_class.eq("island")].copy()
    e["is_island"]=1
    e["log_attempts"]=np.log1p(e.photo_attempts.astype(float))
    require(e[["white_pp","log_light","log_pop","log_rain",
                "bio5","elevation","log_attempts"]].notna().all().all(),
            "Genuine source-matched photographic and environmental values missing")
    return e,{"original_strict_without_environment_target":missing,
              "original_strict_with_incomplete_environment":invalid}

def fit_photo_arm(cells,land,distances,name):
    if cells.land_id.nunique()<10 or len(cells)<80:
        return {"estimable":False,"name":name,
                "all_strict_cell_years":len(cells),
                "all_strict_islands":int(cells.land_id.nunique())}
    allisland,within,g=area.join_area_and_select(cells,land)
    if within.land_id.nunique()<10 or len(within)<80:
        return {"estimable":False,"name":name,
                "all_strict_cell_years":len(cells),
                "all_strict_islands":int(cells.land_id.nunique()),
                "n_within_island_photo_cells":len(within),
                "n_multisite_islands":int(within.land_id.nunique())}
    df=isolation.complete_island_frame(within,distances)
    centre_a=float(np.log10(df.groupby("land_id").area_km2.first()).median())
    centre_d=float(np.log10(1+df.groupby("land_id")
                         .distance_mainland_vertex_km.first()).median())
    models={}
    for weight,eq in (("cell_equal",False),("island_equal",True)):
        models["joint_"+weight]=isolation.fit_light_moderation(
            df,centre_a,centre_d,"joint",eq)
    top=str(df.land_id.value_counts().index[0])
    omit=df.loc[df.land_id.ne(top)]
    if omit.land_id.nunique()>=10 and len(omit)>=80:
        models["leave_largest_photo_island_out_equal_island"]=(
            isolation.fit_light_moderation(
            omit,centre_a,centre_d,"joint",True))
    return {"estimable":True,"name":name,
            "n_strict_all_island_cells":len(allisland),
            "n_strict_all_islands":int(allisland.land_id.nunique()),
            "n_spatial_multisite_islands":int(within.land_id.nunique()),
            "n_spatial_multisite_cells":len(within),
            "largest_island_share_of_multisite_cells":float(
                  df.land_id.value_counts().iloc[0]/len(df)),
            "original_source_years":sorted(map(int,df.year.unique())),
            "models":models}

def per_source_balanced_white(old_photos,new_photos,pooled):
    """Mean historical/new source-specific taxon-equal white values PER SITE.

    This target is NOT equivalent to pooling photo records; it is a
    sensitivity against extra new photos dominating on shared sites.
    """
    components=[]
    for name,source in ((SOURCE_OLD,old_photos),(SOURCE_NEW,new_photos)):
        _,g=photo_cell_stats(source)
        # A cell with fewer than 5 photo attempts could still contain
        # genuine accepted species/ROI9 and is valid for descriptive
        # within-cell source-mean, but cannot enter the source-only strict
        # outcome. Reconstruct the no-strict-gate species-equal value.
        allcells,_=photo_cell_stats(source)
        components.append(allcells[["target_id","white_equal"]].rename(
            columns={"white_equal":name}))
    both=pooled.merge(components[0],on="target_id",how="left",
                      validate="one_to_one").merge(
         components[1],on="target_id",how="left",validate="one_to_one")
    present=both[[SOURCE_OLD,SOURCE_NEW]].notna().sum(axis=1)
    require(present.ge(1).all(),"Pooled strict cell lacks authentic source-specific colour")
    both["white_pp"]=100*both[[SOURCE_OLD,SOURCE_NEW]].mean(axis=1)
    return both,{"pooled_strict_cells_with_both_source_numeric_outcomes":int(
        present.eq(2).sum()),
        "pooled_strict_cells_with_single_source_numeric_outcome":int(
        present.eq(1).sum())}

def run(paths,coast,source06,source07,prior06,bulk06,pilot07,bulk07,out):
    out=Path(out)
    require(not out.exists(), "Refuse overwriting signed true original photo pooled model")
    canonical,wcvp=panel.load_frozen_photos(paths)
    environment=fresh.original_sameyear_environment(paths)
    old,byyear=old_original_photos(canonical,wcvp,paths)
    new,signed=new_original_photos(
       source06,source07,prior06,bulk06,pilot07,bulk07)
    require(not set(old.photo_id).intersection(set(new.photo_id)) and
            not set(old.observation_id).intersection(set(new.observation_id)),
            "historical_main and both new photo source arms share photo/observation IDs")
    # Independently validate exact old and new strict sampled site/year cells.
    old_allsites,old_sites=photo_cell_stats(old)
    new_allsites,new_sites=photo_cell_stats(new)
    prior=historic.by_year_tier(historic.build_frozen_five_years(paths))[
        "strict_5p3r3s"]
    prior=prior.loc[prior.is_island.eq(1)].copy()
    old_env,old_invalid=environment_qualified(old_sites,environment)
    new_env,new_invalid=environment_qualified(new_sites,environment)
    require(len(prior)==STRICT_OLD and
            set(prior.target_id)==set(old_env.target_id),
            "Historical-main 710 true strict ROI photo cell-years NOT reproduced")
    white_diff=old_env[["target_id","white_pp"]].merge(
        prior[["target_id","white_pp"]],on="target_id",
        suffixes=("_reconstructed","_historical"),validate="one_to_one")
    require(np.allclose(white_diff.white_pp_reconstructed,
                        white_diff.white_pp_historical,atol=1e-6,rtol=0),
            "Historical original accepted-species equal white mass has drifted")
    require(len(new_env)==STRICT_NEW and
            old_env.land_id.nunique()==101 and
            new_env.land_id.nunique()==75,
            "Actual new strict 1,137 photo cell-years/75 islands not reproduced")
    shared_sites=set(old_env.target_id)&set(new_env.target_id)
    shared_islands=set(old_env.land_id)&set(new_env.land_id)
    require(len(shared_sites)==SITES_SHARED and
            len(shared_islands)==ISLANDS_SHARED,
            "Exact historical/new overlap must be measured without double-counting")
    pooled_photos=pd.concat([old,new],ignore_index=True)
    all_pooled,eligible_pooled=photo_cell_stats(pooled_photos)
    pooled_env,pooled_invalid=environment_qualified(eligible_pooled,environment)
    expected_union=set(old_env.target_id)|set(new_env.target_id)
    actual_union=set(pooled_env.target_id)
    require(expected_union.issubset(actual_union) and
            pooled_env.target_id.is_unique,
            "Combined original photographs lost previously strict old/new sites")
    added=actual_union-expected_union
    require(len(expected_union)==STRICT_OLD+STRICT_NEW-SITES_SHARED,
            "Combined old/new strict sites must de-duplicate exact island cell/year")
    require(pooled_env.land_id.nunique()>=len(shared_islands),
            "Pooled original photo site island count fell below common islands")
    balanced,balance_audit=per_source_balanced_white(old,new,pooled_env)
    land=isolation.load_signed_land(coast)
    landmeta=land[["land_id","land_class","area_km2"]].copy()
    all_island_ids=set(pooled_env.land_id)
    distances=isolation.isolated_mainland_distance_km(land,all_island_ids)
    models={}
    for name,d in (("historical_main_only",old_env),
                   ("new_36449_only",new_env),
                   ("pooled_all_true_photo_records",pooled_env),
                   ("pooled_equal_source_within_cell",balanced)):
        models[name]=fit_photo_arm(d,landmeta,distances,name)
    old_beta=models["historical_main_only"]["models"]["joint_island_equal"][
        "terms"][AREA_TERM]["beta"]
    new_beta=models["new_36449_only"]["models"]["joint_island_equal"][
        "terms"][AREA_TERM]["beta"]
    require(abs(old_beta-5.383529291)<0.0001 and
            abs(new_beta-(-1.090476086034753))<0.0001,
            "Historical or new independent model no longer numerically reproducible")
    # Common island support (not all shared islands qualify for within-island
    # site contrasts separately): never force a statistically identified fit.
    common_islands=set(old_env.land_id)&set(new_env.land_id)
    matched={}
    for label,frame in (("historical_on_common_islands",old_env),
                        ("new_on_common_islands",new_env),
                        ("pooled_on_common_islands",pooled_env)):
        eligible=frame.loc[frame.land_id.isin(common_islands)]
        matched[label]=fit_photo_arm(eligible,landmeta,distances,label)
    counts={
      "historical_canonical_all_source_subarms_observation_UUIDs":len(canonical),
      "historical_main_original_island_photo_attempts_in_unambiguous_GSHHG_cells":len(old),
      "historical_main_original_numeric_strict9_island_photos":int(old.real_numeric_ROI9.sum()),
      "new_tier06_tier07_separate_observation_and_photo_UUIDs":len(new),
      "new_tier06_tier07_original_numeric_ROI9_photos":int(new.real_numeric_ROI9.sum()),
      "actual_pooled_old_and_new_independent_original_island_photo_UUIDs":len(pooled_photos),
      "actual_pooled_old_and_new_original_numeric_ROI9_photos":int(
          pooled_photos.real_numeric_ROI9.sum()),
      "historical_main_strict_island_photo_cells":len(old_env),
      "new_36449_strict_island_photo_cells":len(new_env),
      "strict_island_photo_sites_in_both_source_arms":len(shared_sites),
      "strict_GSHHG_islands_in_both_source_arms":len(shared_islands),
      "strict_unique_sites_in_source_union_before_photo_recombination":len(expected_union),
      "actual_strict_unique_source_pooled_photo_site_years":len(pooled_env),
      "extra_strict_sites_only_after_genuine_old_new_photo_pooling":len(added),
      "strict_unique_GSHHG_islands_from_photo_pool":int(pooled_env.land_id.nunique()),
      "all_photo_attempt_source_island_site_years":len(all_pooled),
      "old_only_source_photo_attempt_island_site_years":len(old_allsites),
      "new_only_source_photo_attempt_island_site_years":len(new_allsites),
      "old_new_photo_UUID_sets_disjoint":True,
      "same_photo_site_and_year_never_double_counted_in_regression":True,
    }
    response={
      "schema":SCHEMA,
      "actual_original_source_photo_years":list(YEARS),
      "counts":counts,
      "old_historical_main_original_photo_geography_audit_by_year":byyear,
      "old_original_environment_quality":old_invalid,
      "new_original_environment_quality":new_invalid,
      "pooled_original_environment_quality":pooled_invalid,
      "per_site_equal_source_white_sensitivity_support":balance_audit,
      "models":models,
      "common_island_support_models":matched,
      "old_island_area_moderator_reproduced_from_original_photo_records":True,
      "new_island_area_moderator_reproduced_from_original_photo_records":True,
      "all_original_photos_inferred_again":False,
      "no_technical_colour_unknown_imputed":True,
      "original_photo_unique_UUID_no_double_count":True,
      "pooled_photo_estimate_conditions_on_island_site_targeted_sample":True,
      "source_design_not_random_global_native_flora":True,
      "source_photo_white_is_not_botanical_white_flower_prevalence":True,
      "no_ALAN_causality_or_night_pollinator_effect_claim":True,
      "source_hashes":{"source06":union.SHA06,
                       "source07":union.SHA07,
                       "canonical":panel.sha(paths["canonical"]),
                       "wcvp":panel.sha(paths["wcvp"])},
      "model_formula":isolation.BASE_COVARS,
    }
    out.mkdir(parents=True)
    output=out/"ALAN_old_plus_new36449_real_photoROI9_pooled_area_light.json"
    output.write_text(json.dumps(response,sort_keys=True,indent=2)+"\n")
    short={
       "counts":counts,
       "area_and_isolation_joint_by_source":{
          name:({k:model["models"][k]["terms"][AREA_TERM]
                 for k in ("joint_cell_equal","joint_island_equal")}
                if model["estimable"] else {"estimable":False})
          for name,model in models.items()},
       "matched_island_source_coverage":{
          k:{"estimable":v["estimable"],
             "n_islands":v.get("n_spatial_multisite_islands"),
             "n_cells":v.get("n_spatial_multisite_cells")}
          for k,v in matched.items()},
       "equal_source_per_cell_numeric_overlap":balance_audit,
    }
    print("ALAN_TRUE_HISTORICAL_AND_NEW36449_ORIGINAL_PHOTO_POOLED_MODEL",
          json.dumps(short,sort_keys=True),flush=True)
    return response

def main():
    p=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp",
          "geo2019","viirs2019","pop2019",
          "geo2021","viirs2021","pop2021",
          "geo2023","viirs2023","pop2023",
          "geo2024","cn2024","uc2024","uc2024receipt",
          "geo2025","geo2025receipt","viirs2025")
    for key in keys:
        p.add_argument("--"+key.replace("_","-"),dest=key,
                       type=Path,required=True)
    for item in ("coast","source06","source07","prior06",
                 "bulk06","pilot07","bulk07","out"):
        p.add_argument("--"+item,type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in keys},
        a.coast,a.source06,a.source07,a.prior06,a.bulk06,
        a.pilot07,a.bulk07,a.out)

if __name__=="__main__":
    main()
