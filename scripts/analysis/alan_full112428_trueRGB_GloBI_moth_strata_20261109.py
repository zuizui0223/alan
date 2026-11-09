#!/usr/bin/env python3
"""Recalculate GloBI moth-evidence photo ROI strata using all 112428 original photos.

Source-disjoint tier08 24,576 genuine original RGB/ROI9 attempts are SHA-audited
and appended to the original 87,852 unique photo IDs. No technical unknown is
converted into an observed flower. GloBI positive visits != local nighttime
effective pollination; absent relations are NOT absence controls.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from scripts.analysis import alan_pooled_historical_new36449_photoROI9_20261109 as historic
from scripts.analysis import alan_pooled_plus_tier08_910_original_photo_20261109 as prior08
from scripts.analysis import alan_tier08_24576_original_photo_stability_20261109 as wavescript
from scripts.analysis import alan_real_oldnew_moth_GloBI_restricted_20261109 as moth
from scripts.analysis import alan_new36449_photoROI9_island_light_replication_20261109 as fresh
from scripts.analysis import alan_isolation_area_light_20261108 as isolation
from scripts.analysis import alan_island_area_light_20261108 as area
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

SCHEMA="ALAN_true112428_original_photo_RGB_GloBI_moth_guild_stratified_v1"
TERM="log_light:log10_area_c"


def require(v, message):
    if not bool(v):
        raise ValueError("ALAN_MOTH_112428_TRUE_RGB: "+message)


def fit_subset(part, environment, land, distances, name, original97=None):
    """Apply the ROI gates after restricting genuinely measured source photos."""
    if part.empty:
        return {"estimable": False, "original_photo_attempts": 0,
                "true_numeric_ROI9": 0, "strict_islands": 0,
                "strict_site_years": 0, "multisite_islands": 0,
                "multisite_site_years": 0}
    _, strict = historic.photo_cell_stats(part)
    cells, quality = historic.environment_qualified(strict, environment)
    fit = historic.fit_photo_arm(cells, land, distances, name)
    result = {
        "original_photo_attempts": int(len(part)),
        "true_numeric_ROI9": int(part.real_numeric_ROI9.sum()),
        "matched_accepted_plant_taxa": int(part.wcvp_accepted_plant_name_id.dropna().nunique()),
        "strict_islands": int(cells.land_id.nunique()),
        "strict_site_years": int(len(cells)),
        "multisite_islands": int(fit.get("n_spatial_multisite_islands", fit.get("n_multisite_islands", 0))),
        "multisite_site_years": int(fit.get("n_spatial_multisite_cells", fit.get("n_within_island_photo_cells", 0))),
        "estimable": bool(fit.get("estimable")),
        "source_quality": quality,
        "nightlight_is_same_year_of_original_photo": True,
    }
    if fit.get("estimable"):
        result["island_equal_area_ALAN"] = fit["models"]["joint_island_equal"]["terms"][TERM]
        result["cell_equal_area_ALAN"] = fit["models"]["joint_cell_equal"]["terms"][TERM]
        result["leave_largest_island_out"] = (
            fit["models"].get("leave_largest_photo_island_out_equal_island",{})
            .get("terms",{}).get(TERM))
        if original97:
            common=cells.loc[cells.land_id.isin(original97)].copy()
            cf=historic.fit_photo_arm(common,land,distances,name+"_fixed_original_97_islands")
            result["fixed_original_97_island_geographic_sensitivity"] = (
                cf["models"]["joint_island_equal"]["terms"][TERM]
                if cf.get("estimable") else {"estimable":False})
    else:
        result["not_estimable_reason"]=fit.get("reason")
    return result


def run(paths, coast, source06,source07,prior06,bulk06,pilot07,bulk07,
        tier08,tier08_receipt,large_receipt,waves,
        acceptedfile,globifile,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite signed source-only result")
    require(panel.sha(acceptedfile)==moth.moth.SHA["accepted"] and
            panel.sha(globifile)==moth.moth.SHA["globi"],
            "GloBI positive evidence or WCVP accepted-name source SHA drift")
    accepted=pd.read_csv(acceptedfile,compression="gzip",dtype=str,
                         usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    globi=pd.read_csv(globifile,compression="gzip",dtype=str)
    require(len(globi)==61967, "GloBI provider records changed")
    canonical,wcvp=panel.load_frozen_photos(paths)
    old,_=historic.old_original_photos(canonical,wcvp,paths)
    recent,_=historic.new_original_photos(
        source06,source07,prior06,bulk06,pilot07,bulk07)
    shared=set(old.photo_id)&set(recent.photo_id)
    require(len(old)==50829 and len(recent)==36449 and len(shared)==329,
            "Old source/photo overlap changed")
    clean_recent=recent.loc[~recent.photo_id.isin(shared)].copy()
    base=pd.concat([old,clean_recent],ignore_index=True)
    first,_=prior08.load_real_tier08(tier08,tier08_receipt)
    clean_first,check=prior08.same_photo_reconciliation(base,first)
    require(len(clean_first)==903 and
            check["tier08_original_photo_ids_already_in_prior"]==7,
            "Earlier 910 source duplicates changed")
    base=pd.concat([base,clean_first],ignore_index=True)
    require(len(base)==87852 and base.photo_id.is_unique and
            int(base.real_numeric_ROI9.sum())==32111,
            "Frozen 87852 prior actual numeric flower ROI source changed")
    additions=[]
    for i,path in enumerate(waves,start=1):
        x,_=wavescript.verify_wave(path,large_receipt,i)
        require(not x.photo_id.isin(base.photo_id).any() and
                not x.photo_id.isin(pd.concat(additions).photo_id).any() if additions
                else not x.photo_id.isin(base.photo_id).any(),
                "Repeated original photo ID in newly measured wave")
        additions.append(x)
    allphotos=pd.concat([base,*additions],ignore_index=True)
    require(len(allphotos)==112428 and allphotos.photo_id.is_unique and
            int(allphotos.real_numeric_ROI9.sum())==40914,
            "Measured all-photo 112428/40914 source reconciliation failed")
    allphotos=moth.accepted_guilds(allphotos,accepted,globi)
    env=fresh.original_sameyear_environment(paths)
    _,strict=historic.photo_cell_stats(allphotos)
    whole,quality=historic.environment_qualified(strict,env)
    require(len(whole)==3481 and whole.land_id.nunique()==263,
            "Full verified 3481 strict sites/263 GSHHG islands did not reproduce")
    landmass=isolation.load_signed_land(coast)
    land=landmass[["land_id","land_class","area_km2"]].copy()
    distance=isolation.isolated_mainland_distance_km(landmass,set(whole.land_id))
    old97_source=allphotos.loc[allphotos.photo_source.isin((
        "historical_main","new_tier06_tier07",prior08.TIER08_ARM))]
    _,base_strict=historic.photo_cell_stats(old97_source)
    original97, _=historic.environment_qualified(base_strict,env)
    common_sites=historic.fit_photo_arm(
        original97,land,distance,"baseline_before24576_all_flora")
    require(common_sites["estimable"] and
            common_sites["n_spatial_multisite_islands"]==97,
            "Original 97 within-island ALAN geography cannot be reproduced")
    # Use the SAME within-island viable 97 islands as the frozen all-flora model;
    # merely having multiple photo sites does not prove light contrast.
    _,spatial_original=area.join_area_and_select(original97,land)
    base_islands=set(spatial_original.land_id)
    require(len(base_islands)==97, "Fixed-geography 97 island support drift")
    cohorts={
        "historical_main":allphotos.loc[allphotos.photo_source.eq("historical_main")],
        "tier06_07_deduplicated":allphotos.loc[allphotos.photo_source.eq("new_tier06_tier07")],
        "tier08_prior910":allphotos.loc[allphotos.photo_source.eq(prior08.TIER08_ARM)],
        "tier08_additional24576":allphotos.loc[allphotos.photo_source.str.startswith("tier08_new24576_wave")],
        "pooled_all_sources":allphotos,
    }
    require([len(cohorts[k]) for k in ("historical_main","tier06_07_deduplicated",
                "tier08_prior910","tier08_additional24576","pooled_all_sources")] ==
            [50829,36120,903,24576,112428], "Source arms not distinct")
    fitted={}
    for stratum in moth.STRATA:
        g={}
        for name,cohort in cohorts.items():
            subset=cohort.loc[cohort[stratum]].copy()
            g[name]=fit_subset(subset,env,land,distance,name+"_"+stratum,base_islands)
        fitted[stratum]=g
    reference=fitted["all_original_photographed_flora"]["pooled_all_sources"]
    require(reference["estimable"] and
            abs(reference["island_equal_area_ALAN"]["beta"]-0.6031352440580511)<1e-4 and
            reference["multisite_islands"]==133 and
            reference["multisite_site_years"]==3315,
            "All-flora independent 3-wave model did not reproduce")
    summary={k:v["pooled_all_sources"] for k,v in fitted.items()}
    report={
        "schema":SCHEMA,
        "original_observation_years":[2019,2021,2023,2024,2025],
        "new_original_RGB_photo_attempts":24576,
        "all_unique_original_photo_IDs":112428,
        "original_numeric_flower_ROI9":40914,
        "new_original_numeric_flower_ROI9":8803,
        "new_original_technical_unknown":15773,
        "original_112428_photo_IDs_never_duplicated":True,
        "baseline_original_geography_97_islands":True,
        "no_GloBI_record_does_not_imply_true_moth_absence":True,
        "moth_and_butterfly_visitor_evidence_overlap":True,
        "Walton_plant_record_does_not_prove_local_nocturnal_pollination":True,
        "no_causal_ALAN_flower_selection_or_native_flower_frequencies":True,
        "historical_numeric_observationID_new_UUID_crosswalk_unproven":True,
        "same_year_VIIRS_area_mainland_distance_and_island_year_FE":True,
        "source_evidence_hashes":{"accepted":moth.moth.SHA["accepted"],
                                  "globi":moth.moth.SHA["globi"]},
        "source_quality":quality,
        "summary":summary,"models_by_source_and_evidence":fitted,
    }
    out.mkdir(parents=True)
    (out/"ALAN_112428_true_RGB_GloBI_moth_strata.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_TRUE112428_GLOBI_MOTH_GUILD_ACTUAL_ORIGINAL_RGB_RECALC",
          json.dumps(summary,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp","geo2019","viirs2019","pop2019",
          "geo2021","viirs2021","pop2021","geo2023","viirs2023","pop2023",
          "geo2024","cn2024","uc2024","uc2024receipt",
          "geo2025","geo2025receipt","viirs2025")
    extras=("coast","source06","source07","prior06","bulk06","pilot07","bulk07",
            "tier08","tier08-receipt","large-receipt","wave1","wave2","wave3",
            "accepted","globi","out")
    for k in (*keys,*extras):
        p.add_argument("--"+k,type=Path,required=True,dest=k.replace("-","_"))
    a=p.parse_args()
    run({k:getattr(a,k) for k in keys},a.coast,a.source06,a.source07,
        a.prior06,a.bulk06,a.pilot07,a.bulk07,a.tier08,a.tier08_receipt,
        a.large_receipt,[a.wave1,a.wave2,a.wave3],a.accepted,a.globi,a.out)


if __name__=="__main__":
    main()
