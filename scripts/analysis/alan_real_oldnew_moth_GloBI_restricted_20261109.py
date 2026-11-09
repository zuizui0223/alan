#!/usr/bin/env python3
"""Actually fit historical and newly measured original photos restricted by GloBI.

Old full-flora island-area x ALAN was NOT restricted to moth-associated plants.
Read signed old + tier06/07 + tier08 910 true photo numeric ROI9, de-duplicate
exact source photo IDs, and match GloBI literature with accepted WCVP names.
For each positive evidence tier, select photos BEFORE numerical ROI site gates.
Moth family visitor != true local effective nocturnal pollination. No record
!= verified moth absence. Technical-unknown photo ROI is never imputed.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd

from scripts.analysis import alan_pooled_historical_new36449_photoROI9_20261109 as oldnew
from scripts.analysis import alan_pooled_plus_tier08_910_original_photo_20261109 as tier08
from scripts.analysis import alan_tier08_moth_positive_source_support_20261109 as moth
from scripts.analysis import alan_new36449_photoROI9_island_light_replication_20261109 as recent
from scripts.analysis import alan_isolation_area_light_20261108 as iso
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

SCHEMA="ALAN_verified_87852_original_photo_GloBI_moth_evidence_area_light_restrict_v1"
TERM="log_light:log10_area_c"
STRATA=("all_original_photographed_flora",
        "moth_record_positive",
        "moth_database_pollinates",
        "walton_direct_nocturnal_pollen_carriage",
        "moth_record_without_recorded_butterfly",
        "butterfly_record_positive")

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_MOTH_RESTRICTED_ACTUAL_PHOTOS: "+msg)

def accepted_guilds(photos,accepted,globi):
    names=accepted[["wcvp_accepted_plant_name_id","wcvp_accepted_name"]].drop_duplicates()
    require(names.wcvp_accepted_plant_name_id.is_unique,
            "WCVP accepted IDs do not uniquely map to plant names")
    q=photos.merge(names,on="wcvp_accepted_plant_name_id",
                   how="left",validate="many_to_one")
    flags=moth.evidence_flags(globi)
    q=q.merge(flags,left_on="wcvp_accepted_name",right_on="plant_name",
              how="left",validate="many_to_one")
    require(len(q)==len(photos) and q.photo_id.is_unique,
            "GloBI plant-guild join altered true original image IDs")
    for k in moth.FLAGS:
        q[k]=q[k].fillna(False).astype(bool)
    q["all_original_photographed_flora"]=True
    q["moth_record_without_recorded_butterfly"]=(
        q.moth_record_positive & ~q.butterfly_record_positive)
    return q

def site_and_model(part,environment,land,distances,name):
    if part.empty:
        return {"estimable":False,"source_original_photo_attempts":0}
    allsites,strict=oldnew.photo_cell_stats(part)
    site,env_audit=oldnew.environment_qualified(strict,environment)
    fitted=oldnew.fit_photo_arm(site,land,distances,name)
    count={
       "source_original_photo_attempts":len(part),
       "source_true_numeric_photo_ROI9":int(part.real_numeric_ROI9.sum()),
       "source_GloBI_matched_accepted_species":int(
          part.wcvp_accepted_plant_name_id.dropna().nunique()),
       "all_original_photo_island_cell_years_before_ROI_gates":len(allsites),
       "strict_original_ROI9_island_cell_years":len(site),
       "strict_GSHHG_islands":int(site.land_id.nunique()),
       "source_environment_quality_audit":env_audit}
    if not fitted.get("estimable"):
        return {"estimable":False,"source":count,"reason":fitted.get("reason"),
                "n_multisite_islands":fitted.get("n_multisite_islands",0),
                "n_spatial_cells":fitted.get("n_within_island_photo_cells",0)}
    a=fitted["models"]["joint_island_equal"]["terms"][TERM]
    b=fitted["models"]["joint_cell_equal"]["terms"][TERM]
    return {"estimable":True,"source":count,
        "n_multisite_islands":fitted["n_spatial_multisite_islands"],
        "n_spatial_cells":fitted["n_spatial_multisite_cells"],
        "island_equal_area_x_ALAN":a,"cell_equal_area_x_ALAN":b,
        "leave_largest_photo_island_out_area_x_ALAN":(
            fitted["models"].get("leave_largest_photo_island_out_equal_island",{})
            .get("terms",{}).get(TERM))}

def run(paths,coast,source06,source07,prior06,bulk06,pilot07,bulk07,
        tier08file,tier08receipt,acceptedfile,globifile,out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite real moth-guild photo analysis")
    require(panel.sha(acceptedfile)==moth.SHA["accepted"] and
            panel.sha(globifile)==moth.SHA["globi"],
            "Source accepted botanical names or GloBI 61,967 interactions changed")
    accepted=pd.read_csv(acceptedfile,compression="gzip",dtype=str,
        usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"])
    globi=pd.read_csv(globifile,compression="gzip",dtype=str)
    require(len(globi)==61967,"Original documented source evidence changed")
    canonical,wcvp=panel.load_frozen_photos(paths)
    old,_=oldnew.old_original_photos(canonical,wcvp,paths)
    new,_=oldnew.new_original_photos(
        source06,source07,prior06,bulk06,pilot07,bulk07)
    overlap=set(old.photo_id)&set(new.photo_id)
    require(len(old)==50829 and len(new)==36449 and len(overlap)==329,
            "Original historical/new photograph overlap not reproduced")
    clean_new=new.loc[~new.photo_id.isin(overlap)]
    first=pd.concat([old,clean_new],ignore_index=True)
    raw08,_=tier08.load_real_tier08(tier08file,tier08receipt)
    clean08,check=tier08.same_photo_reconciliation(first,raw08)
    full=pd.concat([first,clean08],ignore_index=True)
    require(len(clean08)==903 and len(full)==87852 and
            full.photo_id.is_unique and int(full.real_numeric_ROI9.sum())==32111,
            "Historical plus tier08 signed source counts changed")
    x=accepted_guilds(full,accepted,globi)
    env=recent.original_sameyear_environment(paths)
    allsites,strict=oldnew.photo_cell_stats(full)
    eligible,quality=oldnew.environment_qualified(strict,env)
    require(len(eligible)==2530 and eligible.land_id.nunique()==193,
            "All-flora strict cohort 2530 cells or 193 islands changed")
    landmass=iso.load_signed_land(coast)
    land=landmass[["land_id","land_class","area_km2"]].copy()
    distances=iso.isolated_mainland_distance_km(landmass,set(eligible.land_id))
    source={
      "old_original":x.loc[x.photo_source.eq("historical_main")],
      "new_tier06_07_photoID_clean":x.loc[
          x.photo_source.eq("new_tier06_tier07")],
      "tier08_910_photoID_clean":x.loc[
          x.photo_source.eq(tier08.TIER08_ARM)],
      "pooled_exact_photo_ID_unique":x}
    require(list(map(len,source.values()))==[50829,36120,903,87852],
            "Source group changed while adding positive moth guild evidence")
    results={}
    for stratum in STRATA:
        results[stratum]={name:site_and_model(
            group.loc[group[stratum]],env,land,distances,name+"__"+stratum)
            for name,group in source.items()}
    baseline=results["all_original_photographed_flora"][
        "pooled_exact_photo_ID_unique"]
    require(baseline["estimable"] and
            abs(baseline["island_equal_area_x_ALAN"]["beta"]-
                1.1059197451811407)<0.00001 and
            baseline["n_multisite_islands"]==97 and
            baseline["n_spatial_cells"]==2419,
            "Old all-flora actual photo baseline not numerically identical")
    summary={}
    for group,by_source in results.items():
        t=by_source["pooled_exact_photo_ID_unique"]
        summary[group]={
          "original_photo_attempts":t["source"]["source_original_photo_attempts"],
          "true_numeric_original_ROI9_photos":t["source"]["source_true_numeric_photo_ROI9"],
          "accepted_species":t["source"]["source_GloBI_matched_accepted_species"],
          "strict_site_years":t["source"]["strict_original_ROI9_island_cell_years"],
          "strict_islands":t["source"]["strict_GSHHG_islands"],
          "regression_estimable":t["estimable"],
          "multisite_islands":t.get("n_multisite_islands",0),
          "multisite_cell_years":t.get("n_spatial_cells",0),
          "island_equal_interaction":t.get("island_equal_area_x_ALAN")}
    report={
      "schema":SCHEMA,"years":[2019,2021,2023,2024,2025],
      "historical_old_area_ALAN_main_model_had_GloBI_moth_filter":False,
      "real_unique_original_photo_IDs":87852,
      "true_numeric_FCP_original_photo_ROI9":32111,
      "original_unstratified_base_source":quality,
      "summary":summary,
      "source_model_detail":results,
      "moth_visitor_not_effective_exclusive_nocturnal_pollination":True,
      "GloBI_no_record_not_true_non_moth_control":True,
      "butterfly_and_moth_positive_same_species_may_overlap":True,
      "historical_numeric_ID_new_UUID_obs_independence_unverified":True,
      "unknown_flower_photo_colour_not_imputed":True,
      "source_nonprobabilistic_island_site_image_selection":True,
      "no_causal_night_light_flower_evolution_claim":True,
      "no_image_retrieval_or_photo_RGB_reinference":True,
      "source_sha256":{"globi":moth.SHA["globi"],"accepted":moth.SHA["accepted"]}
    }
    out.mkdir(parents=True)
    (out/"ALAN_87852_original_photo_moth_visitor_restricted_area_light.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_TRUE_ORIGINAL_ALL_FLORA_VS_MOTH_GLOBI_RESTRICTED_OLD_NEW",
          json.dumps(summary,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    keys=("canonical","wcvp",
          "geo2019","viirs2019","pop2019","geo2021","viirs2021","pop2021",
          "geo2023","viirs2023","pop2023","geo2024","cn2024","uc2024","uc2024receipt",
          "geo2025","geo2025receipt","viirs2025")
    for name in (*keys,"coast","source06","source07","prior06","bulk06",
                 "pilot07","bulk07","tier08","tier08-receipt",
                 "accepted","globi","out"):
        p.add_argument("--"+name,type=Path,required=True,
                       dest=name.replace("-","_"))
    a=p.parse_args()
    run({k:getattr(a,k) for k in keys},a.coast,a.source06,a.source07,
        a.prior06,a.bulk06,a.pilot07,a.bulk07,a.tier08,a.tier08_receipt,
        a.accepted,a.globi,a.out)
if __name__=="__main__":main()
