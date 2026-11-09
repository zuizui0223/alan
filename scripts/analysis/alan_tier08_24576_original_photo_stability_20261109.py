#!/usr/bin/env python3
"""Source-authenticated stage-by-stage original-photo ALAN convergence.

Real old historical-main + signed tier06/07 + signed tier08 910 source, then
three genuinely NEW 8192-photo tier08 waves (24,576 FCP ROI9 attempts).
Reconstruct accepted-taxon-equal flower-colour ROI per exact source site/year,
same-year environmental raster, island area & mainland coast distance. All
previous photograph IDs are de-duplicated. Unknown palette is always missing.

No p-value-driven stopping: use prespecified coefficient drift/precision,
island support, common-geographic-frame and leave-dominant-island metrics.
All observed waves are kept in report even if an earlier checkpoint is stable.
No causal ALAN or native botanical white-flower prevalence inference.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis import alan_pooled_historical_new36449_photoROI9_20261109 as historic
from scripts.analysis import alan_pooled_plus_tier08_910_original_photo_20261109 as prior08
from scripts.analysis import alan_new36449_photoROI9_island_light_replication_20261109 as fresh
from scripts.analysis import alan_isolation_area_light_20261108 as isolation
from scripts.analysis import alan_island_area_light_20261108 as area
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

SCHEMA="ALAN_24576_true_RGB_tier08_3wave_island_area_light_stability_v1"
LARGE_SCHEMA="ALAN_real_original_tier08_additional24576_RGB_ROI9_3waves_v1"
ROI_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"
WAVES=(8192,16384,24576)
AREA_TERM="log_light:log10_area_c"

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_24576_CONVERGENCE: "+msg)

def verify_wave(file,receipt,wave):
    r=json.loads(Path(receipt).read_text())
    require(r.get("schema")==LARGE_SCHEMA and
            r.get("real_original_image_photo_attempts")==24576 and
            r.get("original_fixed_FCP_image_measurement_blob")==ROI_BLOB and
            r.get("real_original_measurement_shards")==192 and
            r.get("source_photo_outcomes_unknown_not_imputed") is True,
            "Original FCP photo worker 3-wave collector source schema failed")
    entry=r["waves"][wave-1]
    digest=r["wave_output_SHA256"][str(wave)]
    require(entry["wave"]==wave and
            entry["real_original_photo_attempts"]==8192 and
            entry["photo_SHA256"]==digest and panel.sha(file)==digest,
            "Measured original 8192 photo wave SHA not signed")
    d=pd.read_csv(file,compression="gzip",dtype=str,keep_default_na=False)
    good=d.strict9.str.lower().eq("true")
    frac=pd.to_numeric(d.strict9_white,errors="coerce")
    require(len(d)==8192 and d.photo_id.is_unique and
            d.observation_id.is_unique and d.wave.eq(str(wave)).all() and
            int(good.sum())==entry["actual_numeric_original_ROI9_photo_success"] and
            frac.loc[good].between(0,1).all() and frac.loc[~good].isna().all() and
            d.wcvp_accepted_plant_name_id.ne("").all() and
            d.target_id.eq(d.year+"|"+d.target_id.str.split("|",regex=False).str[1]+
                           "|"+d.land_id).all(),
            "Original wave numerical photo ROI, site/year/accepted species evidence invalid")
    q=d[["photo_id","observation_id","target_id","land_id","year",
          "wcvp_accepted_plant_name_id"]].copy()
    q["year"]=q.year.astype(int)
    q["real_numeric_ROI9"]=good.astype(bool)
    q["source_photo_white_fraction"]=frac
    q["photo_source"]="tier08_new24576_wave"+str(wave)
    return q,{
        "wave":wave,
        "source_actual_original_photo_attempts":len(q),
        "true_numeric_original_flower_ROI9_success":int(good.sum()),
        "true_technical_unknown_photos":int((~good).sum()),
        "source_SHA256":digest
    }

def convergence_gate(stages,threshold=0.30,halfwidth=1.20,
                     minimum_islands=120,minimum_cells=3200,
                     source_disagreement_max=0.50):
    """Predeclared precision/stability rule independent of p-values.

    Three nonoverlapping fixed photo wave acquisitions are performed even if
    convergence occurs at wave2. Stability requires incremental coefficient
    shifts both <=0.30 for the last two 8192-photo additions, nominal CI
    halfwidth <=1.20 pp, >=120 true multisite islands, >=3200 multisite
    site-years and an alternate equal-source model within .50 pp of the
    original-photo pooled model. This is a *diagnostic*, not evidence of
    causality or hypothesis truth. Data already collected are not dropped.
    """
    if len(stages)<4:
        return {"stable":False,"reason":"need baseline and three signed 8192-photo waves"}
    last=stages[-1]
    diffs=[abs(stages[i]["island_equal_area_beta"]-
               stages[i-1]["island_equal_area_beta"]) for i in (-2,-1)]
    ci=(last["island_equal_area_CI_high"]-last["island_equal_area_CI_low"])/2
    nisl=last["multisite_islands"]
    nsite=last["multisite_site_years"]
    source=abs(last["island_equal_area_beta"]-
               last.get("equal_source_area_beta",float("nan")))
    met={
       "last_two_8192_photo_beta_shifts_both_at_most_0p30":
           bool(max(diffs)<=threshold),
       "nominal_CI_halfwidth_at_most_1p20":
           bool(ci<=halfwidth),
       "at_least_120_spatially_estimable_islands":
           bool(nisl>=minimum_islands),
       "at_least_3200_spatially_estimable_siteyears":
           bool(nsite>=minimum_cells),
       "equal_source_model_within_0p50":
           bool(np.isfinite(source) and source<=source_disagreement_max),
    }
    return {"stable":all(met.values()),"gates":met,
            "last_two_absolute_coefficient_changes":diffs,
            "nominal_CI_halfwidth_pp":ci,
            "between_source_weighting_estimand_difference_pp":source,
            "needs_additional_genuine_photo_waves":not all(met.values()),
            "p_value_does_not_control_stopping":True,
            "diagnostic_thresholds_are_operational_not_formal_power_or_equivalence":True}

def run(paths,coast,source06,source07,prior06,bulk06,pilot07,bulk07,
        tier08,tier08_receipt,large_receipt,waves,out):
    out=Path(out)
    require(not out.exists(),"No overwrite of signed ALAN photo convergence")
    canonical,wcvp=panel.load_frozen_photos(paths)
    old,_=historic.old_original_photos(canonical,wcvp,paths)
    new,_=historic.new_original_photos(
        source06,source07,prior06,bulk06,pilot07,bulk07)
    shared=set(old.photo_id)&set(new.photo_id)
    require(len(old)==50829 and len(new)==36449 and len(shared)==329,
            "Historical and tier06/07 source photo ID union changed")
    new_clean=new.loc[~new.photo_id.isin(shared)].copy()
    earlier=pd.concat([old,new_clean],ignore_index=True)
    require(earlier.photo_id.is_unique and len(earlier)==86949,
            "Old pooled source original photos duplicated or lost")
    first,_=prior08.load_real_tier08(tier08,tier08_receipt)
    first_clean,first_overlap=prior08.same_photo_reconciliation(earlier,first)
    require(len(first_clean)==903 and first_overlap[
         "tier08_original_photo_ids_already_in_prior"]==7,
            "Prior actual real tier08 910 image IDs/ROI9 have changed")
    base=pd.concat([earlier,first_clean],ignore_index=True)
    require(len(base)==87852 and int(base.real_numeric_ROI9.sum())==32111 and
            base.photo_id.is_unique,
            "Previous 87,852-photo / 32,111 genuinely measured numeric ROI9 source changed")
    loaded=[]
    receipt=[]
    for wave,file in enumerate(waves,start=1):
        q,r=verify_wave(file,large_receipt,wave)
        require(not set(base.photo_id).intersection(set(q.photo_id)),
                "New bulk photo source reuses prior historical/other new image ID")
        loaded.append(q);receipt.append(r)
    full=pd.concat([base,*loaded],ignore_index=True)
    require(full.photo_id.is_unique and
            len(full)==87852+24576,
            "Full 24,576 original tier08 distinct photo IDs not real independent photos")
    environment=fresh.original_sameyear_environment(paths)
    allframes=[base]
    for q in loaded:
        allframes.append(pd.concat([allframes[-1],q],ignore_index=True))
    island_frames=[]
    support=[]
    for i,p in enumerate(allframes):
        allsites,strict=historic.photo_cell_stats(p)
        a,qual=historic.environment_qualified(strict,environment)
        require(len(a)>0 and a.target_id.is_unique,
                "Cumulative true original photograph island site/year invalid")
        if i>0:
            require(set(island_frames[-1].target_id).issubset(set(a.target_id)),
                    "New genuine photos destroyed prior strict original source sites")
        island_frames.append(a)
        support.append(qual)
    require(len(island_frames[0])==2530 and
            island_frames[0].land_id.nunique()==193,
            "Previous signed 910-overlay strict island cell×year count changed")
    coastland=isolation.load_signed_land(coast)
    meta=coastland[["land_id","land_class","area_km2"]].copy()
    distances=isolation.isolated_mainland_distance_km(
        coastland,set(island_frames[-1].land_id))
    baseline_isl = set(island_frames[0].land_id)
    pre_multisite=area.join_area_and_select(island_frames[0],meta)[1]
    baseline_multisite_islands=set(pre_multisite.land_id)
    report_stages=[]
    for i,frame in enumerate(island_frames):
        label=f"cumulative_wave{i}_plus_{i*8192}_new_RGB_photos"
        fit=historic.fit_photo_arm(frame,meta,distances,label)
        require(fit["estimable"] and
                "joint_island_equal" in fit["models"],
                "True cumulative island contrasts could not be estimated")
        coef=fit["models"]["joint_island_equal"]["terms"][AREA_TERM]
        equal_source=None
        if i==len(island_frames)-1:
            tier08all=pd.concat([first_clean,*loaded],ignore_index=True)
            alternate,equal_support=prior08.three_arm_equal_source_photo_white(
                old,new_clean,tier08all,frame)
            equal_source=historic.fit_photo_arm(
                alternate,meta,distances,"final_equal_photo_source_within_cell")
        restricted=frame.loc[frame.land_id.isin(baseline_multisite_islands)].copy()
        common=historic.fit_photo_arm(
            restricted,meta,distances,"same_original_baseline_97_islands")
        # For robustness compare with the same old source geographic island set.
        stage={
          "checkpoint":i,
          "additional_original_RGB_attempts":i*8192,
          "all_distinct_original_photo_IDs":int(len(allframes[i])),
          "true_numeric_ROI9_original_photo_success":int(
              allframes[i].real_numeric_ROI9.sum()),
          "strict_true_island_site_years":len(frame),
          "strict_real_GSHHG_islands":int(frame.land_id.nunique()),
          "multisite_islands":fit["n_spatial_multisite_islands"],
          "multisite_site_years":fit["n_spatial_multisite_cells"],
          "island_equal_area_beta":coef["beta"],
          "island_equal_area_CI_low":coef["ci_low"],
          "island_equal_area_CI_high":coef["ci_high"],
          "island_equal_area_p":coef["p"],
          "cell_equal_area_beta":fit["models"]["joint_cell_equal"]["terms"][AREA_TERM]["beta"],
          "within_baseline_97_islands":
              common["models"]["joint_island_equal"]["terms"][AREA_TERM]
              if common["estimable"] else {"estimable":False},
          "leave_largest_photo_island_out_area":
              fit["models"].get("leave_largest_photo_island_out_equal_island",{}
                ).get("terms",{}).get(AREA_TERM),
          "photo_based_not_unbiased_botanical_colour":True,
        }
        if equal_source is not None:
            stage["equal_source_area_beta"]=equal_source["models"][
                "joint_island_equal"]["terms"][AREA_TERM]["beta"]
            stage["equal_source_sensitivity_support"]=equal_support
        report_stages.append(stage)
    require(abs(report_stages[0]["island_equal_area_beta"]-
                1.1059197451811407)<1e-5 and
            report_stages[0]["multisite_islands"]==97 and
            report_stages[0]["multisite_site_years"]==2419,
            "Previous truly observed 910-photo stage regression drift")
    gate=convergence_gate(report_stages)
    report={
        "schema":SCHEMA,
        "source_original_photo_image_measurement_git_blob":ROI_BLOB,
        "fixed_8192_photo_original_RGB_wave_SHA256":receipt,
        "all_original_image_attempts_from_clean_prior":len(base)+24576,
        "no_original_photo_ID_reused":True,
        "historical_source_numeric_obs_ID_new_UUID_full_independence_unverified":True,
        "source_sampling_conditional_island_area_light_pair_targeted":True,
        "no_original_RGB_copyrighted_image_export":True,
        "no_colour_unknown_imputation":True,
        "sameyear_raster_years":[2019,2021,2023,2024,2025],
        "baseline_97_multisite_islands_with_fixed_geography":len(baseline_multisite_islands),
        "checkpoints":report_stages,
        "fixed_a_priori_non_pvalue_stability_check":gate,
        "if_not_stable_use_remainder_of_56440_master_after_photo_ID_dedup":not gate["stable"],
        "ALAN_causality_or_native_wild_flower_prevalence_unverified":True,
    }
    out.mkdir(parents=True)
    file=out/"ALAN_tier08_24576_true_RGB_three_wave_island_area_light_convergence.json"
    file.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_TIER08_ADDITIONAL_24576_3_STAGE_REAL_FLOWER_ROI9_CONVERGENCE",
          json.dumps({"checkpoints":[{k:s[k] for k in (
             "checkpoint","additional_original_RGB_attempts",
             "all_distinct_original_photo_IDs",
             "true_numeric_ROI9_original_photo_success",
             "strict_true_island_site_years","strict_real_GSHHG_islands",
             "multisite_islands","multisite_site_years",
             "island_equal_area_beta","island_equal_area_CI_low",
             "island_equal_area_CI_high","island_equal_area_p")}
             for s in report_stages],"stability":gate},sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("canonical","wcvp",
                "geo2019","viirs2019","pop2019",
                "geo2021","viirs2021","pop2021",
                "geo2023","viirs2023","pop2023",
                "geo2024","cn2024","uc2024","uc2024receipt",
                "geo2025","geo2025receipt","viirs2025",
                "coast","source06","source07","prior06","bulk06","pilot07","bulk07",
                "tier08","tier08-receipt","large-receipt","wave1","wave2","wave3","out"):
        p.add_argument("--"+key.replace("_","-"),dest=key.replace("-","_"),type=Path,required=True)
    a=p.parse_args()
    paths={k:getattr(a,k) for k in (
       "canonical","wcvp","geo2019","viirs2019","pop2019",
       "geo2021","viirs2021","pop2021","geo2023","viirs2023","pop2023",
       "geo2024","cn2024","uc2024","uc2024receipt",
       "geo2025","geo2025receipt","viirs2025")}
    run(paths,a.coast,a.source06,a.source07,a.prior06,a.bulk06,a.pilot07,
        a.bulk07,a.tier08,a.tier08_receipt,a.large_receipt,
        [a.wave1,a.wave2,a.wave3],a.out)
if __name__=="__main__":main()
