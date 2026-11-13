#!/usr/bin/env python3
"""Original field day/night fruit/seed function x sunlight exposure and daily
temperature range (DTR) as independent conditions of flowering-time function.

Inputs:
- 301 original Kendall/Nicholson (2025) paired real day/night reproductive
  study×species×site×year×season×metric outcomes from previous verified source;
- authentic Li–Cao corrected 1992–2021 previous-year PCNL exposure MD5 joined;
- independently published author botanicals (bloom night/day);
- original author publication study-site DvN_Site_Daylength.csv and
  DvN_Site_DailyTemperatureRange.csv, matched by exact study, rounded original
  field location, true year and original start/end field MONTHS.

Original article already evaluates daylength and DTR separately. Novel here is
conditional identification of flowering-time correspondence and cross FAMILY/
study+species transfer with nightlight in the same frozen original paired sample.

Do NOT claim manipulated nightlight, effective moth pollen, evolutionary
change, or daylength/DTR as novel independent literature findings. No raw
study×plant×location or direct trait rows in results.
"""
from __future__ import annotations
import argparse
import io
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold
from scipy.stats import t as student_t
import statsmodels.formula.api as smf

from scripts.analysis.alan_real_DvN_flower_function_family_LOFO_purgedCV_20261113 import (
    source_data, BASE_COLS, ALPHA, N_PAIRS
)
from scripts.analysis.alan_flower_multitraits_reproductive_diel_transfer_20261113 import (
    LIGHT, source_fit, BLOOM_FORMULA
)
from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import require

BLOB_DAY="97ace3708bc356d87ce4fcd2e00fc4f1987be949"
BLOB_DTR="f2f16164e96570e6acd5e6ab975295d2bcc6bd1f"
ENV_KEYS=["study_key","geo_lat5","geo_lon5","year_exp","month_start","month_end"]
ENV_COLS=("Daylength","DTR")
MODELS={
    "original_field_geo_year_light":(),
    "bloom_only":("night_bloom",),
    "environment_only":("photoperiod_h","DTR_C"),
    "bloom_plus_environment":("night_bloom","photoperiod_h","DTR_C"),
    "bloom_x_environment":("night_bloom","photoperiod_h","DTR_C",
                            "night_bloom_x_DTR","night_bloom_x_photoperiod"),
}


def sha_git_blob(path:Path):
    b=path.read_bytes()
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()


def normal_study(s:pd.Series):
    """Only two SOURCE-DOCUMENTED mojibake repairs, not generic fuzzy match."""
    return (s.fillna("").astype(str).str.strip()
         .str.replace("JÃ¼rgens_2014","Jürgens_2014",regex=False)
         .str.replace("Larrea-AlcÃ¡zar_2011","Larrea-Alcázar_2011",regex=False)
         .str.replace("JÃÂ¼rgens_2014","Jürgens_2014",regex=False))


def source_environment_file(path:Path,expected_blob:str,value:str):
    require(sha_git_blob(path)==expected_blob,
            "Independent original study site climatic variable file Git SHA changed")
    # Original author data includes pre-Unicode R CSVs with Latin-1 bytes
    # (e.g., Jürgens_2014). Preserve exact raw Git blob verification,
    # decode to Unicode explicitly before study-ID key reconciliation.
    raw=path.read_bytes()
    try:
        text_source=raw.decode("utf-8-sig")
        source_encoding="utf-8-sig"
    except UnicodeDecodeError:
        text_source=raw.decode("latin-1")
        source_encoding="latin-1"
    z=pd.read_csv(io.StringIO(text_source),dtype=str,keep_default_na=False)
    names={"study_ID","coordinates_lat","coordinates_lon","year_start",
           "month_start","month_end",value}
    require(names.issubset(z.columns),"Original author environment columns drift")
    require(len(z)==218,"Original author environmental 219 site season/year rows drift")
    z["study_key"]=normal_study(z.study_ID)
    z["geo_lat5"]=pd.to_numeric(z.coordinates_lat,errors="coerce").round(5)
    z["geo_lon5"]=pd.to_numeric(z.coordinates_lon,errors="coerce").round(5)
    z["year_exp"]=pd.to_numeric(z.year_start,errors="coerce")
    z["month_start"]=z.month_start.astype(str).str.strip()
    z["month_end"]=z.month_end.astype(str).str.strip()
    z[value]=pd.to_numeric(z[value],errors="coerce")
    # Distinct same-site/year/season rows must not be arbitrarily averaged or
    # randomly picked. Conflicted genuine original metadata stay UNKNOWN.
    valid=z[ENV_KEYS].notna().all(axis=1)&np.isfinite(z[value])
    cols=ENV_KEYS+[value]
    z=z.loc[valid,cols].copy()
    groups=z.groupby(ENV_KEYS,dropna=False,observed=True)[value]
    report=groups.agg(["size","min","max"]).reset_index()
    duplicate_count=int(report["size"].gt(1).sum())
    conflict_count=int((report["max"]-report["min"]).abs().gt(1e-6).sum())
    good=report.loc[(report["max"]-report["min"]).abs().le(1e-6),ENV_KEYS+["min"]]
    good=good.rename(columns={"min":value})
    require(not good.duplicated(ENV_KEYS).any(),"Ambiguous original study weather row join")
    return good,{
        "original_source_file_blob":expected_blob,
        "original_source_text_encoding":source_encoding,
        "original_site_season_year_rows":218,
        "numeric_season_rows":int(len(z)),
        "original_distinct_site_season_year_keys":int(len(report)),
        "keys_with_repeated_rows":duplicate_count,
        "keys_with_conflicting_original_values_EXCLUDED":conflict_count,
        "consistent_original_weather_keys_eligible":int(len(good)),
    }


def join_environment(q:pd.DataFrame,day:Path,dtr:Path):
    ph,ph_meta=source_environment_file(day,BLOB_DAY,"Daylength")
    temp,temp_meta=source_environment_file(dtr,BLOB_DTR,"DTR")
    original=q.copy()
    original["geo_lat5"]=original.lat.round(5)
    original["geo_lon5"]=original.lon.round(5)
    original["study_key"]=normal_study(original.study_key)
    original["month_start"]=original.month_start.astype(str).str.strip()
    original["month_end"]=original.month_end.astype(str).str.strip()
    x=original.merge(ph,on=ENV_KEYS,how="left",validate="many_to_one")
    x=x.merge(temp,on=ENV_KEYS,how="left",validate="many_to_one")
    require(len(x)==len(q),"Environmental join expanded original site-year records")
    source_known=x.Daylength.notna() & x.DTR.notna()
    x["photoperiod_h"]=pd.to_numeric(x.Daylength,errors="coerce")
    x["DTR_C"]=pd.to_numeric(x.DTR,errors="coerce")
    x["night_bloom_x_DTR"]=x.night_bloom*x.DTR_C
    x["night_bloom_x_photoperiod"]=x.night_bloom*x.photoperiod_h
    # Source shares honest partial and both climatological coverage.
    q2=x.loc[source_known & x.photoperiod_h.between(4,23) &
             x.DTR_C.between(0,35)].copy()
    return q2, {
      "source_original_bloom_complete_pairs":len(q),
      "original_source_published_studies":int(q.study_key.nunique()),
      "original_botanical_families":int(q.plant_family.nunique()),
      "daylength_author_source":ph_meta,
      "daily_temp_range_author_source":temp_meta,
      "matching_requires_original_study_site_year_and_month_season":True,
      "original_reproductive_pairs_with_published_Daylength":int(x.Daylength.notna().sum()),
      "original_reproductive_pairs_with_published_DTR":int(x.DTR.notna().sum()),
      "original_reproductive_pairs_with_both_same_site_year_season":int(source_known.sum()),
      "joint_eligible_pairs":len(q2),
      "joint_independent_studies":int(q2.study_key.nunique()),
      "joint_botanical_species":int(q2.species.nunique()),
      "joint_botanical_families":int(q2.plant_family.nunique()),
      "joint_independent_nightbloom_studies":int(q2.loc[q2.night_bloom.eq(1),"study_key"].nunique()),
      "joint_independent_daybloom_studies":int(q2.loc[q2.night_bloom.eq(0),"study_key"].nunique()),
      "joint_families_with_both_bloom_types":int(q2.groupby("plant_family").night_bloom.nunique().eq(2).sum()) if len(q2) else 0,
      "photoperiod_hours_range":[float(q2.photoperiod_h.min()),float(q2.photoperiod_h.max())] if len(q2) else None,
      "DTR_C_range":[float(q2.DTR_C.min()),float(q2.DTR_C.max())] if len(q2) else None,
      "original_study_climate_source_year_months_not_general_future_weather":True,
    }


def conditional_models(s:pd.DataFrame):
    if len(s)<85 or s.study_key.nunique()<40 or s.plant_family.nunique()<12 or s.night_bloom.nunique()!=2:
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_REPRODUCTIVE_SOURCE_ENVIRONMENT_OVERLAP"}
    # Center at the ORIGINAL measured source study-season means, with
    # no interpretation as environmental manipulation.
    a=s.copy()
    a["photoperiod_center"]=a.photoperiod_h-float(s.photoperiod_h.mean())
    a["DTR_center"]=a.DTR_C-float(s.DTR_C.mean())
    frozen=("exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric) + "+LIGHT)
    form={
      "bloom_only_comparable_source":"night_share ~ night_bloom + "+frozen,
      "bloom_daylength_DTR_main":"night_share ~ night_bloom + photoperiod_center + DTR_center + "+frozen,
      "nightbloom_x_DTR":"night_share ~ night_bloom*DTR_center + photoperiod_center + "+frozen,
      "nightbloom_x_daylength":"night_share ~ night_bloom*photoperiod_center + DTR_center + "+frozen,
      "nightbloom_x_both_natural_light_temp":"night_share ~ night_bloom*(photoperiod_center + DTR_center) + "+frozen,
      "nightbloom_natural_env_family_FE":"night_share ~ night_bloom + photoperiod_center + DTR_center + "+frozen+" + C(plant_family)",
      "nightbloom_x_historical_ALAN_natural_controls":"night_share ~ night_bloom*"+LIGHT+
          " + photoperiod_center + DTR_center + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)"
    }
    terms={
      "bloom_only_comparable_source":["night_bloom"],
      "bloom_daylength_DTR_main":["night_bloom","photoperiod_center","DTR_center"],
      "nightbloom_x_DTR":["night_bloom","night_bloom:DTR_center"],
      "nightbloom_x_daylength":["night_bloom","night_bloom:photoperiod_center"],
      "nightbloom_x_both_natural_light_temp":["night_bloom","night_bloom:DTR_center","night_bloom:photoperiod_center"],
      "nightbloom_natural_env_family_FE":["night_bloom","photoperiod_center","DTR_center"],
      "nightbloom_x_historical_ALAN_natural_controls":["night_bloom",LIGHT,"night_bloom:"+LIGHT]
    }
    return {k:source_fit(a,formula,terms[k]) for k,formula in form.items()}


def grouped_transfer(s:pd.DataFrame,group:str):
    """Original frozen ridge model family/publication block heldout
    comparison on EXACT same data rows and original outcome scales.
    """
    if len(s)<95 or s.study_key.nunique()<45 or s.plant_family.nunique()<15:
        return {"status":"NOT_ESTIMABLE_ORIGINAL_STUDY_FAMILY_ENVIRONMENT_CV_SUPPORT",
                "cases":len(s),"studies":int(s.study_key.nunique()),
                "families":int(s.plant_family.nunique())}
    d=s.copy()
    d["metric_seed_set"]=d.metric.eq("seed set").astype(int)
    predictors=["exp_year","lat","lat2","lon","lon2","elevation_km",LIGHT,"metric_seed_set"]
    variants={
       "original_context":[],
       "context_plus_night_bloom":["night_bloom"],
       "context_plus_bloom_daylength_only":["night_bloom","photoperiod_h"],
       "context_plus_bloom_DTR_only":["night_bloom","DTR_C"],
       "context_plus_natural_daylength_DTR":["photoperiod_h","DTR_C"],
       "context_plus_night_bloom_daylength_DTR":["night_bloom","photoperiod_h","DTR_C"],
       "context_plus_bloom_x_DTR_daylength":["night_bloom","photoperiod_h","DTR_C",
                                             "night_bloom_x_DTR","night_bloom_x_photoperiod"]
    }
    xx={k:d[predictors+v].to_numpy(float) for k,v in variants.items()}
    require(all(np.isfinite(a).all() for a in xx.values()),
            "Source missing environmental data imputed or treated as zero")
    y=d.night_share.to_numpy(float)
    foldgroups=d[group].astype(str).to_numpy()
    if len(set(foldgroups))<10:return {"status":"NOT_ESTIMABLE_FEWER_THAN_TEN_BLOCKS"}
    if group=="plant_family":
        groups=[(np.flatnonzero(foldgroups!=c),np.flatnonzero(foldgroups==c))
                for c in sorted(set(foldgroups))]
        weights=1/d.groupby("plant_family").plant_family.transform("size").to_numpy(float)
    else:
        groups=list(GroupKFold(n_splits=5).split(d,y,groups=foldgroups))
        weights=1/d.groupby("study_key").study_key.transform("size").to_numpy(float)
    preds={k:np.full(len(d),np.nan) for k in variants}
    for train,test in groups:
        if len(train)<70 or d.iloc[train].night_bloom.nunique()!=2 or d.iloc[train].study_key.nunique()<35:
            return {"status":"NOT_ESTIMABLE_HELDOUT_TRAINING_STUDY_SUPPORT",
                    "source_cases":len(s),"heldout_group":group}
        trweights=1/d.iloc[train].groupby("study_key").study_key.transform("size").to_numpy(float)
        for name,arr in xx.items():
            scale=StandardScaler()
            train_x=scale.fit_transform(arr[train])
            test_x=scale.transform(arr[test])
            m=Ridge(alpha=5.0)
            m.fit(train_x,y[train],sample_weight=trweights/trweights.mean())
            preds[name][test]=np.clip(m.predict(test_x),0,1)
    result={}
    w=weights/weights.sum()
    for k,p in preds.items():
        require(np.isfinite(p).all(),"Original environmental source group test omitted")
        mse=float(np.sum(w*(p-y)**2))
        result[k]={"RMSE_fraction":float(np.sqrt(mse)),
                   "MAE_fraction":float(np.sum(w*np.abs(p-y)))}
    base=result["original_context"]["RMSE_fraction"]
    for val in result.values():
        val["RMSE_improvement_pct_against_original_context"]=100*(base-val["RMSE_fraction"])/base
    # Aggregate family-specified paired prediction errors, not original
    # individual biological source data. Even with LOFO, the source is a
    # nonrandom highly imbalanced 28-family literature sample.
    family_bootstrap=None
    if group=="plant_family":
        fam_names=sorted(d.plant_family.unique())
        source_family_MSE={}
        for variant,pred in preds.items():
            source_family_MSE[variant]=np.array([
                float(np.mean((y[d.plant_family.to_numpy()==fam]-
                               pred[d.plant_family.to_numpy()==fam])**2))
                for fam in fam_names])
        focal=source_family_MSE["context_plus_night_bloom"]
        combined=source_family_MSE["context_plus_night_bloom_daylength_DTR"]
        baseline=source_family_MSE["original_context"]
        rng=np.random.default_rng(20261113)
        repeat=rng.integers(0,len(fam_names),size=(3000,len(fam_names)))
        base_rmse=np.sqrt(np.mean(focal[repeat],axis=1))
        combined_rmse=np.sqrt(np.mean(combined[repeat],axis=1))
        gain=100*(base_rmse-combined_rmse)/base_rmse
        family_bootstrap={
            "family_independent_source_groups":len(fam_names),
            "family_count_where_bloom_plus_climate_lower_MSE_than_bloom_only":
                int((combined<focal).sum()),
            "family_count_where_bloom_plus_climate_lower_MSE_than_original_context":
                int((combined<baseline).sum()),
            "family_equal_RMSE_incremental_gain_pct_climate_OVER_bloom_only":
                float(100*(np.sqrt(np.mean(focal))-np.sqrt(np.mean(combined)))/np.sqrt(np.mean(focal))),
            "family_block_fixed_seed_descriptive_95pct_gain_interval_pct":
                [float(np.quantile(gain,.025)),float(np.quantile(gain,.975))],
            "family_block_bootstrap_draws":3000,
            "not_universal_random_botanical_family_replicates":True
        }
    return {"status":"EXPLORATORY_ORIGINAL_SOURCE_INDEPENDENT_GROUP_TRANSFER",
            "source_holdout_group":group,"n_source_pairs":len(d),
            "original_published_studies":int(d.study_key.nunique()),
            "original_plant_families":int(d.plant_family.nunique()),
            "n_disjoint_heldout_groups":len(groups),
            "all_models_identical_eligible_original_study_site_year_pairs":True,
            "source_selection_for_environment_independent_of_reproductive_outcomes":True,
            "unobserved_weather_never_imputed_from_future_month":True,
            "environment_incremental_family_bootstrap":family_bootstrap,
            "metrics":result}



def study_plus_species_purged_transfer(s:pd.DataFrame):
    """True study GroupKFold, remove all heldout test PLANT SPECIES from
    training. No paper IDs nor species values used as predictive features.
    """
    if len(s)<110 or s.study_key.nunique()<48 or s.plant_family.nunique()<15:
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_STUDY_SPECIES_ENVIRONMENT_SUPPORT"}
    d=s.copy()
    d["metric_seed_set"]=d.metric.eq("seed set").astype(int)
    cols=["exp_year","lat","lat2","lon","lon2","elevation_km",LIGHT,"metric_seed_set"]
    variants={
       "context":[],
       "bloom":["night_bloom"],
       "bloom_and_daylength_DTR":["night_bloom","photoperiod_h","DTR_C"],
       "bloom_x_natural_conditions":["night_bloom","photoperiod_h","DTR_C",
                "night_bloom_x_DTR","night_bloom_x_photoperiod"]
    }
    xx={key:d[cols+v].to_numpy(float) for key,v in variants.items()}
    require(all(np.isfinite(x).all() for x in xx.values()),"Incomplete original study environmental features")
    y=d.night_share.to_numpy(float)
    tests=list(GroupKFold(n_splits=5).split(d,y,groups=d.study_key))
    prediction={k:np.full(len(d),np.nan) for k in variants}
    deleted=[]
    for train,test in tests:
        forbidden=set(d.iloc[test].species)
        retained=train[~d.iloc[train].species.isin(forbidden).to_numpy()]
        deleted.append(int(len(train)-len(retained)))
        if (len(retained)<85 or d.iloc[retained].study_key.nunique()<38 or
            d.iloc[retained].night_bloom.nunique()!=2):
            return {"status":"NOT_ESTIMABLE_AFTER_STUDY_AND_SPECIES_PURGE",
                    "source_samples":len(d),"purged_training_rows_per_fold":deleted}
        w=1/d.iloc[retained].groupby("study_key").study_key.transform("size").to_numpy(float)
        for name,arr in xx.items():
            st=StandardScaler()
            trainxx=st.fit_transform(arr[retained])
            testxx=st.transform(arr[test])
            model=Ridge(alpha=5.0)
            model.fit(trainxx,y[retained],sample_weight=w/w.mean())
            prediction[name][test]=np.clip(model.predict(testxx),0,1)
    weight=1/d.groupby("study_key").study_key.transform("size").to_numpy(float)
    weight=weight/weight.sum()
    out={}
    for name,p in prediction.items():
        require(np.isfinite(p).all(),"Unpredicted independent original source test fold")
        out[name]={
            "RMSE_fraction":float(np.sqrt(np.sum(weight*(p-y)**2))),
            "MAE_fraction":float(np.sum(weight*np.abs(p-y)))}
    for name,a in out.items():
        a["RMSE_improvement_pct_vs_context"]=float(
            100*(out["context"]["RMSE_fraction"]-a["RMSE_fraction"])/out["context"]["RMSE_fraction"])
        a["RMSE_change_pct_vs_bloom_only"]=float(
            100*(out["bloom"]["RMSE_fraction"]-a["RMSE_fraction"])/out["bloom"]["RMSE_fraction"])
    return {"status":"EXPLORATORY_STUDY_HELDOUT_ALL_TEST_SPECIES_PURGED",
            "n_original_pairs":len(d),
            "n_heldout_study_folds":5,
            "n_original_botanical_species":int(d.species.nunique()),
            "n_original_independent_studies":int(d.study_key.nunique()),
            "study_rows_excluded_for_species_overlap_per_fold":deleted,
            "training_study_weights_recomputed_after_species_exclusion":True,
            "metrics":out,
            "study_source_species_and_georeferenced_locations_private":True}


def run(joined,receipt,author_trait,day,dtr,out):
    out=Path(out)
    if out.exists():raise FileExistsError("Refuse overwriting original scientific report")
    original,rec,trcov,original_fit=source_data(joined,receipt,author_trait)
    alljoint,support=join_environment(original,day,dtr)
    env_models=conditional_models(alljoint)
    transfers={
       "leave_one_whole_botanical_family_out":grouped_transfer(alljoint,"plant_family"),
       "held_out_published_study_fivefold":grouped_transfer(alljoint,"study_key"),
       "heldout_study_and_all_shared_species_purged":study_plus_species_purged_transfer(alljoint)
    }
    report={
        "schema":"ALAN_real301_original_DvN_reproductive_daylength_DTR_bloom_priorPCNL_environment_v1",
        "status":"AUTHENTIC_INDEPENDENT_PUBLISHED_ENVIRONMENT_PAIRS_GLOBAL_FIELD_BLOOM_ANALYSIS",
        "base":{
            "source_301_original_daynight_reproduction_study_site_year_original_PCNL_join_SHA":rec["joined_file_sha256"],
            "independent_botanical_original_bloom_regression_unchanged":original_fit,
            "original_complete_flower_bloom_pairs":len(original),
            "true_PCNL_raster_annual_years":[1992,2021],
        },
        "environmental_source_support":support,
        "source_environment_and_nightbloom_effect_models":env_models,
        "out_of_study_and_family_true_reproductive_transfer":transfers,
        "conclusion_limits":{
            "article_already_examined_daylength_DTR_separately":True,
            "source_date_weather_joins_original_study_site_year_and_season_only":True,
            "night_bloom_and_temperature_not_manipulated_field_treatment":True,
            "light_exposure_is_not_randomized_and_not_at_flower_lux":True,
            "no_direct_moth_pollinator_identity_at_experiment":True,
            "family_coverage_uneven_and_literature_selected":True,
            "no_island_only_or_colonization_claim":True,
            "multiple_condition_terms_are_exploratory":True,
            "no_raw_field_study_species_site_location_or_trait_rows_exported":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_TRUE_DAYNIGHT_DTR_DAYLENGTH_BLOOM_PRIORLIGHT",
        json.dumps({"support":support,
                    "models":{k:{"status":v.get("status"),"terms":v.get("terms")} for k,v in env_models.items()},
                    "transfer":transfers},ensure_ascii=False),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("joined","receipt","trait","day","dtr","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    run(a.joined,a.receipt,a.trait,a.day,a.dtr,a.out)
