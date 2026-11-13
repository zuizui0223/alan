#!/usr/bin/env python3
"""ALAN: actual day/night fruit/seed reproductive function versus independent
floral syndrome traits, and transfer beyond original study/plant family.

Unit: original Kendall/Nicholson 2025 published study × plant × field site ×
experiment year and season × fruit/seed outcome with TRUE night/day exclusion
arms. Light: ORIGINAL MD5-audited Li–Cao corrected prior 1992–year-1 PCNL
(mean annual log1p). The paired reproduction IS experimentally measured,
but flower traits/ALAN are externally observed and not manipulated.

Scientific questions: Is published plant bloom time associated with ACTUAL
nighttime reproductive benefit independent of scent, nectar, corolla shape,
white color and flower length? Does adding morphology improve generalization
when whole original studies or plant families are held out? Is an ALAN × bloom
effect consistently supported? Fail closed on missing/collinear trait support.
No source plant×study coordinate records or raw trait rows exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import t as student_t
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, LIGHT, require
)

AUTHOR_TRAIT_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
BASE_BLOOM_BETA=.21474990651679526
MAIN_FIXED=("exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)")
BLOOM_FORMULA=("night_share ~ night_bloom + "+LIGHT+" + "+MAIN_FIXED)
# Fully matched sample for ALL incremental comparisons. Source missingness
# is never automatically assigned the negative/diurnal class.
MORPH_NAMES=("odour_yes","nectar_yes","white","tube","log_flower_length")
BLOOM_NAMES=("night_bloom",)
CV_MODELS={
    "original_geography_year_metric_only":(),
    "bloom_time_only":BLOOM_NAMES,
    "other_flower_morphology_only":MORPH_NAMES,
    "bloom_plus_other_flower_morphology":BLOOM_NAMES+MORPH_NAMES
}


def original_traits(path:Path):
    raw=path.read_bytes()
    sha=hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
    require(sha==AUTHOR_TRAIT_BLOB,"Original independent botanical author trait blob changed")
    x=pd.read_csv(path,dtype=str,keep_default_na=False)
    require(len(x)==142 and x.plant_species.is_unique,"Author original 142 species records drift")
    expected={"plant_species","bloom_period_simple","odour","nectar","color",
              "flower_shape","flower_length_midpoint_mm","plant_family"}
    require(expected.issubset(x.columns),"Independent botanical trait columns unavailable")
    return x[list(expected)].rename(columns={"plant_species":"species"})


def qualified_source(original:pd.DataFrame,traits:pd.DataFrame):
    require(len(original)==301 and original.study_key.nunique()==118 and
            original.species.nunique()==123,
            "Original source from 30 validated prior-light rasters changed")
    full=original.merge(traits,on="species",how="left",validate="many_to_one")
    for c in ("odour","nectar","color","flower_shape",
              "bloom_period_simple","plant_family"):
        full[c]=full[c].replace({"NA":pd.NA,"":pd.NA,"nan":pd.NA})
    full["night_bloom"]=full.bloom_period_simple.eq("night")
    full["day_bloom"]=full.bloom_period_simple.eq("day")
    full["odour_yes"]=full.odour.eq("present")
    full["nectar_yes"]=full.nectar.eq("present")
    full["white"]=full.color.eq("white")
    full["tube"]=full.flower_shape.eq("tube")
    full["flower_length_mm"]=pd.to_numeric(full.flower_length_midpoint_mm,
                                           errors="coerce")
    full["log_flower_length"]=np.where(
        full.flower_length_mm.gt(0),np.log(full.flower_length_mm.clip(lower=.001)),np.nan)
    full["morph_bloom_known"]=full.bloom_period_simple.isin(("night","day")) & (
        full[["odour","nectar","color","flower_shape"]].notna().all(axis=1)) & (
        full.flower_length_mm.gt(0))
    full["night_bloom"]=full.night_bloom.astype(int)
    for c in ("odour_yes","nectar_yes","white","tube"):
        full[c]=full[c].astype("boolean")
    q=full.loc[full.morph_bloom_known].copy()
    for c in ("odour_yes","nectar_yes","white","tube"):
        require(q[c].notna().all(),"Incomplete direct floral trait erroneously negative")
        q[c]=q[c].astype(int)
    q["plant_family"]=q.plant_family.fillna("UNKNOWN_FAMILY")
    q["study_equal_weight"]=1/q.groupby("study_key").study_key.transform("size")
    bloom_only=full.loc[full.bloom_period_simple.isin(("night","day"))].copy()
    require((len(bloom_only),bloom_only.study_key.nunique(),
             bloom_only.species.nunique())==(196,73,74),
            "Original independent day/night bloom functional source drift")
    support={
       "original_field_study_site_year_reproductive_pairs":len(full),
       "author_bloom_day_night_source_pairs":len(bloom_only),
       "author_bloom_source_studies":int(bloom_only.study_key.nunique()),
       "author_bloom_source_species":int(bloom_only.species.nunique()),
       "all_other_morphology_and_bloom_direct_complete_pairs":len(q),
       "complete_original_publications":int(q.study_key.nunique()),
       "complete_botanical_species":int(q.species.nunique()),
       "complete_original_plant_families":int(q.plant_family.nunique()),
       "original_night_bloom_pairs_complete":int(q.night_bloom.sum()),
       "original_day_bloom_pairs_complete":int(len(q)-q.night_bloom.sum()),
       "source_study_keys_with_night_and_day_bloom_comparison":int(
           q.groupby("study_key").night_bloom.nunique().eq(2).sum()) if len(q) else 0,
       "source_families_with_night_and_day_bloom_comparison":int(
           q.groupby("plant_family").night_bloom.nunique().eq(2).sum()) if len(q) else 0,
       "source_trait_availability_by_original_variable":{
           col:{"known_source_pairs":int(full[col].notna().sum()),
                "source_studies":int(full.loc[full[col].notna(),"study_key"].nunique()),
                "source_species":int(full.loc[full[col].notna(),"species"].nunique())}
           for col in ("odour","nectar","color","flower_shape",
                       "flower_length_midpoint_mm","bloom_period_simple")},
       "known_flower_morphology_is_literature_sourced_not_individual_site_assay":True,
       "study_and_family_source_nonrandom":True
    }
    return q,bloom_only,support


def source_fit(d,formula,terms):
    """Study-equal WLS, study cluster and study×plant two-way covariance;
    reject rank-deficient cross-trait functional claims."""
    details={"source_pairs":len(d),
             "source_studies":int(d.study_key.nunique()),
             "source_species":int(d.species.nunique()),
             "source_families":int(d.plant_family.nunique())}
    if (len(d)<85 or d.study_key.nunique()<40 or
        d.night_bloom.nunique()!=2 or
        min(d.loc[d.night_bloom.eq(z),"study_key"].nunique() for z in (0,1))<18):
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_STUDY_TRAIT_SUPPORT",**details}
    s=d.copy()
    s["study_equal_weight"]=1/s.groupby("study_key").study_key.transform("size")
    r=smf.wls(formula,data=s,weights=s.study_equal_weight).fit(
       cov_type="cluster",cov_kwds={"groups":s.study_key},use_t=False)
    rank=int(np.linalg.matrix_rank(r.model.exog))
    if rank!=len(r.params):
        return {"status":"NOT_ESTIMABLE_EXACT_OR_NEAR_COLLINEARITY",
                "n_parameters":len(r.params),"rank":rank,**details}
    result={}
    for term in terms:
        if term not in r.params or not np.isfinite(r.bse[term]):
            result[term]={"status":"NOT_ESTIMABLE_TERM"};continue
        ci=r.conf_int().loc[term]
        result[term]={
            "coefficient_fraction":float(r.params[term]),
            "95CI_source_publication_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
            "nominal_p":float(r.pvalues[term]),
        }
    # Same source original study and botanical species are nonindependent.
    try:
        vcov,_,_=cov_cluster_2groups(
            r,pd.factorize(s.study_key)[0],pd.factorize(s.species)[0])
        df=min(s.study_key.nunique(),s.species.nunique())-1
        critical=float(student_t.ppf(.975,df))
        for term in terms:
            if term not in r.params:continue
            j=list(r.params.index).index(term)
            variance=float(vcov[j,j])
            if variance<=0 or not np.isfinite(variance):
                result[term]["two_way"]={"status":"NOT_ESTIMABLE_TWOWAY_VARIANCE"}
                continue
            se=float(np.sqrt(variance));v=float(r.params[term])
            result[term]["two_way"]={
                "95CI_student_t":[float(v-critical*se),float(v+critical*se)],
                "p_two_sided_t":float(2*student_t.sf(abs(v/se),df)),"df":df}
    except (ValueError,np.linalg.LinAlgError) as ex:
        for term in terms:
            if term in result:result[term]["two_way"]={
                "status":"NOT_ESTIMABLE_TWOWAY","reason":str(ex)[:150]}
    return {"status":"EXPLORATORY_DIRECT_REPRODUCTIVE_FUNCTION_TRAIT_MODEL",
            "formula":formula,"terms":result,
            "n_parameters":int(len(r.params)),"model_R2_descriptive":float(r.rsquared),
            **details,"not_a_causal_flower_trait_or_ALAN_response":True}


def fivefold_blocked_transfer(q:pd.DataFrame,groupkey:str):
    """True out-of-source group CV for four frozen models with identical cases.

    Five disjoint held-out original studies or BOTANICAL FAMILIES per split.
    Training-standardized ridge alpha=5, without hyperparameter search.
    No family identity/species/site encoded as a feature. All models have
    same outcome sample and splits; prevent near-duplicate published plant
    identities from being treated as independent family replication.
    """
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    if (len(q)<90 or q.study_key.nunique()<45 or
        q.plant_family.nunique()<12 or q.night_bloom.nunique()!=2):
        return {"status":"NOT_ESTIMABLE_EXPERIMENT_CV_STUDY_FAMILY_SUPPORT"}
    groups=q[groupkey].astype(str).to_numpy()
    if len(set(groups))<10:
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_CV_GROUPS"}
    baseline=["exp_year","lat","lat2","lon","lon2",
              "elevation_km",LIGHT,"metric_seed_set"]
    data=q.copy()
    data["metric_seed_set"]=data.metric.eq("seed set").astype(int)
    y=data.night_share.to_numpy(float)
    wg=1.0/data.groupby("study_key").study_key.transform("size").to_numpy(float)
    wg=wg/wg.sum()
    splits=list(GroupKFold(n_splits=5).split(data,y,groups=groups))
    scores={}
    all_pred={}
    for name,additional in CV_MODELS.items():
        cols=baseline+list(additional)
        X=data[cols].to_numpy(float)
        require(np.isfinite(X).all(),"CV features contain missing/imputed unknowns")
        pred=np.full(len(data),np.nan)
        for train,test in splits:
            scaler=StandardScaler()
            tr=scaler.fit_transform(X[train])
            ts=scaler.transform(X[test])
            fold_weight=wg[train]/np.mean(wg[train])
            fit=Ridge(alpha=5.0,fit_intercept=True)
            fit.fit(tr,y[train],sample_weight=fold_weight)
            pred[test]=np.clip(fit.predict(ts),0,1)
        require(np.isfinite(pred).all(),"CV incomplete prediction")
        err=pred-y
        mse=float(np.sum(wg*err*err))
        mae=float(np.sum(wg*np.abs(err)))
        scores[name]={
            "weighted_out_of_group_RMSE_fraction":float(np.sqrt(mse)),
            "weighted_out_of_group_MAE_fraction":mae,
            "source_groups_heldout_each_once":True}
        all_pred[name]=pred
    zero=scores["original_geography_year_metric_only"]["weighted_out_of_group_RMSE_fraction"]
    for a in scores:
        scores[a]["RMSE_change_pct_vs_same_fold_geography_baseline"]=float(
            100*(zero-scores[a]["weighted_out_of_group_RMSE_fraction"])/zero)
    return {
        "status":"EXPLORATORY_HELD_OUT_ORIGINAL_SOURCE_FUNCTIONAL_TRANSFER",
        "held_out_group":"original_publication" if groupkey=="study_key" else "botanical_family",
        "n_original_field_reproductive_pairs":len(q),
        "n_independent_original_studies":int(q.study_key.nunique()),
        "n_botanical_families":int(q.plant_family.nunique()),
        "fivefold_disjoint_group_count":5,
        "fixed_ridge_alpha":5.,
        "comparison_on_exact_same_original_valid_trait_cases":True,
        "source_held_out_prediction_clipped_0_to_1":True,
        "weighted_MAE_and_RMSE_original_night_reproductive_share_fraction":scores,
        "not_targeted_pollen_limitation_or_causal_intervention":True,
    }


def run(joined,receipt,trait_path,out):
    out=Path(out)
    if out.exists():raise FileExistsError("Refuse overwriting original source scientific result")
    original,orig_receipt=load_source(joined,receipt)
    q,bloom,support=qualified_source(original,original_traits(trait_path))
    source={"original_field_site_year_reproductive_pairs":len(original),
            "original_independent_publications":int(original.study_key.nunique()),
            "original_accepted_botanical_taxa":int(original.species.nunique()),
            "author_independent_trait_blob":AUTHOR_TRAIT_BLOB,
            "PCNL_original_30year_joined_sha256":orig_receipt["joined_file_sha256"]}
    # Original independently reproduced morphological time baseline on full
    # published night/day bloom evidence BEFORE complete-case reduction.
    main=bloom.copy()
    main["night_bloom"]=main.night_bloom.astype(int)
    require((len(main),main.study_key.nunique(),main.species.nunique())==
            (196,73,74),"Original night/day bloom study population changed")
    check=source_fit(main,BLOOM_FORMULA,["night_bloom",LIGHT])
    require(check["status"]=="EXPLORATORY_DIRECT_REPRODUCTIVE_FUNCTION_TRAIT_MODEL" and
            abs(check["terms"]["night_bloom"]["coefficient_fraction"]-BASE_BLOOM_BETA)<1e-7,
            "Original +21.475pp bloom and actual fruit/seed source coefficient not reconstructed")
    models={"baseline_original_bloom_only_full_support":check}
    complete_possible=(len(q)>=110 and q.study_key.nunique()>=45 and
                       q.night_bloom.nunique()==2)
    if complete_possible:
        formula=BLOOM_FORMULA
        plus=formula+" + "+" + ".join(MORPH_NAMES)
        morph=("night_share ~ "+LIGHT+" + "+MAIN_FIXED+
               " + "+" + ".join(MORPH_NAMES))
        models["same_source_completecase_bloom_only"]=source_fit(
            q,formula,["night_bloom",LIGHT])
        models["same_source_completecase_other_morphology_only"]=source_fit(
            q,morph,list(MORPH_NAMES)+[LIGHT])
        models["completecase_bloom_plus_other_flower_morphology"]=source_fit(
            q,plus,["night_bloom"]+list(MORPH_NAMES)+[LIGHT])
        models["completecase_family_FE_bloom_plus_morph"]=source_fit(
            q,plus+" + C(plant_family)",
            ["night_bloom"]+list(MORPH_NAMES)+[LIGHT])
        models["same_completecases_ALAN_x_bloom_separate"]=source_fit(
            q,plus.replace("night_bloom + "+LIGHT,
                           "night_bloom*"+LIGHT),
            ["night_bloom",LIGHT,"night_bloom:"+LIGHT])
    else:
        models["morphology_joint"]={
            "status":"NOT_ESTIMABLE_COMPLETECASE_STUDY_SUPPORT","source_support":support}
    cv={
      "fivefold_independent_published_study_holdout":fivefold_blocked_transfer(
          q,"study_key"),
      "fivefold_independent_plant_family_holdout":fivefold_blocked_transfer(
          q,"plant_family")
    }
    result={
        "schema":"ALAN_301_true_reproductive_diel_function_multi_floral_trait_study_family_transfer_v1",
        "status":"AUTHENTIC_PAIR_DAYNIGHT_FIELD_REPRODUCTION_FLORAL_TRAIT_EXPLORATION",
        "source":source,"direct_trait_coverage":support,
        "models":models,
        "out_of_source_functional_prediction_transfer":cv,
        "strict_scientific_limits":{
           "original_day_night_paired_breeding_metric_really_measured":True,
           "source_floral_traits_literature_not_same_field_assay":True,
           "bloom_flower_trait_effect_NOT_caused_by_experimental_manipulation":True,
           "ALAN_previous_PCNL_observational_NOT_at_flower_lux":True,
           "daynight_function_NOT_moth_specific":True,
           "family_study_holdout_prediction_NOT_causal_validation":True,
           "same_study_night_day_bloom_phenotype_common_support_limited":True,
           "not_island_only_and_not_measured_population_recovery":True,
           "multiple_posthoc_original_trait_models_not_preregistered":True,
           "no_individual_source_taxon_field_coordinates_or_restricted_trait_rows_exported":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False,default=str)+"\n")
    print("ALAN_ORIGINAL_FLORAL_TRAIT_DIEL_FUNCTION_MODEL",
          json.dumps({"source":source,"coverage":support,
                      "models":{k:{"status":v.get("status"),
                                   "terms":v.get("terms")} for k,v in models.items()},
                      "transfer":cv},ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("joined","receipt","trait","out"):
        p.add_argument("--"+name,required=True,type=Path)
    x=p.parse_args()
    run(x.joined,x.receipt,x.trait,x.out)
