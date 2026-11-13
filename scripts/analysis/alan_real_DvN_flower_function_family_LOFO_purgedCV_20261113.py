#!/usr/bin/env python3
"""Diagnose how independently measured floral traits transfer to unseen plant
families and species from REAL day/night fruit/seed field experiments.

Original source: 301 source-locked paired night-only/day-only reproductive
outcomes from 118 published studies and 123 plants, independently original
MD5-verified 1992–2021 prior-PCNL at actual experiment site/year. Restrict
to identically measured 196 complete-case reproductive responses, 73 studies,
74 plants, 28 families with externally documented day or night bloom, scent,
nectar, white color, tube/non-tube and measured positive flower length.

Never treat nondocumented blooms as 'day' nor infer moth-only service. No raw
species, site IDs, named study or restricted trait-by-geography rows exported.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, require
)
from scripts.analysis.alan_flower_multitraits_reproductive_diel_transfer_20261113 import (
    qualified_source, original_traits, source_fit, BLOOM_FORMULA,
    BASE_BLOOM_BETA, CV_MODELS, LIGHT, MORPH_NAMES
)

N_PAIRS=196
N_STUDIES=73
N_SPECIES=74
N_FAMILIES=28
ALPHA=5.0
BASE_COLS=("exp_year","lat","lat2","lon","lon2",
           "elevation_km",LIGHT,"metric_seed_set")
BASE_KEY="original_geography_year_metric_only"


def source_data(joined:Path,receipt:Path,traits:Path):
    s,signature=load_source(joined,receipt)
    q,bloom,summary=qualified_source(s,original_traits(traits))
    require((len(q),q.study_key.nunique(),q.species.nunique(),
             q.plant_family.nunique())==
            (N_PAIRS,N_STUDIES,N_SPECIES,N_FAMILIES),
            "Original 196/73/74/28 complete trait and experimental corpus drift")
    base=source_fit(q,BLOOM_FORMULA,["night_bloom"])
    require(base["status"]=="EXPLORATORY_DIRECT_REPRODUCTIVE_FUNCTION_TRAIT_MODEL" and
            abs(base["terms"]["night_bloom"]["coefficient_fraction"]-BASE_BLOOM_BETA)<1e-7,
            "Parent +21.5pp original actual reproductive functional comparison drift")
    q=q.copy()
    q["genus"]=q.species.str.split().str[0]
    q["metric_seed_set"]=q.metric.eq("seed set").astype(int)
    q["source_study_weight"]=1/q.groupby("study_key").study_key.transform("size")
    q["source_family_weight"]=1/q.groupby("plant_family").plant_family.transform("size")
    return q,signature,summary,base


def lofo_models(q:pd.DataFrame, group="plant_family"):
    """Every botanical family/genus omitted completely from training.

    No hyperparameter search: same original complete cases, identical held-out
    groups, controls and source study-weight for all four candidate trait sets.
    """
    require(group in ("plant_family","genus"),"Unexpected source holdout dimension")
    require(len(q)>=150 and q.plant_family.nunique()>=20,"Too few original source families")
    y=q.night_share.to_numpy(float)
    base_w=q.source_study_weight.to_numpy(float)
    family_w=q.source_family_weight.to_numpy(float)
    feature_arrays={}
    for name, extra in CV_MODELS.items():
        cols=list(BASE_COLS)+list(extra)
        feature_arrays[name]=q[cols].to_numpy(float)
        require(np.isfinite(feature_arrays[name]).all(),
                "Missing original trait mistakenly coded as false or imputed")
    families=q.plant_family.astype(str).to_numpy()
    groups=q[group].astype(str).to_numpy()
    names=list(CV_MODELS)
    n_groups=int(q[group].nunique())
    require(n_groups>=20,"Insufficient independent leave-one-group-out replication")
    pred={k:np.full(len(q),np.nan) for k in names}
    total_per_group=[]
    for category in sorted(np.unique(groups)):
        test=np.flatnonzero(groups==category)
        train=np.flatnonzero(groups!=category)
        if (len(train)<100 or q.iloc[train].study_key.nunique()<40 or
            q.iloc[train].plant_family.nunique()<16 or
            q.iloc[train].night_bloom.nunique()!=2):
            return {"status":"NOT_ESTIMABLE_LEAVEOUT_FAMILY_UNBALANCED_TRAINING",
                    "group_type":group,"blocked_test_groups":n_groups}
        for name,xx in feature_arrays.items():
            scale=StandardScaler()
            tr=scale.fit_transform(xx[train])
            te=scale.transform(xx[test])
            # Preserve original published study equal weighting; independent
            # studies' test records never enter the scaling/fitting.
            tw=base_w[train]/np.mean(base_w[train])
            fit=Ridge(alpha=ALPHA,fit_intercept=True)
            fit.fit(tr,y[train],sample_weight=tw)
            pred[name][test]=np.clip(fit.predict(te),0.,1.)
        total_per_group.append({
            "n_original_pairs":len(test),
            "original_independent_published_studies":int(q.iloc[test].study_key.nunique()),
            "night_bloom_fraction":float(q.iloc[test].night_bloom.mean()),
            "has_both_bloom_classes":bool(q.iloc[test].night_bloom.nunique()==2),
            "group_name_private":True
        })
    for name in pred:require(np.isfinite(pred[name]).all(),
                              "Incomplete original held-out functional transfer")
    # Original independent unit aggregation for family bootstrap.
    family_ids=sorted(np.unique(families))
    fam_score={}
    for family in family_ids:
        ix=np.flatnonzero(families==family)
        w=base_w[ix]/base_w[ix].sum()
        fam_score[family]={
            "n_original_field_outcomes":len(ix),
            "n_published_studies":int(q.iloc[ix].study_key.nunique()),
            "n_botanical_species":int(q.iloc[ix].species.nunique()),
            "both_bloom_classes":bool(q.iloc[ix].night_bloom.nunique()==2),
            "errors":{name:float(np.sum(w*(y[ix]-pred[name][ix])**2))
                      for name in names},
        }
    result={}
    for setting,weights in (("study_equal",base_w),
                            ("botanical_family_equal",family_w)):
        norm=weights/weights.sum()
        metrics={}
        for name in names:
            residual=y-pred[name]
            mse=float(np.sum(norm*residual**2))
            metrics[name]={
                "holdout_RMSE_nightshare_fraction":float(np.sqrt(mse)),
                "holdout_MAE_nightshare_fraction":float(np.sum(norm*np.abs(residual))),
            }
        benchmark=metrics[BASE_KEY]["holdout_RMSE_nightshare_fraction"]
        for name in names:
            metrics[name]["RMSE_improvement_percent_vs_geography_only"]=float(
                100*(benchmark-metrics[name]["holdout_RMSE_nightshare_fraction"])/benchmark)
        result[setting]=metrics

    # Family-level difference of squared errors is the scientifically relevant
    # transfer comparison, NOT in-sample significance of coefficients.
    delta={}
    for name in names:
        diffs=np.array([fam_score[f]["errors"][BASE_KEY]-
                        fam_score[f]["errors"][name] for f in family_ids])
        better=float(np.mean(diffs>0))
        mixed=np.array([fam_score[f]["both_bloom_classes"] for f in family_ids])
        # 2,500 fixed deterministic independent-FAMILY resamples. These
        # intervals reflect the variability of THIS selected family set, not
        # population- or phylogeny-wide confirmatory inference.
        rng=np.random.default_rng(20261113)
        baseline=np.array([fam_score[f]["errors"][BASE_KEY] for f in family_ids])
        candidate=np.array([fam_score[f]["errors"][name] for f in family_ids])
        boots=rng.integers(0,len(family_ids),size=(2500,len(family_ids)))
        change=100*(np.sqrt(np.mean(baseline[boots],axis=1))-
                    np.sqrt(np.mean(candidate[boots],axis=1))) / np.sqrt(
                    np.mean(baseline[boots],axis=1))
        delta[name]={
            "source_family_count":len(family_ids),
            "family_level_mse_gain_mean":float(np.mean(diffs)),
            "fraction_families_with_lower_MSE_than_geography_baseline":better,
            "source_families_with_both_bloom_categories":int(mixed.sum()),
            "mean_family_MSE_gain_in_within_family_mixed_bloom_sources":
                float(np.mean(diffs[mixed])) if mixed.any() else None,
            "family_equal_RMSE_improvement_bootstrap_95pct_interval_pct":
                [float(np.quantile(change,.025)),float(np.quantile(change,.975))],
            "bootstrap_draws":2500,
            "bootstrap_is_descriptive_over_original_selected_families":True
        }
    table_n=np.array([v["n_original_field_outcomes"] for v in fam_score.values()])
    return {
        "status":"EXPLORATORY_ORIGINAL_DATA_LEAVE_ONE_WHOLE_FAMILY_OR_GENUS_OUT",
        "heldout":"published_botanical_family" if group=="plant_family" else "published_botanical_genus",
        "original_source_paired_reproductive_outcomes":len(q),
        "original_published_studies":int(q.study_key.nunique()),
        "original_botanical_species":int(q.species.nunique()),
        "independent_heldout_groups":n_groups,
        "heldout_source_group_size_min":int(table_n.min()),
        "heldout_source_group_size_median":float(np.median(table_n)),
        "heldout_source_group_size_max":int(table_n.max()),
        "source_families_with_day_and_night_bloom_species":int(sum(
            v["both_bloom_classes"] for v in fam_score.values())),
        "fixed_features":list(BASE_COLS),
        "fixed_ridge_alpha":ALPHA,
        "all_models_same_exact_source_sample":True,
        "no_original_field_study_or_species_used_as_feature":True,
        "metrics":result,
        "individual_family_source_error_improvement_distribution":delta,
        "individual_species_and_family_identifiers_never_exported":True
    }


def purged_study_out(q:pd.DataFrame):
    """Blocked study-fold CV with additionally excluded TEST plant species
    from TRAIN; conventional leave-publication-out can share a species across
    different papers and overstate ecologically independent generalization."""
    require((len(q),q.study_key.nunique(),q.species.nunique())==
            (N_PAIRS,N_STUDIES,N_SPECIES),"Study/species CV source altered")
    from sklearn.model_selection import GroupKFold
    y=q.night_share.to_numpy(float)
    groups=q.study_key.to_numpy()
    w=q.source_study_weight.to_numpy(float)
    features={a:q[list(BASE_COLS)+list(extra)].to_numpy(float)
              for a,extra in CV_MODELS.items()}
    splits=list(GroupKFold(n_splits=5).split(q,y,groups))
    pred={name:np.full(len(q),np.nan) for name in CV_MODELS}
    detail=[]
    for k,(train,test) in enumerate(splits):
        held_species=set(q.iloc[test].species)
        purged=train[~q.iloc[train].species.isin(held_species).to_numpy()]
        info={
            "fold":k+1,
            "source_test_studies":int(q.iloc[test].study_key.nunique()),
            "source_test_pairs":len(test),
            "source_train_pairs_before_species_purge":len(train),
            "source_train_pairs_after_species_purge":len(purged),
            "original_train_studies_after_purge":int(q.iloc[purged].study_key.nunique()),
            "purged_source_training_pairs":int(len(train)-len(purged)),
        }
        detail.append(info)
        if (len(purged)<85 or info["original_train_studies_after_purge"]<38 or
            q.iloc[purged].night_bloom.nunique()!=2):
            return {
                "status":"NOT_ESTIMABLE_INDEPENDENT_STUDY_AND_TAXON_CV_TRAINING",
                "fold_support":detail
            }
        for name,XX in features.items():
            scaler=StandardScaler()
            tr=scaler.fit_transform(XX[purged]);ts=scaler.transform(XX[test])
            r=Ridge(alpha=ALPHA,fit_intercept=True)
            r.fit(tr,y[purged],sample_weight=w[purged]/w[purged].mean())
            pred[name][test]=np.clip(r.predict(ts),0,1)
    for name in pred:require(np.isfinite(pred[name]).all(),
                              "Missing source fold from independent study-species purge")
    scores={}
    wn=w/w.sum()
    for name in CV_MODELS:
        e=y-pred[name]
        scores[name]={
            "study_and_species_purged_RMSE_nightshare_fraction":float(
                np.sqrt(np.sum(wn*e**2))),
            "study_and_species_purged_MAE_nightshare_fraction":float(
                np.sum(wn*np.abs(e)))}
    b=scores[BASE_KEY]["study_and_species_purged_RMSE_nightshare_fraction"]
    for name in scores:
        scores[name]["RMSE_improvement_percent_vs_geography_only"]=float(
            100*(b-scores[name]["study_and_species_purged_RMSE_nightshare_fraction"])/b)
    return {
        "status":"EXPLORATORY_TRUE_STUDY_PLUS_SHARED_PLANT_PURGED_TRANSFER",
        "n_study_holdout_folds":len(splits),
        "n_source_pairs":len(q),
        "all_heldout_studies_and_heldout_species_removed_from_training":True,
        "all_folds_train_sufficient_independent_studies":True,
        "removed_source_train_pairs_range":[int(min(x["purged_source_training_pairs"] for x in detail)),
                                             int(max(x["purged_source_training_pairs"] for x in detail))],
        "scores":scores,
        "no_plant_taxon_or_paper_identifier_output":True
    }


def family_omission_coefficient(q:pd.DataFrame):
    """Test known nighttime bloom association under ALL whole-family
    deletions, without citing the 28 overlapping trials as replications."""
    estimates=[]
    not_estimable=0
    for group in q.plant_family.unique():
        a=q.loc[q.plant_family.ne(group)].copy()
        fit=source_fit(a,BLOOM_FORMULA,["night_bloom"])
        if (fit["status"]!="EXPLORATORY_DIRECT_REPRODUCTIVE_FUNCTION_TRAIT_MODEL"
            or "night_bloom" not in fit["terms"]):
            not_estimable+=1
            continue
        estimates.append({
            "beta":fit["terms"]["night_bloom"]["coefficient_fraction"],
            "ci":fit["terms"]["night_bloom"]["95CI_source_publication_cluster"]
        })
    if not estimates:return {"status":"NOT_ESTIMABLE_AFTER_FAMILY_DELETIONS"}
    x=np.array([z["beta"] for z in estimates])
    return {
       "status":"EXPLORATORY_PAIRED_REPRODUCTIVE_BLOOM_ASSOCIATION_FAMILY_LEVERAGE",
       "independent_botanical_family_deletions":int(q.plant_family.nunique()),
       "estimable_deletions":len(estimates),
       "not_estimable_deletions":not_estimable,
       "coefficient_fraction_range":[float(x.min()),float(x.max())],
       "median":float(np.median(x)),
       "fraction_positive":float((x>0).mean()),
       "n_nominal_study_cluster_CI_above_zero":int(sum(
            r["ci"][0]>0 for r in estimates)),
       "all_family_names_suppressed":True,
       "not_28_independent_experimental_replicates":True
    }


def run(joined:Path,receipt:Path,traits:Path,out:Path):
    if out.exists():raise FileExistsError("Refuse to overwrite audited source data")
    q,meta,support,baseline=source_data(joined,receipt,traits)
    result={
        "schema":"ALAN_original_DvN_flower_trait_true_reproductive_family_genus_transfer_20261113_v1",
        "status":"EXPLORATORY_ORIGINAL_301_FIELD_NIGHT_REPRO_LEAVE_FAMILY_OUT_AND_PURGED_STUDY",
        "source":{
            "original_301_paired_experiment_source_SHA256":meta["joined_file_sha256"],
            "signed_PCNL_years":[1992,2021],
            "original_301_comparisons":301,
            "eligible_same_completecase_morphology_pairs":len(q),
            "eligible_independent_publications":int(q.study_key.nunique()),
            "eligible_plant_species":int(q.species.nunique()),
            "eligible_botanical_families":int(q.plant_family.nunique()),
            "parent_source_bloom_repro_beta":baseline["terms"]["night_bloom"],
            "original_flower_trait_coverage":support
        },
        "whole_family_out_of_sample_predictions":lofo_models(q,"plant_family"),
        "whole_genus_out_of_sample_predictions":lofo_models(q,"genus"),
        "independent_original_publication_with_all_test_species_purged":
            purged_study_out(q),
        "whole_family_bloom_point_estimate_influence":family_omission_coefficient(q),
        "scientific_limits":{
            "not_local_effective_moth_pollination":True,
            "no_actual_ALAN_field_treatment_randomized":True,
            "year_prior_pcNL_remote_satellite_not_at_flower_lux":True,
            "source_true_day_night_reproduction_within_same_experiment":True,
            "source_bloom_morphology_from_literature_not_manipulated":True,
            "cross_family_generalization_not_a_test_of_biological_evolution":True,
            "original_field_author_species_and_study_selection_nonrandom":True,
            "CV_model_complexity_and_hyperparameter_frozen_before_fit":True,
            "family_cluster_bootstrap_descriptive_not_global_sampling_distribution":True,
            "all_studies_species_coordinates_privately_retained":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_ACTUAL_DAYNIGHT_FAMILY_TRANSFER_DIAGNOSTIC",
          json.dumps({
            "family":result["whole_family_out_of_sample_predictions"],
            "genus":result["whole_genus_out_of_sample_predictions"],
            "study_purged":result["independent_original_publication_with_all_test_species_purged"],
            "family_omission":result["whole_family_bloom_point_estimate_influence"]
          },ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("joined","receipt","traits","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.joined,a.receipt,a.traits,a.out)
