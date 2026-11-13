#!/usr/bin/env python3
"""Genuine published pollen-supplementation results × preceding-year nightlight
  moderated by independently curated species-DIRECT autonomous selfing.

Source population: source-SHA-locked 1,419 original GloPL experimental cases
with TRUE experiment years 1993-2015, true original experimental site coordinates
and actual 1992 to experiment-year-minus-one Li-Cao corrected PCNL satellite.
Independent botanical trait source: fixed Island 2026-09 source-scale direct
species high/medium and high-only ledger. Unknown does NOT mean no selfing.

Primary contrast is the original published GloPL model's geography/method/year
control + prior light × autonomous-selfing. Preexisting original 10deg region FE
is a sensitivity, not chosen posthoc for significance. Full study+species
two-way cluster, quality, outlier, supplementation type, lit-only,
genus-rank and exhaustive leave-one-study/plant sensitivity are all reported.

Biological inferential ceiling: source-ecological correlation, NOT light-caused
pollen limitation, confirmed local moth reliance, or local rescue by autogamy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
import statsmodels.formula.api as smf
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups

from scripts.analysis.alan_glopl_real_PL_priorPCNL_exploratory_20261110 import (
    validate, gate_source, METH_COLS, EXPOSURE, require
)
from scripts.analysis.alan_real_daynight_reproductive_PCNL_562_20261113 import (
    experimental_daynight
)

TRAIT_LEDGER_SHA="7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577"
N_LEDGER_ROWS=134268
N_LEDGER_SPECIES=62423
BENCHMARK=0.010715792632233193
TARGET="pcnl_prior_mean_log1p:autonomy"
REPRO_CLASSES=("absent","autonomous")
METHOD_FORMULA=" + ".join("C("+k+")" for k in METH_COLS)
CONTROLS="exp_year + lat + lat2 + lon + lon2 + "+METHOD_FORMULA
METHOD_DIFF_TERM=EXPOSURE+":autonomy:C(PL_Effect_Size_Type2)[T.Sup]"
METHOD_DESIGN_CONTROLS=(
    "exp_year + lat + lat2 + lon + lon2 + "+
    " + ".join("C("+k+")" for k in METH_COLS if k!="PL_Effect_Size_Type2"))
METHOD_DESIGN_MAIN="PL_log_ratio ~ "+EXPOSURE+"*autonomy*C(PL_Effect_Size_Type2) + "+METHOD_DESIGN_CONTROLS

PRIMARY="PL_log_ratio ~ "+EXPOSURE+"*autonomy + "+CONTROLS
BASE="PL_log_ratio ~ "+EXPOSURE+" + "+CONTROLS


def original_source(joined, receipt):
    g, rec=validate(joined,receipt)
    g=gate_source(g)
    require(len(g)==1419 and g.original_row_id.is_unique,
            "Original 1419 site/experiment-light source changed")
    require(g.pcnl_historical_complete.astype(str).str.lower().eq("true").all(),
            "False original historic PCNL complete flags")
    require(g.pcnl_pre_experiment_expected_years.eq(
        g.pcnl_pre_experiment_observed_years).all(),
        "Incomplete independent real annual historical satellite exposure")
    require(g.year_exp.between(1993,2015).all(),"Experiment year out of source range")
    g["PL_log_ratio"]=pd.to_numeric(g.PL_Effect_Size,errors="coerce")
    require(g.PL_log_ratio.notna().all(),"Missing original pollen-addition outcome")
    g["exp_year"]=g.year_exp.astype(int)
    g["lat2"]=g.lat.astype(float)**2
    g["lon2"]=g.lon.astype(float)**2
    for col in METH_COLS:
        g[col]=g[col].astype(str).fillna("unknown").replace(
            {"nan":"unknown","":"unknown"})
    g["geo10"]=((g.lat+90)//10).astype(int).astype(str).str.cat(
         ((g.lon+180)//10).astype(int).astype(str),sep="_")
    g["geo20"]=((g.lat+90)//20).astype(int).astype(str).str.cat(
         ((g.lon+180)//20).astype(int).astype(str),sep="_")
    g["genus"]=g.species.str.split().str[0]
    g["study_equal_weight"]=1/g.groupby("study_key").study_key.transform("size")
    require(g.study_key.nunique()==341 and g.species.nunique()==507,
            "True GloPL source publication/species support changed")
    benchmark=smf.wls(BASE,data=g,weights=g.study_equal_weight).fit(
        cov_type="cluster",cov_kwds={"groups":g.study_key},use_t=False)
    require(abs(float(benchmark.params[EXPOSURE])-BENCHMARK)<1e-8,
            "Prior original 1419 experimental pollen-limitation coefficient not reproduced")
    return g,rec,benchmark


def source_traits(path, high_only=False):
    require(hashlib.sha256(path.read_bytes()).hexdigest()==TRAIT_LEDGER_SHA,
            "Island source species-direct selfing trait ledger changed")
    d=pd.read_csv(path,dtype=str,keep_default_na=False)
    require((len(d),d.accepted_species.nunique())==
            (N_LEDGER_ROWS,N_LEDGER_SPECIES),
            "Island exact direct-species source database dimensions changed")
    q=d.loc[d.trait_name.eq("autonomous_selfing_capacity") &
        d.resolution_status.eq("resolved") &
        d.quality.isin(["high"] if high_only else ["high","medium"]) &
        d.multistate.eq("False") &
        d.normalized_value.isin(REPRO_CLASSES)].copy()
    require(q.accepted_species.is_unique,"Conflicting independent direct autonomous selfing traits")
    return q[["accepted_species","normalized_value"]]


def prepared_moderator(g, trait):
    x=g.merge(trait,left_on="species",right_on="accepted_species",
              how="inner",validate="many_to_one")
    x["autonomy"]=x.normalized_value.eq("autonomous").astype(int)
    x["study_equal_weight"]=1/x.groupby("study_key").study_key.transform("size")
    require(x.original_row_id.is_unique,
            "Duplicate original pollen supplementation case from direct trait join")
    return x


def fit_source(g, formula, term=TARGET,*,twoway=True):
    """No inferred p-values for rank-deficient fixed effects."""
    if (len(g)<100 or g.study_key.nunique()<30 or
        min(g.groupby("autonomy").species.nunique().to_dict().values(),default=0)<10):
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_STUDY_TRAIT_SUPPORT",
                "cases":int(len(g)),"studies":int(g.study_key.nunique()),
                "species":int(g.species.nunique())}
    gg=g.copy()
    gg["study_equal_weight"]=1/gg.groupby("study_key").study_key.transform("size")
    m=smf.wls(formula,data=gg,weights=gg.study_equal_weight).fit(
        cov_type="cluster",cov_kwds={"groups":gg.study_key},use_t=False)
    rank=int(np.linalg.matrix_rank(m.model.exog))
    if rank<len(m.params) or term not in m.params:
        return {"status":"NOT_ESTIMABLE_DESIGN_COLLINEARITY",
                "rank":rank,"parameters":int(len(m.params)),
                "cases":len(gg),"studies":int(gg.study_key.nunique())}
    ci=m.conf_int().loc[term]
    res={
        "status":"EXPLORATORY_ORIGINAL_EXPERIMENT_PL_PRIOR_PCNL_TRAIT_MODERATION",
        "beta":float(m.params[term]),
        "CI95_study_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
        "p_nominal_study_cluster":float(m.pvalues[term]),
        "SE_study_cluster":float(m.bse[term]),
        "cases":int(len(gg)),
        "studies":int(gg.study_key.nunique()),
        "species":int(gg.species.nunique()),
        "geo10_blocks":int(gg.geo10.nunique()),
        "geographic_blocks_with_both_trait_states":int(
            gg.groupby("geo10").autonomy.nunique().eq(2).sum()),
        "parameters":len(m.params),
        "design_rank":rank,
        "response":"published actual hand-supplementation log reproductive-output ratio",
        "covariates_formula":formula,
        "coefficient":term
    }
    if twoway:
        a=pd.factorize(gg.study_key,sort=False)[0]
        b=pd.factorize(gg.species,sort=False)[0]
        var,_,_=cov_cluster_2groups(m,a,b)
        i=list(m.params.index).index(term)
        variance=float(var[i,i])
        if variance>0 and np.isfinite(variance):
            se=float(np.sqrt(variance))
            dof=min(gg.study_key.nunique()-1,gg.species.nunique()-1)
            q=float(student_t.ppf(.975,dof))
            pv=float(2*student_t.sf(abs(float(m.params[term])/se),dof))
            res["two_way_study_and_species_cluster"]={
                "SE":se,
                "CI95_student_t":[float(m.params[term]-q*se),
                                   float(m.params[term]+q*se)],
                "p_two_sided_t":pv,
                "df":int(dof)}
        else:
            res["two_way_study_and_species_cluster"]={
                "status":"COVARIANCE_NOT_POSITIVE_SEMIDEFINITE"}
    return res


def general_source_variation(g):
    cells=g.groupby("geo10",observed=True).agg(
        cases=("autonomy","size"),trait_states=("autonomy","nunique"),
        sd_light=(EXPOSURE,"std"))
    both=cells.trait_states.eq(2)
    both_lit=both & cells.sd_light.gt(.03)
    return {
        "original_experiment_cases":int(len(g)),
        "plant_species":int(g.species.nunique()),
        "autonomous_species":int(g.loc[g.autonomy.eq(1),"species"].nunique()),
        "no_autonomous_species":int(g.loc[g.autonomy.eq(0),"species"].nunique()),
        "source_experiments_autonomous":int(g.autonomy.sum()),
        "source_experiments_non_autonomous":int((1-g.autonomy).sum()),
        "independent_original_publications":int(g.study_key.nunique()),
        "publications_with_both_autonomy_classes":int(
            g.groupby("study_key").autonomy.nunique().eq(2).sum()),
        "published_source_plant_taxa_with_multiple_publications":int(
            g.groupby("species").study_key.nunique().gt(1).sum()),
        "original_geography_10deg_blocks":int(len(cells)),
        "blocks_with_both_autonomous_states":int(both.sum()),
        "blocks_with_both_autonomy_and_nightlight_variation_sd_over_003":int(both_lit.sum()),
        "experiments_in_those_blocks":int(g.geo10.isin(set(cells.index[both_lit])).sum()),
        "studies_in_those_blocks":int(g.loc[g.geo10.isin(set(cells.index[both_lit])),"study_key"].nunique()),
        "fraction_zero_light":float(g[EXPOSURE].eq(0).mean()),
        "fraction_outcome_abs_gt2":float(g.PL_log_ratio.abs().gt(2).mean()),
        "no_plant_specific_moth_study_in_this_PL_experiment":True
    }


def complete_leaveouts(g,formula,baseline,max_cases=1000):
    """All one-at-a-time plant and study omissions of same *original* 10deg FE.
    This is a robustness probe of ONE posthoc estimator, not 100s of new tests.
    """
    estimates={}
    for label,field in (("study","study_key"),("plant_species","species")):
        vals=[]
        rank_fail=0
        for item in g[field].unique():
            sub=g[g[field].ne(item)]
            if (sub.autonomy.nunique()!=2 or sub.study_key.nunique()<30):
                rank_fail+=1
                continue
            sub=sub.copy()
            sub["study_equal_weight"]=1/sub.groupby("study_key").study_key.transform("size")
            fit=smf.wls(formula,data=sub,weights=sub.study_equal_weight).fit()
            if TARGET not in fit.params or np.linalg.matrix_rank(fit.model.exog)!=len(fit.params):
                rank_fail+=1
                continue
            vals.append(float(fit.params[TARGET]))
        if vals:
            arr=np.array(vals)
            estimates[label]={
                "n_original_deletions":int(g[field].nunique()),
                "n_rank_estimable":len(vals),
                "n_inestimable":rank_fail,
                "beta_range":[float(arr.min()),float(arr.max())],
                "median":float(np.median(arr)),
                "fraction_negative":float((arr<0).mean()),
                "n_nonnegative":int((arr>=0).sum()),
                "max_abs_delta_from_full":float(np.max(np.abs(arr-baseline))),
                "taxon_and_publication_identities_suppressed":True
            }
        else:
            estimates[label]={"status":"NOT_ESTIMABLE_AFTER_DELETIONS","n_deletions":int(g[field].nunique())}
    return estimates


def daynight_overlay(g,source_study,source_trait):
    """Source support only unless each independent night/day and selfing class
    has 10+ accepted plants and 20 independent publications on original GloPL.
    Not every nighttime reproductive experiment is moth pollination.
    """
    experiments,meta=experimental_daynight(source_study,source_trait)
    match=g.merge(experiments[["species","night_stronger","day_stronger",
                               "night_share","n_studies"]],on="species",
                  how="inner",validate="many_to_one")
    support={
      "all_author_experimental_reproductive_species":len(experiments),
      "GloPL_original_1419_source_species_matched_experiment_reproductive_class":int(
          g[g.species.isin(experiments.species)].species.nunique()),
      "autonomous_selfing_strict_source_species_with_valid_diel_experiment":int(match.species.nunique()),
      "GloPL_outcome_cases_with_both_direct_autonomy_and_diel_experiment":len(match),
      "independent_GloPL_publications":int(match.study_key.nunique()),
      "autonomy_x_night_study_species_counts":{
          str(k):int(v) for k,v in match.groupby(
              ["autonomy","night_stronger"]).species.nunique().items()},
      "biological_identification":"Daytime/nighttime exclusion can include non-moth pollinators; no local ALAN-at-flower treatment"
    }
    if (len(match)<200 or match.study_key.nunique()<60 or
        (match.groupby(["autonomy","night_stronger"]).species.nunique().min()
         if len(match) else 0)<10):
        support["ALAN_x_autonomy_x_night_study"]="NOT_ESTIMABLE_FUNCTIONAL_COMPARISON_SUPPORT"
    else:
        support["ALAN_x_autonomy_x_night_study"]="SUPPORT_EXISTS_NOT_FITTED_UNREGISTERED_3WAY"
    return support


def run(joined,receipt,ledger,author_study,author_trait,out):
    out=Path(out)
    if out.exists():
        raise FileExistsError("Refuse overwriting original science output")
    base,source,benchmark=original_source(joined,receipt)
    results={}
    cohorts={}
    main=None
    for tier in ("high_medium","high_only"):
        g=prepared_moderator(base,source_traits(ledger,high_only=tier=="high_only"))
        require(len(g)>500 and g.study_key.nunique()>125,
                "Original global original-experiment/selfing overlap unexpectedly sparse")
        if tier=="high_medium":
            require((len(g),g.species.nunique(),g.study_key.nunique())==(725,242,179),
                    "Original high/medium source data changed; cannot reinterpret sign")
        cohorts[tier]=general_source_variation(g)
        formulas={
            "original_geography_year_and_method_controls":PRIMARY,
            "preexisting_10deg_spatial_fixed_effects":PRIMARY+" + C(geo10)",
            "alternative_20deg_spatial_fixed_effects":PRIMARY+" + C(geo20)",
            "10deg_remove_abs_PL_logratio_gt2":PRIMARY+" + C(geo10)",
            "10deg_only_actual_positive_prior_PCNL":PRIMARY+" + C(geo10)",
            "10deg_supplementation_only_SourceType2_Sup":PRIMARY+" + C(geo10)",
            "10deg_add_botanical_genus_fixed_effects":PRIMARY+" + C(geo10) + C(genus)"
        }
        subsets={
            "10deg_remove_abs_PL_logratio_gt2":g.loc[g.PL_log_ratio.abs().le(2)],
            "10deg_only_actual_positive_prior_PCNL":g.loc[g[EXPOSURE].gt(0)],
            "10deg_supplementation_only_SourceType2_Sup":g.loc[
                g.PL_Effect_Size_Type2.eq("Sup")]
        }
        fits={}
        for name,formula in formulas.items():
            gsource=subsets.get(name,g)
            print("ALAN_GLOPL_DIRECT_TRAIT_MODEL",tier,name,len(gsource),flush=True)
            fits[name]=fit_source(gsource,formula)
        method_contrasts={
            "GloPL_Bagout_vs_Sup_10deg_source_design": (
                 METHOD_DESIGN_MAIN+" + C(geo10)",g),
            "GloPL_Bagout_vs_Sup_20deg_source_design": (
                 METHOD_DESIGN_MAIN+" + C(geo20)",g),
            "GloPL_Bagout_vs_Sup_10deg_no_absPL_gt2": (
                 METHOD_DESIGN_MAIN+" + C(geo10)",
                 g.loc[g.PL_log_ratio.abs().le(2)]),
        }
        method_result={}
        for model_name,(formula,subset) in method_contrasts.items():
            method_result[model_name]=fit_source(
                subset,formula,term=METHOD_DIFF_TERM)
        summaries=[]
        for src_type,v in g.groupby("PL_Effect_Size_Type2",observed=True):
            summaries.append({
                "original_GloPL_method_type":src_type,
                "original_genuine_PL_experiment_rows":int(len(v)),
                "source_plant_species":int(v.species.nunique()),
                "source_publications":int(v.study_key.nunique()),
                "source_autonomous_cases":int(v.autonomy.sum()),
                "source_no_autonomous_cases":int((1-v.autonomy).sum()),
                "source_PL_abs_gt2_cases":int(v.PL_log_ratio.abs().gt(2).sum()),
                "source_dark_PCNL_zero_case_count":int(v[EXPOSURE].eq(0).sum()),
            })
        fits["source_pollination_treatment_design_method_contrast"]={
            "posthoc_diagnostic_NO_CLAIM_OF_PHYSIOLOGICAL_MECHANISM":True,
            "GloPL_Bagout_differs_from_Sup_in_excluding_natural_pollination":True,
            "source_method_groups":summaries,
            "Sup_minus_Bagout_PCNL_autonomy_moderation":method_result,
        }
        results[tier]=fits
        if tier=="high_medium":
            main=g
            results[tier]["nighttime_functional_external_source_support"]=daynight_overlay(
                g,author_study,author_trait)
            if (fits["preexisting_10deg_spatial_fixed_effects"]["status"]==
                "EXPLORATORY_ORIGINAL_EXPERIMENT_PL_PRIOR_PCNL_TRAIT_MODERATION"):
                results[tier]["all_study_and_plant_leaveouts_10deg_source"]=complete_leaveouts(
                   g,PRIMARY+" + C(geo10)",
                   fits["preexisting_10deg_spatial_fixed_effects"]["beta"])
    result={
        "schema":"ALAN_GLOPL_1419_orig_experiment_datelocked_priorPCNL_direct_autonomous_selfing_moderation_v1",
        "status":"VERIFIED_ORIGINAL_EXPERIMENT_OUTCOME_PRIOR_PCNL_TRAIT_EXPLORATORY",
        "source":{
            "true_experiment_year_range":[1993,2015],
            "real_preexperiment_PCNL_years":[1992,2014],
            "original_GloPL_experiments":1419,
            "original_GloPL_publications":341,
            "original_benchmark_beta":float(benchmark.params[EXPOSURE]),
            "original_benchmark_CI95":benchmark.conf_int().loc[EXPOSURE].tolist(),
            "original_GloPL_plus_PCNL_sha256":source["output_sha256"],
            "original_island_direct_species_trait_sha256":TRAIT_LEDGER_SHA,
            "original_author_exclusion_study_blob_sha":"39a82d966afcdde670603add1c881ace8791ffe9",
            "excluded_multistate_and_unmeasured_selfing_as_unknown":True,
        },
        "cohort_support":cohorts,
        "all_models":results,
        "limitations":{
            "not_effective_moth_pollination":True,
            "selfing_direct_trait_external_not_measured_at_GloPL_study_population":True,
            "GloPL_log_pollen_supplementation_advantage_not_local_moth_abundance":True,
            "previous_PCNL_sensing_not_ground_illumination":True,
            "nonrandom_external_autonomous_selfing_trait_coverage":True,
            "geography_control_changes_interaction_conclusion":True,
            "preplanned_in_an_existing_source_exploratory_after_many_analyses_NOT_CONFIRMATORY":True,
            "no_PL_outlier_exclusion_without_reporting_full":True,
            "unobserved_guild_not_biologically_absent":True,
            "independent_local_plant_and_taxon_experimental_source_heterogeneity":True,
            "no_original_restricted_species_trait_or_experimental_coordinates_redistributed":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    summary={k:{
         "support":cohorts[k],
         "models":{z:{"status":v.get("status"),"beta":v.get("beta"),
                      "CI":v.get("CI95_study_cluster"),
                      "p":v.get("p_nominal_study_cluster"),
                      "two_way":v.get("two_way_study_and_species_cluster")}
                   for z,v in result["all_models"][k].items()
                   if isinstance(v,dict) and z!="all_study_and_plant_leaveouts_10deg_source"},
         "night":result["all_models"][k].get("nighttime_functional_external_source_support"),
         "LOO":result["all_models"][k].get("all_study_and_plant_leaveouts_10deg_source")}
             for k in cohorts}
    print("ALAN_VERIFIED_TRUE_GLOPL_SELFING_PRIORLIGHT_MODELS",
          json.dumps(summary,ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for x in ("joined","receipt","ledger","author_study","author_trait","out"):
        parser.add_argument("--"+x,type=Path,required=True)
    a=parser.parse_args()
    run(a.joined,a.receipt,a.ledger,a.author_study,a.author_trait,a.out)
