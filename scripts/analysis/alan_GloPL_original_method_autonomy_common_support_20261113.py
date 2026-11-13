#!/usr/bin/env python3
"""Audit identification of Sup vs Bagout x autonomous selfing x historical ALAN.

Original GloPL true field-experiment years 1993-2015, 1992..year-1 PCNL,
external source-direct botanical autonomous selfing, fixed source SHA
and old 1,419-experiment benchmark. Identify whether the apparent
method contrast can be reconstructed WITHIN comparable original studies,
species, years and location, instead of between disjoint research designs.

Scientific priority: Source support/positivity and regular-Sup-only outcome
association. Never interpret the old Sup minus Bagout β as causal method
response, local moth-pollination change or plant selfing compensation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from scripts.analysis.alan_glopl_direct_autonomous_selfing_PCNLPrior_20261113 import (
    original_source, source_traits, prepared_moderator,
    fit_source, METHOD_DESIGN_MAIN, METHOD_DIFF_TERM,
    PRIMARY, EXPOSURE, METH_COLS, require
)

SOURCE_COUNTS = (725, 242, 179)
CATEGORIES = ("Bagout", "Sup")
RESPONSE = "PL_log_ratio"
SELFS = (0,1)


def distinct_coverage(g:pd.DataFrame):
    """Comprehensive cross-level support, not 725 independent treatment trials."""
    q=g.copy()
    q["method"]=q.PL_Effect_Size_Type2.astype(str)
    require(set(q.method)==set(CATEGORIES) and set(q.autonomy)==set(SELFS),
            "Original method or directly verified selfing categories changed")
    method_blocks={}
    for m in CATEGORIES:
        c=q[q.method.eq(m)]
        method_blocks[m]={
            "experiment_cases":len(c),
            "species":int(c.species.nunique()),
            "publications":int(c.study_key.nunique()),
            "autonomous_cases":int(c.autonomy.sum()),
            "non_autonomous_cases":int((1-c.autonomy).sum()),
            "PL_absolute_gt2":int(c.PL_log_ratio.abs().gt(2).sum()),
            "zero_prior_PCNL_cases":int(c[EXPOSURE].eq(0).sum()),
        }
    count=lambda group,col: q.groupby(group,observed=True)[col].nunique()
    study_methods=count(["study_key"],"method")
    study_autonomy=count(["study_key"],"autonomy")
    plant_methods=count(["species"],"method")
    pair_methods=count(["study_key","species"],"method")
    original_year_methods=count(["study_key","species","exp_year"],"method")
    geo_methods=count(["study_key","species","exp_year","lat","lon"],"method")
    unique_study=q[["study_key","method","autonomy"]].drop_duplicates()
    research_pub_summary={
       "n_publications":int(len(study_methods)),
       "studies_only_one_method":int(study_methods.eq(1).sum()),
       "studies_with_BOTH_methods":int(study_methods.eq(2).sum()),
       "studies_with_BOTH_selfing_classes":int(study_autonomy.eq(2).sum()),
       "studies_with_BOTH_methods_AND_BOTH_selfing_classes":int(
           ((study_methods==2)&(study_autonomy==2)).sum()),
       "species_with_BOTH_methods":int(plant_methods.eq(2).sum()),
       "study_species_pairs_with_BOTH_methods":int(pair_methods.eq(2).sum()),
       "study_species_experimentyear_groups_with_BOTH_methods":int(original_year_methods.eq(2).sum()),
       "study_species_year_exact_fieldpoint_groups_with_BOTH_methods":int(geo_methods.eq(2).sum()),
       "publications_with_multiple_species_selfing_classes":int(study_autonomy.eq(2).sum()),
    }
    require((len(q),q.species.nunique(),research_pub_summary["n_publications"]) in
            ((725,242,179),(688,227,170)) and
            research_pub_summary["studies_with_BOTH_methods"]==3 and
            research_pub_summary["studies_with_BOTH_selfing_classes"]==10 and
            research_pub_summary["study_species_pairs_with_BOTH_methods"]==1 and
            research_pub_summary["study_species_experimentyear_groups_with_BOTH_methods"]==0,
            "Source actual study-level method overlap changed")
    # An apparent pooled three-way interaction cannot identify treatment effects
    # from only one study×plant pair, none within the same real experiment year.
    research_pub_summary["within_same_published_experiment_method_comparison"] = (
       "NOT_IDENTIFIABLE_FROM_SOURCE_MATCHED_DESIGN")
    return {"source_method_groups":method_blocks,"research_unit_overlap":research_pub_summary}


def spatial_joint_positivity(g:pd.DataFrame, degrees=(5,10,20,30)):
    """Four source-positive combinations are needed per region to compare
    method × selfing, PLUS overlapping measured prior light exposure.
    A geographic FE alone cannot manufacture nonexistent within-region cells.
    """
    res={}
    for degree in degrees:
        q=g.copy()
        q["region"]=((q.lat+90)//degree).astype(int).astype(str)+"_"+(
                     ((q.lon+180)//degree).astype(int).astype(str))
        groups=q.groupby("region",observed=True)
        n4=[]
        max_overlap=[]
        for _,d in groups:
            subgroup=d.groupby(["PL_Effect_Size_Type2","autonomy"],observed=True)
            if subgroup.ngroups!=4:
                continue
            n4.append(d)
            ends=subgroup[EXPOSURE].agg(["min","max"])
            max_overlap.append(float(ends["max"].min()-ends["min"].max()))
        both_method=int(groups.PL_Effect_Size_Type2.nunique().eq(2).sum())
        both_self=int(groups.autonomy.nunique().eq(2).sum())
        total4=int(len(n4))
        if n4:
            row_count=int(sum(len(d) for d in n4))
            papers_count=int(pd.concat(n4,ignore_index=True).study_key.nunique())
        else:
            row_count=0
            papers_count=0
        res[str(degree)+"degree"]={
           "all_regions":int(q.region.nunique()),
           "regions_with_both_methods":both_method,
           "regions_with_both_selfing_states":both_self,
           "regions_with_all_four_method_x_selfing_source_cells":total4,
           "experiment_cases_in_four_source_category_regions":row_count,
           "publications_in_four_source_category_regions":papers_count,
           "regions_with_all_four_AND_positive_exposure_interval_overlap":int(
               sum(v>.01 for v in max_overlap)),
           "source_coverage_not_independent_field_replicates":True,
           "degree_width":degree,
        }
    require(res["10degree"]["all_regions"]==(75 if len(g)==725 else 72) and
            res["10degree"]["regions_with_all_four_method_x_selfing_source_cells"]==1 and
            res["10degree"]["regions_with_all_four_AND_positive_exposure_interval_overlap"]==0,
            "Original source 10deg method-selfing-priorlight identification support changed")
    return res


def within_research_light_support(g:pd.DataFrame, minimum_study_light_range=.03):
    """For autonomous-selfing contrast a single study with both selfing types
    and real light variation is more interpretable than unrelated papers.
    """
    s=g.groupby("study_key",observed=True).agg(
       method_classes=("PL_Effect_Size_Type2","nunique"),
       self_classes=("autonomy","nunique"),
       light_range=(EXPOSURE,lambda x:float(x.max()-x.min())),
       cases=("original_row_id","size"))
    both=s.self_classes.eq(2)
    genuine=both & s.light_range.ge(minimum_study_light_range)
    n_both=(s.method_classes.eq(2)&s.self_classes.eq(2))
    return {
        "source_publications_with_both_selfing_types":int(both.sum()),
        "source_publications_both_selfing_types_AND_light_range_at_least_003":int(genuine.sum()),
        "source_experiments_in_both_selfing_and_light_varying_publications":int(g.loc[
             g.study_key.isin(s.index[genuine])].shape[0]),
        "publications_with_both_methods_both_selfing":int(n_both.sum()),
        "publications_with_both_methods_both_selfing_and_light_range_at_least_003":int(
             (n_both&s.light_range.ge(minimum_study_light_range)).sum()),
        "independent_study_support_cannot_identify_method_specific_ALAN_causality":True,
    }


def sup_only_selfing_direct(g:pd.DataFrame,trim=False):
    """Primary valid hand-pollen SUPPLEMENTATION protocol, direct actual
    pollen-limitation outcome. No ecological interaction sign search.

    A selfing coefficient is a cross-plant observational comparison of actual
    measured pollen-limitation output, not an intervention on selfing.
    """
    a=g.loc[g.PL_Effect_Size_Type2.eq("Sup")].copy()
    if trim:
        a=a.loc[a.PL_log_ratio.abs().le(2)].copy()
    if (len(a)<400 or a.study_key.nunique()<100 or
         a.loc[a.autonomy.eq(1),"species"].nunique()<20):
        return {"status":"NOT_ESTIMABLE_SUP_STUDY_SOURCE_SUPPORT",
                "cases":len(a),"studies":int(a.study_key.nunique())}
    a["w"]=1.0/a.groupby("study_key",observed=True).study_key.transform("size")
    typ=" + ".join("C("+k+")" for k in METH_COLS)
    formula=(RESPONSE+" ~ autonomy + "+EXPOSURE+
      " + exp_year + lat + lat2 + lon + lon2 + "+typ+" + C(geo10)")
    m=smf.wls(formula,data=a,weights=a.w).fit(
       cov_type="cluster",cov_kwds={"groups":a.study_key},use_t=False)
    rank=int(np.linalg.matrix_rank(m.model.exog))
    require(rank==len(m.params),"Sup-only high-dim source fixed effects rank deficient")
    beta=float(m.params["autonomy"])
    ci=m.conf_int().loc["autonomy"]
    return {
       "status":"EXPLORATORY_ORIGINAL_SUP_EXPERIMENT_POLLEN_LIMITATION_BY_DIRECT_SELFING",
       "beta_PL_log_reproductive_supplementation_ratio_for_autonomous_vs_absent":beta,
       "CI95_original_publication_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
       "p_nominal":float(m.pvalues["autonomy"]),
       "n_experiments":int(len(a)),
       "n_publications":int(a.study_key.nunique()),
       "n_plant_species":int(a.species.nunique()),
       "autonomous_botanical_species":int(a.loc[a.autonomy.eq(1),"species"].nunique()),
       "nonautonomous_botanical_species":int(a.loc[a.autonomy.eq(0),"species"].nunique()),
       "original_publication_equal_weight":True,
       "formula":formula,
       "prior_light_is_not_experimental_ALAN":True,
    }


def run(joined:Path,receipt:Path,ledger:Path,out:Path):
    if out.exists():
        raise FileExistsError("Refusing to overwrite existing original source science")
    original,receipts,benchmark=original_source(joined,receipt)
    results={}
    for tier in ("high_medium","high_only"):
        d=prepared_moderator(original,source_traits(ledger,high_only=tier=="high_only"))
        if tier=="high_medium":
            require((len(d),d.species.nunique(),d.study_key.nunique())==SOURCE_COUNTS,
                    "Original 725 source overlap changed")
        source=distinct_coverage(d)
        spatial=spatial_joint_positivity(d)
        study=within_research_light_support(d)
        parent=fit_source(
           d,METHOD_DESIGN_MAIN+" + C(geo10)",term=METHOD_DIFF_TERM)
        sup=sup_only_selfing_direct(d)
        sup_trim=sup_only_selfing_direct(d,trim=True)
        results[tier]={
          "source_study_and_protocol_overlap":source,
          "joint_spatial_and_light_exposure_common_support":spatial,
          "same_publication_trait_and_illuminance_source_support":study,
          "repeated_original_DESCRIPTIVE_method_heterogeneity_coefficient":parent,
          "actual_nonbagged_supplementation_source_selfing_association":sup,
          "same_sup_response_after_abs_PL_gt2_exclusion":sup_trim,
          "adjudication":{
            "Sup_Bagout_historicallight_selfing_difference_causal_effect_identified":False,
            "method_specific_floristic_autogamy_response_identifiable_within_same_real_experiment":False,
            "geographic_FE_cannot_restore_unmeasured_method_selfing_exposure_quartets":True,
            "method_aligned_source_regions_are_not_random_method_assignments":True,
            "normal_Sup_only_biological_outcome_measured":True,
            "study_population_selfing_measured_locally":False,
            "actual_moth_effective_pollination_measured":False
          }
        }
    main=results["high_medium"]
    check=main["repeated_original_DESCRIPTIVE_method_heterogeneity_coefficient"]
    require(check["status"]=="EXPLORATORY_ORIGINAL_EXPERIMENT_PL_PRIOR_PCNL_TRAIT_MODERATION" and
            abs(check["beta"]-.86994)<.00005,
            "Historical original descriptive method coefficient changed")
    stats=main["source_study_and_protocol_overlap"]["research_unit_overlap"]
    require(stats["studies_with_BOTH_methods"]<10 and
            stats["study_species_experimentyear_groups_with_BOTH_methods"]==0,
            "New source replication unexpected; review hard stop")
    data={
        "schema":"ALAN_GloPL_725_selfing_priorlight_treatment_method_joint_positivity_audit_20261113_v1",
        "status":"AUTHENTIC_SOURCE_METHOD_HETEROGENEITY_NOT_WITHIN_STUDY_IDENTIFIED",
        "source":{"GloPL_true_experiment_years":[1993,2015],
                  "original_PCNL_annual_years":[1992,2014],
                  "original_1419_reproduced_light_beta":float(benchmark.params[EXPOSURE]),
                  "source_full_GloPL_original_SHA256":receipts["output_sha256"],
                  "source_DIRECT_selfing_ledger_SHA256":"7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
                  "published_GloPL_outcome":"hand-supplementation advantage log reproductive success ratio"},
        "results":results,
        "conclusion":"Original Sup−Bagout difference is a predominantly BETWEEN-STUDY/BETWEEN-REGION descriptive result; no source-supported within-study same-species same-year comparison or four-class region light common support to infer Bagout vs Sup biology.",
        "limitations":{
          "source_not_native_wild_island_flora_only":True,
          "external_species_level_autonomy_not_measured_in_source_field_population":True,
          "rolling_prior_PCNL_not_local_manipulated_light":True,
          "treatment_design_is_nonrandom_across_articles_and_geography":True,
          "not_moth_specific_pollination":True,
          "not_causal_method_or_selfing_compensation":True,
          "all_treatment_effect_and_overlap_results_posthoc":True,
          "no_species_or_study_identifiers_exported":True,
          "no_location_by_trait_or_restricted_source_rows_exported":True,
        }}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("VERIFIED_GLOPL_725_JOINT_DESIGN_SUPPORT",
          json.dumps({"source":main["source_study_and_protocol_overlap"],
                      "spatial":main["joint_spatial_and_light_exposure_common_support"],
                      "within_study":main["same_publication_trait_and_illuminance_source_support"],
                      "Sup":main["actual_nonbagged_supplementation_source_selfing_association"],
                      "SupTrim":main["same_sup_response_after_abs_PL_gt2_exclusion"]},
                     ensure_ascii=False),flush=True)
    return data


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for item in ("joined","receipt","ledger","out"):
        ap.add_argument("--"+item,required=True,type=Path)
    a=ap.parse_args()
    run(a.joined,a.receipt,a.ledger,a.out)
