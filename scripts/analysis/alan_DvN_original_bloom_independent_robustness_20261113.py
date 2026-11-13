#!/usr/bin/env python3
"""Independent author botanical bloom time × TRUE paired night-vs-day
FRUIT/SEED reproductive function in the MD5 verified original 301 field
site/season/year sample; evaluate genus/family, geography, exposure and
publication/species dependence BEFORE ecological interpretation.

The day/night treatments are paired within each published study plant
site/year; comparisons of different plant blooming habits remain mostly
between studies. Not a causal ALAN mechanism, and not moth-specific service.
"""
from __future__ import annotations

import argparse
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
from scripts.analysis.alan_DvN_real_field_light_study_influence_bloom_20261113 import (
    original_trait_bloom
)

TARGET="night_bloom"
BENCHMARK=.21474990651679526
SOURCE_MODEL="night_share ~ "+LIGHT+" + night_bloom + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)"
SOURCE_INTERACTION="night_share ~ "+LIGHT+"*night_bloom + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)"


def valid_morphological_groups(x):
    s=x.copy()
    s=s.loc[s.bloom_period_simple.isin(("night","day"))].copy()
    s[TARGET]=s.bloom_period_simple.eq("night").astype(int)
    s["source_study_weight"]=1/s.groupby("study_key").study_key.transform("size")
    s["geo20"]=((s.lat+90)//20).astype(int).astype(str)+"_"+(
               ((s.lon+180)//20).astype(int).astype(str))
    s["family"]=s.plant_family.fillna("unknown").replace(
        {"NA":"unknown","":"unknown"}).astype(str)
    both=s.groupby("study_key")[TARGET].nunique().eq(2)
    mixed_family=s.groupby("family")[TARGET].nunique().eq(2)
    return s,{
       "source_pairs":len(s),
       "independent_studies":int(s.study_key.nunique()),
       "species":int(s.species.nunique()),
       "night_bloom_pairs":int(s.night_bloom.sum()),
       "day_bloom_pairs":int((1-s.night_bloom).sum()),
       "night_bloom_studies":int(s.loc[s.night_bloom.eq(1),"study_key"].nunique()),
       "day_bloom_studies":int(s.loc[s.night_bloom.eq(0),"study_key"].nunique()),
       "source_studies_with_both_bloom_classes":int(both.sum()),
       "source_study_pairs_in_both_class_papers":int(
           s.study_key.isin(set(both.index[both])).sum()),
       "source_families":int(s.family.nunique()),
       "families_with_both_day_and_night_blooming_species":int(mixed_family.sum()),
       "source_bloom_traits_from_original_independent_author_species_table":True,
    }


def fitting(d,formula,term=TARGET,twoway=True):
    if (len(d)<75 or d.study_key.nunique()<35 or
         min(d.loc[d.night_bloom.eq(k),"study_key"].nunique() for k in (0,1))<15):
        return {"status":"NOT_ESTIMABLE_BOTANICAL_BLOOM_STUDY_SUPPORT",
                "pairs":len(d),"studies":int(d.study_key.nunique())}
    x=d.copy()
    x["source_study_weight"]=1/x.groupby("study_key").study_key.transform("size")
    fit=smf.wls(formula,data=x,weights=x.source_study_weight).fit(
        cov_type="cluster",cov_kwds={"groups":x.study_key},use_t=False)
    mat=fit.model.exog
    rank=int(np.linalg.matrix_rank(mat))
    if rank!=len(fit.params) or term not in fit.params or not np.isfinite(fit.bse[term]):
        return {"status":"NOT_ESTIMABLE_FIXED_EFFECT_RANK",
                "pairs":len(x),"n_covariates":len(fit.params),"design_rank":rank}
    ci=fit.conf_int().loc[term]
    out={"status":"EXPLORATORY_ACTUAL_FIELD_REPRODUCTION_BLOOM_TIME_ASSOCIATION",
         "beta_night_reproductive_share_fraction":float(fit.params[term]),
         "CI95_publication_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
         "p_nominal":float(fit.pvalues[term]),
         "term":term,"source_model":formula,
         "source_original_study_pairs":len(x),
         "source_original_publications":int(x.study_key.nunique()),
         "source_original_botanical_species":int(x.species.nunique()),
         "n_parameters":len(fit.params)}
    if twoway:
        aa=pd.factorize(x.study_key)[0]
        bb=pd.factorize(x.species)[0]
        try:
            cov,_,_=cov_cluster_2groups(fit,aa,bb)
            idx=list(fit.params.index).index(term)
            var=float(cov[idx,idx])
            if var>0 and np.isfinite(var):
                se=float(np.sqrt(var))
                df=min(x.study_key.nunique(),x.species.nunique())-1
                crit=student_t.ppf(.975,df)
                out["two_way_publication_and_plant_cluster"]={
                   "SE":se,
                   "CI95_student_t":[float(fit.params[term]-crit*se),
                                     float(fit.params[term]+crit*se)],
                   "p_two_sided_t":float(2*student_t.sf(abs(float(fit.params[term])/se),df)),
                   "df":int(df)}
            else:out["two_way_publication_and_plant_cluster"]={
                "status":"NOT_ESTIMABLE_VARIANCE"}
        except (ValueError,np.linalg.LinAlgError) as e:
            out["two_way_publication_and_plant_cluster"]={
                "status":"NOT_ESTIMABLE","reason":str(e)[:180]}
    return out


def leverage(d,field):
    require(field in ("study_key","species"),"Unexpected entire-source deletion")
    vals=[];unestimable=0
    for label in d[field].unique():
        subset=d.loc[d[field].ne(label)]
        v=fitting(subset,SOURCE_MODEL,twoway=False)
        if v["status"]!="EXPLORATORY_ACTUAL_FIELD_REPRODUCTION_BLOOM_TIME_ASSOCIATION":
            unestimable+=1
            continue
        vals.append(v["beta_night_reproductive_share_fraction"])
    if not vals:return {"status":"NOT_ESTIMABLE_AFTER_ALL_SOURCE_DELETIONS"}
    arr=np.array(vals)
    return {
      "n_full_source_deletions":int(d[field].nunique()),
      "n_eligible_recomputed_source_deletions":len(vals),
      "n_not_estimable_deletions":unestimable,
      "beta_fraction_range":[float(arr.min()),float(arr.max())],
      "median":float(np.median(arr)),
      "fraction_positive":float((arr>0).mean()),
      "fraction_negative":float((arr<0).mean()),
      "no_publication_or_taxon_names_released":True,
      "highly_overlapping_deletions_not_independent_confirmatory_tests":True,
    }


def all_models(s):
    main=fitting(s,SOURCE_MODEL)
    require(main["status"]=="EXPLORATORY_ACTUAL_FIELD_REPRODUCTION_BLOOM_TIME_ASSOCIATION" and
            abs(main["beta_night_reproductive_share_fraction"]-BENCHMARK)<1e-7,
            "Original source 196 author-published paired reproduction bloom phenology baseline changed")
    # Keep the originally published study metadata and pair source unchanged:
    # These comparisons are observational across studies, not experimental bloom assignment.
    fit={
        "source_main_study_year_geography_elevation_metric":main,
        "source_20deg_geo_FE":fitting(s,SOURCE_MODEL+" + C(geo20)"),
        "source_only_native_studied_species":fitting(
            s[s.plant_native_source.eq("native")],SOURCE_MODEL),
        "source_original_coordinates_from_text_only":fitting(
            s[s.coordinate_source.eq("provided in text")],SOURCE_MODEL),
        "source_original_fruit_set_only":fitting(
            s[s.metric.eq("fruit set")],SOURCE_MODEL),
        "source_original_seed_set_only":fitting(
            s[s.metric.eq("seed set")],SOURCE_MODEL),
        "source_daynight_original_treatment_hours_reported":fitting(
            s[s.night_share_per_hour.notna()],SOURCE_MODEL),
        "source_bloom_x_original_priorPCNL":fitting(
            s,SOURCE_INTERACTION,term=LIGHT+":night_bloom"),
    }
    fam_support=int(s.groupby("family").night_bloom.nunique().eq(2).sum())
    if fam_support>=5:
        fit["source_botanical_family_FE"]=fitting(s,SOURCE_MODEL+" + C(family)")
    else:
        fit["source_botanical_family_FE"]={
           "status":"NOT_ESTIMABLE_INDEPENDENT_DAY_NIGHT_BLOOMED_FAMILIES",
           "mixed_bloom_families":fam_support}
    samepaper=s.groupby("study_key").night_bloom.nunique().eq(2)
    if int(samepaper.sum())>=20:
        fit["within_actual_published_study_FE"]=fitting(
            s,SOURCE_MODEL+" + C(study_key)")
    else:
        fit["within_actual_published_study_FE"]={
            "status":"NOT_ESTIMABLE_INDEPENDENT_MIXED_BLOOM_STUDIES",
            "studies_with_both_bloom_categories":int(samepaper.sum())}
    return fit


def run(joined,receipt,trait,out):
    out=Path(out)
    if out.exists():raise FileExistsError("Refuse altering immutable source-reproduction output")
    original,rec=load_source(joined,receipt)
    require(len(original)==301 and original.study_key.nunique()==118,
            "Real original source PCNL field pairing altered")
    merged=original.merge(original_trait_bloom(trait),on="species",how="left",
                          validate="many_to_one")
    cohort,support=valid_morphological_groups(merged)
    require((len(cohort),cohort.study_key.nunique(),cohort.species.nunique())==
             (196,73,74),"Original nightly and diurnal bloom group source changed")
    models=all_models(cohort)
    returnresult={
        "schema":"ALAN_true301_daynight_field_reproductive_bloom_phenology_independent_confirmation_v1",
        "source":{
           "original_field_pair_count":301,
           "original_studies":118,
           "original_real_prior_PCNl_years":[1992,2021],
           "original_source_sha256":rec["joined_file_sha256"],
           "author_trait_gitblob":"18381e550a7e405fe4f198909f0bb381bbecf1e0",
           "source_bloom_trial_populations":support,
        },
        "models":models,
        "full_source_publication_omissions":leverage(cohort,"study_key"),
        "full_source_plant_omissions":leverage(cohort,"species"),
        "interpretation_limits":{
           "experimental_day_night_reproductive_output_but_not_blooming_induced":True,
           "source_bloom_trait_study_compilation_not_manipulated":True,
           "family_or_study_fixed_effects_may_lack_support":True,
           "field_ALAN_not_experimentally_assigned":True,
           "not_night_moth_specific":True,
           "not_island_only_or_local_longitudinal_population":True,
           "actual_satellite_zero_not_ground_lux":True,
           "posthoc_multiple_stress_tests_not_confirmatory":True,
           "source_export_only_aggregate_no_rare_species_site_coordinates":True
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(returnresult,indent=2,ensure_ascii=False,
                              allow_nan=False)+"\n")
    print("ACTUAL_DAYNIGHT_BLOOM_FUNCTION_AND_ALAN_SOURCE_AUDIT",
        json.dumps({"support":support,
                    "models":{k:{"status":v.get("status"),
                        "beta":v.get("beta_night_reproductive_share_fraction"),
                        "ci":v.get("CI95_publication_cluster"),
                        "p":v.get("p_nominal"),
                        "two_way":v.get("two_way_publication_and_plant_cluster")}
                        for k,v in models.items()},
                    "study_loo":returnresult["full_source_publication_omissions"],
                    "plant_loo":returnresult["full_source_plant_omissions"]
                    },ensure_ascii=False),flush=True)
    return returnresult


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("joined","receipt","trait","out"):
        p.add_argument("--"+a,type=Path,required=True)
    args=p.parse_args()
    run(args.joined,args.receipt,args.trait,args.out)
