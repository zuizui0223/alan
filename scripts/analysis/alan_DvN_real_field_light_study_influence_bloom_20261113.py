#!/usr/bin/env python3
"""Source-authenticated ORIGINAL 301 day/night reproduction vs historic ALAN:
  (1) exhaustive publication/species source influence;
  (2) real lit/dark geographic common support within independent papers;
  (3) independent author species bloom-time trait and PCNL cross-context gate.

This is NOT experimental ALAN, not moth-only effective service, not an island
population study. The published day/night-reproduction outcome IS experimentally
compared within study×species×location×season×year×metric, whereas historical
PCNL differs primarily BETWEEN observational field study locations.
No original per-species or field-site records are exported.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import t

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, estimate, LIGHT, require
)

TRAITS_BLOB="18381e550a7e405fe4f198909f0bb381bbecf1e0"
BASELINE=.0030676
BASELINE_FORMULA=("night_share ~ "+LIGHT+
   " + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)")
ORIGINAL_STUDIES=118
ORIGINAL_PAIRS=301
ORIGINAL_PLANTS=123


def original_trait_bloom(trait_path:Path):
    data=Path(trait_path).read_bytes()
    blob=hashlib.sha1(("blob "+str(len(data))+"\0").encode()+data).hexdigest()
    require(blob==TRAITS_BLOB,"Original author plant trait file changed")
    trait=pd.read_csv(trait_path,dtype=str,keep_default_na=False)
    require(len(trait)==142 and trait.plant_species.is_unique,
            "Author 142 species independent bloom traits source changed")
    subset=trait[["plant_species","bloom_period_simple","breeding_system","plant_family"]].copy()
    subset=subset.rename(columns={"plant_species":"species"})
    subset["species"]=subset.species.astype(str).str.split().str.join(" ")
    subset["bloom_period_simple"]=subset.bloom_period_simple.replace(
        {"NA":"unknown","":"unknown"}).str.lower()
    return subset


def exhaustively_remove(d:pd.DataFrame,group:str):
    """Each case deletes an entire original published study or plant species.
    Do not select only omissions that enhance the apparent ALAN signal."""
    require(group in ("study_key","species"),"Source omission group invalid")
    base=estimate(d,"frozen_original",BASELINE_FORMULA)
    require(base.get("status")=="EXPLORATORY_FIELD_STUDY_NIGHT_VS_DAY_REPRODUCTIVE_OUTPUT_PRIOR_PCNL" and
            abs(base["pcnl_prior_light_coefficient_fraction_per_mean_log1p"]-BASELINE)<.00003,
            "Original 301 source and study-equal baseline light estimate drift")
    rows=[]
    failed={}
    for tag in d[group].drop_duplicates():
        subset=d.loc[d[group].ne(tag)].copy()
        model=estimate(subset,"exclude_one_"+group,BASELINE_FORMULA,two_way=False)
        if model["status"]!="EXPLORATORY_FIELD_STUDY_NIGHT_VS_DAY_REPRODUCTIVE_OUTPUT_PRIOR_PCNL":
            failed[model["status"]]=failed.get(model["status"],0)+1
            continue
        rows.append({"beta":float(model["pcnl_prior_light_coefficient_fraction_per_mean_log1p"]),
                     "p_nominal":float(model["p_nominal_study_cluster"]),
                     "ci":model["95CI_study_cluster"]})
    beta=np.array([z["beta"] for z in rows],float)
    require(len(beta)>0,"No original study/type deletion estimable")
    return {
       "original_omission_units":int(d[group].nunique()),
       "n_estimable_after_deletion":len(rows),
       "n_design_failure":int(sum(failed.values())),
       "n_no_source_support_failure_by_status":failed,
       "beta_fraction_per_prior_logPCNL_range":[float(beta.min()),float(beta.max())],
       "beta_median":float(np.median(beta)),
       "fraction_negative":float((beta<0).mean()),
       "fraction_positive":float((beta>0).mean()),
       "n_nominal_CI_excluding_zero":int(sum(
           (v["ci"][0]>0 or v["ci"][1]<0) for v in rows)),
       "smallest_nominal_p":float(min(v["p_nominal"] for v in rows)),
       "all_publication_species_identities_suppressed":True,
       "interpretation":"Repeated highly correlated source deletions; not independent validations",
    }


def actual_comparable_light(d:pd.DataFrame):
    """Exposure positivity in ACTUAL independent experiments rather than
    pretending observed dark PCNL means absence of lights at individual flowers."""
    x=d.copy()
    x["zero_prior_PCNl"]=x[LIGHT].eq(0)
    x["lit_prior_PCNl"]=x[LIGHT].gt(0)
    x["geo20"]=((x.lat+90)//20).astype(int).astype(str)+"_"+(
                ((x.lon+180)//20).astype(int).astype(str))
    x["geo10"]=((x.lat+90)//10).astype(int).astype(str)+"_"+(
                ((x.lon+180)//10).astype(int).astype(str))
    def spatial(col):
        block=[]
        for _,g in x.groupby(col,observed=True):
            dark=g.loc[g.zero_prior_PCNl]
            light=g.loc[g.lit_prior_PCNl]
            if dark.empty or light.empty:continue
            studies0=dark.study_key.nunique();studies1=light.study_key.nunique()
            block.append({
                "dark_studies":studies0,"light_studies":studies1,
                "dark_pairs":len(dark),"light_pairs":len(light),
                "both_arms_at_least_2_published_studies":studies0>=2 and studies1>=2,
                "both_arms_at_least_4_published_studies":studies0>=4 and studies1>=4,
            })
        return {
           "all_source_regions":int(x[col].nunique()),
           "regions_containing_both_zero_and_positive_prior_PCNL":len(block),
           "regions_both_arms_at_least_2_independent_publications":sum(
               z["both_arms_at_least_2_published_studies"] for z in block),
           "regions_both_arms_at_least_4_independent_publications":sum(
               z["both_arms_at_least_4_published_studies"] for z in block),
           "study_pair_count_in_both_strata_regions":sum(
               z["dark_pairs"]+z["light_pairs"] for z in block),
           "zero_is_satellite_source_mean_not_measured_garden_lux":True}
    def withinstudy():
        info=x.groupby("study_key",observed=True).agg(
           n_pairs=("source_pair_id","size"),
           n_orig_sites=("lat","nunique"),
           exposure_span=(LIGHT,lambda z:float(z.max()-z.min())),
           zero_states=("zero_prior_PCNl","nunique"))
        return {
           "total_source_studies":int(len(info)),
           "studies_with_repeated_same_source_paper_outcomes":int(info.n_pairs.ge(2).sum()),
           "studies_with_measurable_real_light_span_greater_than_003":int(
               info.exposure_span.gt(.03).sum()),
           "studies_with_both_zero_and_positive_prior_satellite_PCNL":int(
               info.zero_states.eq(2).sum()),
           "source_studies_sufficient_for_within_study_PCNL_effect":False,
           "no_within_study_light_manipulation":True}
    return {
       "reproductive_experiment_pairs_original":len(x),
       "original_independent_studies":int(x.study_key.nunique()),
       "source_light_zero_pairs":int(x.zero_prior_PCNl.sum()),
       "source_light_positive_pairs":int(x.lit_prior_PCNl.sum()),
       "source_light_zero_studies":int(x.loc[x.zero_prior_PCNl,"study_key"].nunique()),
       "source_light_positive_studies":int(x.loc[x.lit_prior_PCNl,"study_key"].nunique()),
       "10degree":spatial("geo10"),
       "20degree":spatial("geo20"),
       "within_real_published_study":withinstudy(),
    }


def source_bloom_species_and_ALAN(d:pd.DataFrame,traits:pd.DataFrame):
    """A priori biological negative-control opportunity: independent original
    botanical bloom-time category should correspond to measured relative
    night-only reproductive benefit, while ALAN×bloom remains observational.

    The author's species bloom trait is NOT a night moth identity.
    """
    x=d.merge(traits,on="species",how="left",validate="many_to_one")
    x["bloom_period_simple"]=x.bloom_period_simple.fillna("unknown")
    tab=x.groupby("bloom_period_simple",observed=True).agg(
       paired_original_reproductive_outcomes=("night_share","size"),
       studies=("study_key","nunique"),
       species=("species","nunique"))
    support={str(k):{
         "paired_original_reproductive_outcomes":int(v.paired_original_reproductive_outcomes),
         "original_studies":int(v.studies),
         "plant_species":int(v.species)}
      for k,v in tab.iterrows()}
    eligible=x[x.bloom_period_simple.isin(("day","night"))].copy()
    night=eligible.bloom_period_simple.eq("night")
    est={"status":"NOT_ESTIMABLE_INDEPENDENT_BLOOM_PHENOTYPE_STUDIES",
         "n_pairs":len(eligible),"n_studies":int(eligible.study_key.nunique()),
         "night_bloom_source_pairs":int(night.sum()),
         "day_bloom_source_pairs":int((~night).sum())}
    studies_by=eligible.groupby("bloom_period_simple").study_key.nunique()
    if (len(eligible)>=95 and eligible.study_key.nunique()>=60 and
        len(studies_by)==2 and int(studies_by.min())>=25 and
        eligible[LIGHT].std()>.02):
        eligible["night_bloom"]=night.astype(int)
        eligible["study_equal_weight"]=1/eligible.groupby("study_key").study_key.transform("size")
        # Main-source covariates frozen from original 301-site-year model.
        no_interaction=("night_share ~ "+LIGHT+" + night_bloom"+
            " + exp_year + lat + lat2 + lon + lon2 + elevation_km + C(metric)")
        with_interaction=no_interaction.replace(LIGHT+" + night_bloom",
                            LIGHT+"*night_bloom")
        fit={}
        for key,formula,term in [
           ("external_night_bloom_vs_day_bloom_reproductive_main",
            no_interaction,"night_bloom"),
           ("night_bloom_differential_historical_ALAN_slope",
            with_interaction,LIGHT+":night_bloom")]:
            try:
                q=smf.wls(formula,data=eligible,
                     weights=eligible.study_equal_weight).fit(
                     cov_type="cluster",cov_kwds={"groups":eligible.study_key},use_t=False)
                if (np.linalg.matrix_rank(q.model.exog)!=len(q.params) or
                    term not in q.params or not np.isfinite(q.bse[term])):
                    fit[key]={"status":"NOT_ESTIMABLE_DESIGN_RANK"}
                    continue
                ci=q.conf_int().loc[term]
                fit[key]={
                    "status":"EXPLORATORY_SOURCE_DAY_NIGHT_BLOOMED_PLANT_STUDY_COMPARISON",
                    "beta_night_share_fraction":float(q.params[term]),
                    "95CI_publication_cluster":[float(ci.iloc[0]),float(ci.iloc[1])],
                    "nominal_p":float(q.pvalues[term]),
                    "n_original_reproductive_pairs":len(eligible),
                    "independent_studies":int(eligible.study_key.nunique()),
                    "n_botanical_species":int(eligible.species.nunique()),
                    "term":term,
                    "original_source_studies_not_randomly_sampled":True
                }
            except (ValueError,np.linalg.LinAlgError) as e:
                fit[key]={"status":"NOT_ESTIMABLE_MODEL_FAILURE",
                          "reason":str(e)[:180]}
        est={"status":"SOURCE_MEASURED_BLOOM_AND_REAL_REPRODUCTIVE_COMPARISON",
             "source_study_support":support,"fits":fit}
    else:est["source_study_support"]=support
    return est


def run(joined:Path,receipt:Path,trait:Path,out:Path):
    if out.exists():raise FileExistsError("Do not overwrite authentic source audit")
    x,src=load_source(joined,receipt)
    require(len(x)==ORIGINAL_PAIRS and
            x.study_key.nunique()==ORIGINAL_STUDIES and
            x.species.nunique()==ORIGINAL_PLANTS,
            "Original genuine paired field reproductive source drift")
    original=estimate(x,"ORIGINAL_BENCHMARK",BASELINE_FORMULA)
    require(original["status"]=="EXPLORATORY_FIELD_STUDY_NIGHT_VS_DAY_REPRODUCTIVE_OUTPUT_PRIOR_PCNL" and
            abs(float(original["pcnl_prior_light_coefficient_fraction_per_mean_log1p"])-BASELINE)<.00003,
            "Parent PR51 original 301 actual field outcome regression changed")
    outputs={
        "schema":"ALAN_global301_daynight_experiment_field_light_study_influence_bloom_v1",
        "status":"AUTHENTIC_PREEXPERIMENT_PCNL_ORIGINAL_PAIRED_FIELD_STUDY_AUDIT",
        "source":{
           "doi":"10.1111/ele.70036",
           "study_treatment_original_gitblob":"39a82d966afcdde670603add1c881ace8791ffe9",
           "plant_traits_original_gitblob":TRAITS_BLOB,
           "original_field_reproductive_pairs":ORIGINAL_PAIRS,
           "original_independent_studies":ORIGINAL_STUDIES,
           "original_botanical_species":ORIGINAL_PLANTS,
           "original_30_MD5_annual_raster_years":[1992,2021],
           "original_parent_joined_source_SHA256":src["joined_file_sha256"],
           "original_parent_model_slope":original,
        },
        "original_research_design_exposure_common_support":actual_comparable_light(x),
        "all_one_at_a_time_publication_omissions":exhaustively_remove(x,"study_key"),
        "all_one_at_a_time_plant_species_omissions":exhaustively_remove(x,"species"),
        "original_author_botanical_bloomtime_validation":source_bloom_species_and_ALAN(
             x,original_trait_bloom(trait)),
        "limits":{
            "does_not_establish_an_effect_of_manipulated_ALAN":True,
            "source_PCNl_1992_prior_year_is_not_flower_lux":True,
            "night_exclusion_allows_non_moth_pollinators":True,
            "study_cluster_not_301_independent_published_experiments":True,
            "whole_source_LOO_point_signs_not_multiple_independent_studies":True,
            "original_botanical_bloomtime_not_local_field_moth_identification":True,
            "selected_source_study_populations_not_wild_island_only":True,
            "no_original_author_study_site_species_rows_exported":True,
        }
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(outputs,sort_keys=True,indent=2,ensure_ascii=False,
                              allow_nan=False)+"\n")
    compact={
       "study_deletion":outputs["all_one_at_a_time_publication_omissions"],
       "plant_deletion":outputs["all_one_at_a_time_plant_species_omissions"],
       "source_exposure":outputs["original_research_design_exposure_common_support"],
       "bloom":outputs["original_author_botanical_bloomtime_validation"]}
    print("ALAN_AUTHENTIC_DIEL_301_STUDY_INFLUENCE_AND_BLOOM",
          json.dumps(compact,ensure_ascii=False),flush=True)
    return outputs


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("joined","receipt","trait","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.joined,a.receipt,a.trait,a.out)
