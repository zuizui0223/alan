#!/usr/bin/env python3
"""2016 Swiss original paired plant×site×survey lit/dark visit allocation
controlling for recorded floral cover on the exact same survey occasion.

Within each matched opportunity (site,pair,plant), under a descriptive null of
equal visits per flower-cover unit: n_T | n = Binom(n,cover_T/(cover_T+cover_C)).
Aggregate OBSERVED-EXPECTED lit visitor counts across true six site pairs.
This is a conditionally exposure-standardized plant-opportunity association,
NOT a randomized-light total effect, pollen efficacy, seed, or phenology.
Never invent missing plants as zeros. No site/plant ID leaves the artifact.
"""
from __future__ import annotations
import argparse,json,itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_Giavi2021_within_site_plant_visit_flower_temporal_20261115 import (
 source_grouping, observed_common_events, site_bootstrap, exact_site_sign_flip, holm_by_family
)
from scripts.analysis.alan_Giavi2021_Swiss_original_daytime_guild_reanalysis_20261114 import (
  load_original_data,require
)

GUILDS={"all_three_documented_diurnal_orders":"all3_visits",
        "Hymenoptera":"Hymenoptera",
        "Diptera":"Diptera",
        "Coleoptera":"Coleoptera"}
MIN_COMMON_TIMES=(1,2,3,4)

def conditional_exposure_rows(common,guild):
    """One plant×real paired survey event, with both arms' actual floral cover.
    Conditional expected light visits E_L=(L+C)*f_L/(f_L+f_C).
    Zero visitor events are actual matched zero not NA.
    """
    c=common.copy()
    fL=c["floral_cover_T"].to_numpy(float)
    fC=c["floral_cover_C"].to_numpy(float)
    vL=c[f"{guild}_T"].to_numpy(float)
    vC=c[f"{guild}_C"].to_numpy(float)
    if len(c)==0 or not (np.isfinite(fL).all() and np.isfinite(fC).all()
        and np.isfinite(vL).all() and np.isfinite(vC).all()):
        raise ValueError("source nonfinite coverage")
    if (fL<=0).any() or (fC<=0).any() or (vL<0).any() or (vC<0).any():
        raise ValueError("zero florals or negative visits; source scope changed")
    n=vL+vC; share=fL/(fL+fC)
    c["total_visits_exact_pair"]=n
    c["actual_lit_visits"]=vL
    c["expected_lit_visits_cover_proportional"]=n*share
    c["actual_dark_visits"]=vC
    c["expected_dark_visits_cover_proportional"]=n*(1-share)
    c["exposure_balance_lit_share"]=share
    c["lit_excess_visits_beyond_flower_share"]=vL-n*share
    c["zero_visit_pair"]=n==0
    return c

def per_site(q,threshold):
    """Keep within plant×site series having k matched real survey occasions.
    Not all missing plant×site×occasion = true flower absence.
    """
    n=q.groupby(["site","plant"],observed=True).pair.nunique()
    eligible=n[n>=threshold].index
    mask=pd.MultiIndex.from_frame(q[["site","plant"]]).isin(eligible)
    e=q.loc[mask]
    rows=[]
    for site,g in e.groupby("site",observed=True):
        visits=float(g.total_visits_exact_pair.sum())
        if visits<=0:continue
        expected=float(g.expected_lit_visits_cover_proportional.sum())
        observed=float(g.actual_lit_visits.sum())
        a=observed/visits
        x=expected/visits
        rows.append({"site":site,"n_matched_plant_occasions":int(len(g)),
          "n_matched_plant_series":int(g.plant.nunique()),
          "n_recorded_visits_over_both_arms":visits,
          "n_observed_lit_visits":observed,"n_expected_lit_visits_from_flower_cover":expected,
          "observed_lit_share":a,"cover_proportional_expected_lit_share":x,
          "observed_minus_expected_lit_visit_fraction":a-x,
          "matched_zero_visit_pairs":int(g.zero_visit_pair.sum())})
    return rows,e

def summarize(rows,matched,seed):
    if not rows:return {
        "sites":0,"eligible_source_matched_plant_opportunity_rows":int(len(matched)),
        "no_lighting_inference_fitted":True}
    series=pd.DataFrame(rows)
    d=series.observed_minus_expected_lit_visit_fraction.to_numpy(float)
    loocv=[float(np.delete(d,i).mean()) for i in range(len(d))] if len(d)>=2 else []
    overall_visits=float(series.n_recorded_visits_over_both_arms.sum())
    observed_pooled=float(series.n_observed_lit_visits.sum())
    expected_pooled=float(series.n_expected_lit_visits_from_flower_cover.sum())
    return {
      "sites":len(series),
      "source_matched_plant_opportunity_rows":len(matched),
      "source_matched_plant_series":int(matched.groupby(["site","plant"],observed=True).ngroups),
      "source_plant_taxa":int(matched.plant.nunique()),
      "source_matched_zero_visitation_rows":int(matched.zero_visit_pair.sum()),
      "source_total_real_visits_both_arms":float(matched.total_visits_exact_pair.sum()),
      "source_true_lit_visits":float(matched.actual_lit_visits.sum()),
      "source_true_dark_visits":float(matched.actual_dark_visits.sum()),
      "site_equal_mean_excess_lit_visitor_share_over_flower_cover":float(d.mean()),
      "site_equal_median_excess_lit_visitor_share":float(np.median(d)),
      "leave_one_true_site_out_site_equal_means":loocv,
      "leave_one_out_changes_mean_effect_sign":bool(
            loocv and min(loocv)<0<max(loocv)),
      "pooled_matched_event_expected_lit_share":float(expected_pooled/overall_visits),
      "pooled_matched_event_observed_lit_share":float(observed_pooled/overall_visits),
      "pooled_matched_event_excess_lit_visit_fraction":float(
            (observed_pooled-expected_pooled)/overall_visits),
      "sites_with_excess_lit_visits_positive":int(np.sum(d>0)),
      "sites_with_excess_lit_visits_negative":int(np.sum(d<0)),
      "site_equal_bootstrap95_descriptive":site_bootstrap(d,seed),
      "site_exact_signflip_descriptive":exact_site_sign_flip(d),
      "true_site_deltas_anonymized":[float(v) for v in d],
      "visitor_count_conditions_on_recorded_matched_plant_events":True,
      "flower_cover_is_post_treatment_proxy_and_may_mediate_total_ALAN_effect":True,
      "null_assumption":"each recorded insect visit assigned to lit arm in proportion to its same-plant same-event measured floral cover, conditional on total visits (Binomial null)",
      "null_is_not_original_randomized_lamp_assignment":True
    }

def analyze():
    raw,receipt=load_original_data()
    all_source,design=source_grouping(raw)
    common=observed_common_events(all_source)
    require(len(common)==193,"exact source 193 matched plant×site×occasion records changed")
    count_by_arm=all_source.groupby("treatment",observed=True).agg(
        n_plant_events=("plant","size"),observed_visits=("all3_visits","sum"),
        floral_cover=("floral_cover","sum")).to_dict("index")
    matched_arms={
        arm:{"n_plant_events":int(len(common)),
         "visits_all_three_orders":float(common[f"all3_visits_{arm}"].sum()),
         "floral_cover_proxy_sum":float(common[f"floral_cover_{arm}"].sum())}
        for arm in ["C","T"]}
    # A pooled ratio of visit counts / sum of floral cover is an intentionally
    # different estimand from the within plant×pair conditional expected share.
    def pooled_rate_ratio(ar,controls):
        if min(float(ar["T"]["floral_cover_proxy_sum"]),float(ar["C"]["floral_cover_proxy_sum"]))<=0:
            raise ValueError("original flower cover proxy cannot have nonpositive totals")
        return float(
            (ar["T"]["visits_all_three_orders"]/ar["T"]["floral_cover_proxy_sum"])/
            (ar["C"]["visits_all_three_orders"]/ar["C"]["floral_cover_proxy_sum"]))
    pooled_matched=pooled_rate_ratio(matched_arms,common)
    pooled_all=(float(count_by_arm["T"]["observed_visits"])/float(count_by_arm["T"]["floral_cover"])) / (float(count_by_arm["C"]["observed_visits"])/float(count_by_arm["C"]["floral_cover"]))
    summary={}
    for i,(name,col) in enumerate(GUILDS.items()):
        v=conditional_exposure_rows(common,col)
        blocks={}
        for j,k in enumerate(MIN_COMMON_TIMES):
            site,eligible=per_site(v,k)
            blocks[str(k)]=summarize(site,eligible,20261115+i*50+j)
        summary[name]=blocks
    primary={k:v["1"]["site_exact_signflip_descriptive"]
               ["two_sided_site_signflip_fraction_DESCRIPTIVE"]
             for k,v in summary.items() if v["1"].get("sites",0)>=3}
    # Real 6 paired original light sites; pair-specific logged statistics
    # are diagnostic, and four guilds plus four thresholds remain exploratory.
    return {
      "schema":"ALAN_Giavi2021_exact_plant_event_matched_light_cover_conditioned_visits_v1",
      "paper_DOI":"10.1038/s41467-021-22011-8",
      "author_Zenodo_DOI":"10.5281/zenodo.4540407",
      "source_author_original_MD5":receipt["original_md5"],
      "source_author_original_SHA256":receipt["sha256"],
      "study_design":{
        **design,
        "independent_lit_dark_true_site_pairs":6,
        "paired_survey_occasions":37,
        "both_treatment_same_plant_same_survey_matched_plant_occasions":193,
        "source_observed_all_plant_opportunity_by_arm":count_by_arm,
        "source_matched_plant_opportunity_by_arm":matched_arms,
        "mismatched_plant_event_opportunities_not_true_absence":True},
      "visitation_rate_per_measured_flower_cover_pooled_diagnostic":{
         "all_628_plant_opportunity_records_lit_over_dark":pooled_all,
         "matched_193_plant_opportunities_lit_over_dark":pooled_matched,
         "both_are_aggregate_count_per_cover_not_within_event_conditional_effect":True,
         "pooled_matched_excess_lit_share_vs_6site_equal_excess_lit_share_sign_mismatch":
           bool((pooled_matched>1 and summary["all_three_documented_diurnal_orders"]["1"]["site_equal_mean_excess_lit_visitor_share_over_flower_cover"]<0) or
                (pooled_matched<1 and summary["all_three_documented_diurnal_orders"]["1"]["site_equal_mean_excess_lit_visitor_share_over_flower_cover"]>0)),
         "source_matched_fraction_by_treatment_of_observed_visit_counts":{
            arm:matched_arms[arm]["visits_all_three_orders"]/float(count_by_arm[arm]["observed_visits"]) for arm in ["C","T"]},
         "source_matched_fraction_by_treatment_of_observed_flower_cover":{
            arm:matched_arms[arm]["floral_cover_proxy_sum"]/float(count_by_arm[arm]["floral_cover"]) for arm in ["C","T"]},
         "observed_lit_dark_rates_not_causal_total_ALAN_effect":True
      },
      "conditional_visit_share_exposure_adjustment":summary,
      "four_guild_primary_threshold1_Holm_exploratory":holm_by_family(primary),
      "four_matched_occasion_thresholds_exploratory_not_multiplicity_corrected":True,
      "no_hourly_pollen_vector_flower_timing_recorded":True,
      "no_seed_or_effective_conspecific_pollen_measured":True,
      "not_a_plant_total_fitness_ALAN_response":True,
      "flower_cover_post_treatment_adjustment_changes_estimand_from_total_effect":True,
      "site_means_not_independent_pub_count":True,
      "true_h3_phenology_mismatch_identified":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    a=ap.parse_args();p=Path(a.out)
    if p.exists():raise FileExistsError("No overwriting source evidence")
    z=analyze()
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(z,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    print("GIAVI2021_REAL_COVER_CONDITIONED_LIT_VISIT_SHARE",json.dumps({
      "matched":z["study_design"]["both_treatment_same_plant_same_survey_matched_plant_occasions"],
      "full":z["study_design"]["source_observed_all_plant_opportunity_by_arm"],
      "matched_arms":z["study_design"]["source_matched_plant_opportunity_by_arm"],
      "estimates":{k:{j:{key:x.get(key) for key in ["sites","source_matched_plant_opportunity_rows","source_matched_zero_visitation_rows","site_equal_mean_excess_lit_visitor_share_over_flower_cover","site_equal_bootstrap95_descriptive","sites_with_excess_lit_visits_positive"]}
                       for j,x in v.items()} for k,v in z["conditional_visit_share_exposure_adjustment"].items()},
      "Holm":z["four_guild_primary_threshold1_Holm_exploratory"],
      "pooled_and_site_estimates":z["visitation_rate_per_measured_flower_cover_pooled_diagnostic"],
      "site_L1_leaveout":z["conditional_visit_share_exposure_adjustment"]["all_three_documented_diurnal_orders"]["1"]["leave_one_true_site_out_site_equal_means"]
      },ensure_ascii=False))
if __name__=="__main__":main()
