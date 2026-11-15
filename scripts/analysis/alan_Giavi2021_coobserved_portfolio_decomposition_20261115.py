#!/usr/bin/env python3
"""Decompose real Swiss-2016 ALAN pooled floral-cover visitation ratios by
co-recorded vs arm-only observed plant×site-pair×sampling-event opportunities.

Original 2016 Giavi et al. (2021) Zenodo 4540407 author data.Rdata:
MD5 f60b21d82c163cbb53520fda0bc61985.
Six true paired field assignments, 37 survey occasions, 21 taxa.
Unmatched observations are real source observations, not actual biological
absence in the other meadow. NEVER calculate a causal paired treatment
contrast from arm-only records or combine their flora as a matched taxon panel.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_Giavi2021_within_site_plant_visit_flower_temporal_20261115 import (
    source_grouping,site_bootstrap
)
from scripts.analysis.alan_Giavi2021_Swiss_original_daytime_guild_reanalysis_20261114 import (
    load_original_data,require
)
KEYS=["site","pair","plant"]
ARMS=["C","T"]

def snapshot(g):
    """Source row grouping exact, no invented missing source key to fill."""
    require(not g.duplicated(KEYS+["treatment"]).any(),
            "original plant x survey x meadow key duplicated")
    counts=g.groupby(KEYS,observed=True).treatment.nunique()
    require(counts.between(1,2).all(),"unrecognized treatment arms")
    key=counts.rename("n_source_arms").reset_index()
    v=g.merge(key,on=KEYS,how="inner",validate="many_to_one")
    require(len(v)==len(g),"source plant×occasion row loss")
    v["opportunity_support"]=v.n_source_arms.map({1:"arm_only_observed",2:"both_arms_observed"})
    require(v.opportunity_support.notna().all(),"invalid arm presence")
    return v,key

def totals(g):
    stats=g.groupby(["opportunity_support","treatment"],observed=True).agg(
        n_source_plant_observations=("plant","size"),
        source_flower_cover_proxy=("floral_cover","sum"),
        source_visit_events=("all3_visits","sum"),
        source_positive_visit_plant_events=("all3_visits",lambda x:int((x>0).sum()))
        ).reset_index()
    return {
      str(r.opportunity_support):{} for r in stats.itertuples(index=False)
    } | {
      support:{
        str(r.treatment):{
          "n_source_plant_observations":int(r.n_source_plant_observations),
          "source_flower_cover_proxy":float(r.source_flower_cover_proxy),
          "source_visit_events":float(r.source_visit_events),
          "source_positive_visit_plant_events":int(r.source_positive_visit_plant_events)
        }
        for r in stats[stats.opportunity_support.eq(support)].itertuples(index=False)
      } for support in stats.opportunity_support.unique()
    }

def density(t):
    C,T=t["C"],t["T"]
    require(C["source_flower_cover_proxy"]>0 and T["source_flower_cover_proxy"]>0,
            "positive observed floral cover for rate")
    rateC=C["source_visit_events"]/C["source_flower_cover_proxy"]
    rateT=T["source_visit_events"]/T["source_flower_cover_proxy"]
    return {"control_visits_per_cover_proxy":rateC,
            "lit_visits_per_cover_proxy":rateT,
            "lit_over_control_visitation_per_cover_proxy":rateT/rateC if rateC>0 else None}

def summarize_source(g):
    v,key=snapshot(g)
    by=totals(v)
    require(set(by)=={"both_arms_observed","arm_only_observed"},"source coverage classes drifted")
    main=by["both_arms_observed"];missing=by["arm_only_observed"]
    require(set(main)==set(ARMS) and set(missing)==set(ARMS),
            "source treatment arms missing")
    common=main["C"]["n_source_plant_observations"]
    require(common==main["T"]["n_source_plant_observations"],
            "paired observed plant opportunity counts differ")
    union=len(key)
    require(union==len(v)-common,"observed union key arithmetic")
    pairset=key.groupby(["site","pair"],observed=True).agg(
        unique_plant_observation_keys=("plant","size"),
        both_arms_observation_keys=("n_source_arms",lambda x:int((x==2).sum()))
    ).reset_index()
    pairset["source_key_jaccard"]=pairset.both_arms_observation_keys/pairset.unique_plant_observation_keys
    require(pairset.pair.nunique()==37 and pairset.site.nunique()==6,
            "actual 37 original paired surveys and six experimental sites required")
    per_site=pairset.groupby("site").source_key_jaccard.mean().sort_index()
    require(len(per_site)==6,"need six site pairs")
    base={
      arm:{
         "n_source_plant_observations":main[arm]["n_source_plant_observations"]+missing[arm]["n_source_plant_observations"],
         "source_flower_cover_proxy":main[arm]["source_flower_cover_proxy"]+missing[arm]["source_flower_cover_proxy"],
         "source_visit_events":main[arm]["source_visit_events"]+missing[arm]["source_visit_events"]
      } for arm in ARMS
    }
    matched_ratio=density(main)["lit_over_control_visitation_per_cover_proxy"]
    all_ratio=density(base)["lit_over_control_visitation_per_cover_proxy"]
    unmatched_ratio=density(missing)["lit_over_control_visitation_per_cover_proxy"]
    retention={arm:{
       "source_observed_plant_record_retention":main[arm]["n_source_plant_observations"]/base[arm]["n_source_plant_observations"],
       "source_flower_cover_retention":main[arm]["source_flower_cover_proxy"]/base[arm]["source_flower_cover_proxy"],
       "source_recorded_visits_retention":main[arm]["source_visit_events"]/base[arm]["source_visit_events"]
    } for arm in ARMS}
    # Exact algebra:
    # R_all/R_matched =
    # (flower retained in T / visits retained in T)
    # divided by (flower retained in C / visits retained in C).
    selection_multiplier=(
       retention["T"]["source_flower_cover_retention"]/retention["T"]["source_recorded_visits_retention"] /
       (retention["C"]["source_flower_cover_retention"]/retention["C"]["source_recorded_visits_retention"]))
    require(math.isclose(all_ratio,matched_ratio*selection_multiplier,rel_tol=1e-12),
            "exact multiplicative attrition decomposition failed")
    # Within each physical paired site, show change in both/all source
    # visitation-per-cover ratios, not fake independent experiment replicates.
    site_effect={}
    for site,w in v.groupby("site",observed=True):
        q=totals(w)
        require("both_arms_observed" in q,"site has no co-recorded flower sample")
        both=q["both_arms_observed"];armonly=q.get("arm_only_observed",{})
        allsite={arm:{
            "source_flower_cover_proxy":both[arm]["source_flower_cover_proxy"]+armonly.get(arm,{}).get("source_flower_cover_proxy",0.),
            "source_visit_events":both[arm]["source_visit_events"]+armonly.get(arm,{}).get("source_visit_events",0.)
        } for arm in ARMS}
        if min(allsite[a]["source_visit_events"] for a in ARMS)==0 or min(both[a]["source_visit_events"] for a in ARMS)==0:
            continue
        x=density(both)["lit_over_control_visitation_per_cover_proxy"]
        y=density(allsite)["lit_over_control_visitation_per_cover_proxy"]
        site_effect[site]={"paired_source_plant_opportunity_lit_dark_ratio":x,
                           "all_source_plant_opportunity_lit_dark_ratio":y,
                           "sign_change_due_to_source_coverage":bool((x-1)*(y-1)<0),
                           "log_ratio_difference_all_minus_matched":float(math.log(y/x))}
    require(len(site_effect)==6,"true paired site source ratio missing")
    sign_chg=sum(v["sign_change_due_to_source_coverage"] for v in site_effect.values())
    return {
      "schema":"ALAN_Giavi2021_Swiss_original_flora_source_support_compositional_decomposition_v1",
      "source_reproduction":{
        "article_DOI":"10.1038/s41467-021-22011-8",
        "original_Zenodo_DOI":"10.5281/zenodo.4540407",
        "source_data_Rdata_expected_md5":"f60b21d82c163cbb53520fda0bc61985",
        "original_insect_group_rows":1884,
        "source_plant_event_arm_rows":int(len(v)),
        "real_plant_species_count":int(v.plant.nunique()),
        "actual_paired_source_sampling_sessions":int(pairset.pair.nunique()),
        "actual_independent_light_dark_site_pairs":6},
      "source_support":{
        "exact_site_survey_plant_union":int(union),
        "both_arms_source_plant_survey_keys":int(common),
        "dark_only_observed_plant_survey_keys":int(missing["C"]["n_source_plant_observations"]),
        "lit_only_observed_plant_survey_keys":int(missing["T"]["n_source_plant_observations"]),
        "co_recorded_keys_over_union":float(common/union),
        "survey_event_source_co_recording_jaccard":{
          "n_actual_paired_surveys":int(len(pairset)),
          "mean_over_37_surveys":float(pairset.source_key_jaccard.mean()),
          "median_over_37_surveys":float(pairset.source_key_jaccard.median()),
          "min_over_37_surveys":float(pairset.source_key_jaccard.min()),
          "max_over_37_surveys":float(pairset.source_key_jaccard.max()),
          "mean_of_six_site_mean_jaccard":float(per_site.mean()),
          "site_bootstrap95_site_equal_jaccard_DESCRIPTIVE":site_bootstrap(per_site.to_numpy(float),20261115),
          "same_pair_survey_flowering_source_recording_not_census_of_true_plant_absence":True}},
      "source_observed_group_totals":{
        "both_arms_observed":main,
        "arm_only_observed":missing,
        "all_plant_source_observed":base},
      "sample_restriction_retention_by_arm":retention,
      "visit_density_ratios":{
        "all_source_observed":density(base),
        "matched_both_arm_observed":density(main),
        "arm_only_NOT_paired":density(missing),
        "all_over_matched_visitation_density_ratio_multiplier":float(all_ratio/matched_ratio),
        "exact_multiplicative_source_selection_factor":float(selection_multiplier),
        "all_ratio_reconstruction_error":float(all_ratio-matched_ratio*selection_multiplier),
        "co_observed_matched_ratio_above_one":bool(matched_ratio>1),
        "arm_only_unpaired_ratio_below_one":bool(unmatched_ratio<1),
        "all_source_ratio_above_one":bool(all_ratio>1)},
      "site_support_heterogeneity":{
        "independent_sites_compared":int(len(site_effect)),
        "site_pairs_where_ratio_sign_switches_all_vs_matched":int(sign_chg),
        "site_pairs_matched_density_positive":int(sum(x["paired_source_plant_opportunity_lit_dark_ratio"]>1 for x in site_effect.values())),
        "site_pairs_all_density_positive":int(sum(x["all_source_plant_opportunity_lit_dark_ratio"]>1 for x in site_effect.values())),
        "site_pair_log_ratio_diff_all_minus_matched_ANONYMIZED":[round(x["log_ratio_difference_all_minus_matched"],8) for x in site_effect.values()]},
      "biological_inference_gates":{
        "arm_only_observations_HAVE_NO_matched_treatment_counterfactual":True,
        "source_no_record_NE_not_observed_biological_absence":True,
        "co_recording_jaccard_NOT_flora_true_occurrence_similarity":True,
        "flower_cover_proxy_NOT_true_open_flower_or_pollen_supply":True,
        "post_treatment_floral_cover_adjustment_NOT_total_causal_light_effect":True,
        "seasonal_visitor_activity_hour_is_not_available":True,
        "P1xP3_fruit_seed_effect_NOT_identified":True,
        "independent_study_count":1,
        "no_nighttime_pollinator_services_or_island_shift_tested":True}
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    a=ap.parse_args();p=Path(a.out)
    if p.exists():raise FileExistsError("Refuse to overwrite study source receipt")
    source,receipt=load_original_data()
    g,design=source_grouping(source)
    z=summarize_source(g)
    require(receipt["original_md5"]==z["source_reproduction"]["source_data_Rdata_expected_md5"],
            "original Zenodo source md5 fail")
    require(design["original_source_plant_cover_event_rows"]==628 and
            z["source_reproduction"]["source_plant_event_arm_rows"]==628,
            "original plant opportunity count changed")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(z,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    key=z["visit_density_ratios"]
    print("GIAVI2016_REAL_FLORA_RECORD_SUPPORT_DECOMPOSITION",json.dumps({
        "source_support":z["source_support"],
        "retention":z["sample_restriction_retention_by_arm"],
        "ratios":key,
        "site_heterogeneity":z["site_support_heterogeneity"]
    },ensure_ascii=False))
if __name__=="__main__":main()
