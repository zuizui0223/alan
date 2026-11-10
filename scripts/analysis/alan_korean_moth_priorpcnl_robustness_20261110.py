#!/usr/bin/env python3
"""Source-audited temporal, weight, leave-one-island and multiplicity checks.

No additional sampling or model fishing. Uses the same 377 original Korean
island light-trap samples and exactly year-minus-one PCNL. The original
primary 20-individual expected rarefied richness remains primary even if
10-individual sensitivity has a nominal p < 0.05.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

import numpy as np
from statsmodels.stats.multitest import multipletests

from scripts.analysis.alan_korean_moth_count_rarefied_priorpcnl_20261110 import (
    source_frame,fit,require,SCHEMA as MODEL_SCHEMA)
from scripts.analysis.alan_korean_moth_prior_year_pcnl_17year_20261110 import checksum

SCHEMA="ALAN_true_Korean_30island_moth_historic_PCNL_leaveislandout_and_FDR_v1"
MAIN_OUTCOMES=("rarefied20_expected_species",
               "log1p_moths_captured","log1p_species_richness",
               "rarefied10_expected_species","hill1_effective_species")

def stat_summary(vals):
    x=np.asarray(vals,dtype=float)
    x=x[np.isfinite(x)]
    if not len(x):return None
    return {"n_estimated":len(x),"minimum":float(x.min()),
            "median":float(np.median(x)),"maximum":float(x.max()),
            "fraction_negative":float(np.mean(x<0)),
            "fraction_positive":float(np.mean(x>0))}

def run(source,receipt,out):
    d,provenance=source_frame(source,receipt)
    primary={}
    for name in MAIN_OUTCOMES:
        primary[name]=fit(d,name,"island_equal",True)
    # Correction over fixed five source outcomes, not the large number of
    # alternative weights, thresholds and leave-one-island sensitivities.
    ps=[primary[v]["nominal_unadjusted_two_sided_p"] if
        primary[v]["status"]=="estimated_exploratory_original_field_data"
        else np.nan for v in MAIN_OUTCOMES]
    adjust={}
    if all(np.isfinite(ps)):
        _,q_bh,_,_=multipletests(ps,method="fdr_bh")
        _,q_holm,_,_=multipletests(ps,method="holm")
        adjust={v:{"nominal_p":float(p),"Benjamini_Hochberg_5_endpoints_q":float(q),
                   "Holm_familywise_5_endpoints_p":float(h)}
                for v,p,q,h in zip(MAIN_OUTCOMES,ps,q_bh,q_holm)}
    sensitivity={}
    for start in (2012,2017,2020):
        filtered=d.loc[d.year.ge(start)].copy()
        results={}
        for name in ("rarefied20_expected_species","rarefied10_expected_species",
                     "log1p_moths_captured"):
            v=fit(filtered,name,"island_equal",True)
            results[name]={
                "status":v["status"],
                "coefficient":v.get("beta_per_log1p_past_year_PCNL"),
                "95CI":v.get("95CI_nominal_cluster_islands"),
                "nominal_p":v.get("nominal_unadjusted_two_sided_p"),
                "real_source_islands":v.get("n_original_islands_clustered",
                                           v.get("support",{}).get("within_island_year_contrast_islands")),
                "real_source_site_years":v.get("support",{}).get("within_island_year_contrast_site_years"),
                "real_island_year_groups":v.get("support",{}).get("within_island_year_contrast_groups")}
        sensitivity["years_ge_"+str(start)]=results
    loo={}
    for name in ("rarefied20_expected_species","rarefied10_expected_species",
                 "log1p_moths_captured"):
        fits=[]
        for island in sorted(d.island.unique()):
            rem=d.loc[~d.island.eq(island)]
            res=fit(rem,name,"island_equal",True)
            if res["status"]=="estimated_exploratory_original_field_data":
                fits.append(res)
        loo[name]={
            "n_original_islands_deleted_one_at_a_time":int(d.island.nunique()),
            "n_estimable_omissions":len(fits),
            "coefficient_range":stat_summary([x["beta_per_log1p_past_year_PCNL"] for x in fits]),
            "number_nominal_CI_below_zero":int(sum(x["95CI_nominal_cluster_islands"][1]<0 for x in fits)),
            "number_nominal_CI_above_zero":int(sum(x["95CI_nominal_cluster_islands"][0]>0 for x in fits)),
            "do_not_select_omission_by_p_value":True}
    report={
        "schema":SCHEMA,"real_original_source_sha256":checksum(source),
        "source_observed_moth_survey_site_years":len(d),
        "source_independent_islands":int(d.island.nunique()),
        "source_actual_moth_specimens":int(d.original_moths_captured.sum()),
        "original_primary_20individual_result_preserved":True,
        "main_predeclared_five_original_field_outcomes_FDR_and_Holm":adjust,
        "historical_cohort_sensitivity":sensitivity,
        "all_island_leave_one_out":loo,
        "no_multiple_test_significance_selected":True,
        "source_real_moth_light_trap_count_not_true_population_density":True,
        "no_flower_visitation_or_effective_pollen_measurement":True,
        "spatial_and_temporal_confounding_unresolved":True,
        "trap_nights_and_month_missing":True}
    out=Path(out)
    require(not out.exists(),"Source-locked sensitivity report must not overwrite")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_ORIGINAL_KOREAN_MOTH_LIGHT_SOURCE_SENSITIVITY",
          json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--receipt",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.receipt,a.out)
