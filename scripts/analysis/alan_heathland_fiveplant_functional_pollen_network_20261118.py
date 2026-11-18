#!/usr/bin/env python3
"""Source-pinned test: visitation network vs functional pollen network in heath.

Ballantyne et al. (2015), original 1338 single visits to virgin stigmas in
one Dorset heathland, five plants including EIV=8.1 ALAN-selected Erica
tetralix (from the frozen alan/main wet flora case list).

Each observed visit has a source-verified conspecific pollen deposition count.
Repeated visits and pollinator taxa are correlated within source survey
dates. No ALAN manipulation, no moth visit in this source, no seed set.

Calculate for each SAME plant:
 - whether the most frequent visitor was also the taxon contributing the most
   *actually measured* conspecific pollen.
 - which rare taxa (<5% of flower visits) deliver >=10% of pollen deposited.
 - how different the first-visitor share is from first-pollen supplier share.
 - the hypothetical one-taxon loss given a frozen observed pollen network.
 - posthoc day-cluster bootstrap support for those source ranking differences,
   not independent treatment or meta trials.
 - hypothetical one-for-one 10% BEES-to-FLIES substitution only when
   within-plant fly source has >=10 visits and >=4 original calendar days.
No fabricated moth/stigma calibration, ALAN coefficients or pollinator
abundances from visitation frequency.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
import random
import statistics
from pathlib import Path

from scripts.analysis.alan_ballantyne2015_heathland_visit_pollen_effectiveness_20261118 import (
  load,normalize,PLANTS,SEED
)

SCHEMA="alan_heathland_2015_original_1338_functional_network_stability_v1"
BOOT=3999

def summarize_plant(records,plant):
    source=[x for x in records if x["plant"]==plant]
    if not source:return None
    n=len(source)
    by=defaultdict(list)
    for r in source:by[r["taxon"]].append(r)
    dep={k:sum(r["conspecific_pollen"] for r in vals) for k,vals in by.items()}
    visits={k:len(vals) for k,vals in by.items()}
    total=sum(dep.values())
    assert total>0
    main_visit=max(visits,key=lambda x:(visits[x],x))
    main_func=max(dep,key=lambda x:(dep[x],x))
    rare=[{"source_taxon":key,
           "visits":count,
           "source_visit_fraction":count/n,
           "original_total_conspecific_pollen":dep[key],
           "source_pollen_fraction":dep[key]/total,
           "pollen_share_to_visitation_share":(dep[key]/total)/(count/n)}
          for key,count in visits.items() if count/n <0.05 and dep[key]/total>=0.10]
    rare.sort(key=lambda x:-x["source_pollen_fraction"])
    function={
      "plant":plant,
      "source_virgin_stigma_visits":n,
      "source_unique_visit_dates":len({r["date"] for r in source}),
      "source_visitor_taxa":len(visits),
      "source_stigma_conspecific_pollen_total":total,
      "top_source_visitation_taxon":main_visit,
      "top_source_functional_pollen_depositor":main_func,
      "source_same_taxon_tops_both":main_visit==main_func,
      "top_visitor_visitation_fraction":visits[main_visit]/n,
      "top_visitor_functional_pollen_fraction":dep[main_visit]/total,
      "top_functional_taxon_visitation_fraction":visits[main_func]/n,
      "top_functional_taxon_pollen_fraction":dep[main_func]/total,
      "source_pollen_miss_if_top_visitation_taxon_missing":
          dep[main_visit]/total,
      "source_pollen_miss_if_top_pollen_taxon_missing":
          dep[main_func]/total,
      "original_rare_under5pct_visits_over10pct_pollen_taxa":rare,
      "source_loss_counterfactual_with_frozen_other_visitation":True,
      "not_plant_fitness_and_not_manipulated_ALAN":True,
    }
    for comparator in ("bees","flies"):
        subset=[x for x in source if x["guild"]==comparator]
        function[comparator+"_original_visits"]=len(subset)
        function[comparator+"_visit_source_days"]=len({x["date"] for x in subset})
        function[comparator+"_per_visit_mean_conspecific_stigma_pollen"]=(
           statistics.fmean(x["conspecific_pollen"] for x in subset)
           if subset else None)
    bee=function["bees_per_visit_mean_conspecific_stigma_pollen"]
    fly=function["flies_per_visit_mean_conspecific_stigma_pollen"]
    nfly=function["flies_original_visits"]
    nflydays=function["flies_visit_source_days"]
    if (nfly>=10 and nflydays>=4 and bee is not None and
        bee>0 and fly is not None):
        ratio=fly/bee
        function["hypothetical_10pct_bee_visits_replaced_with_fly_visits"]={
           "status":"source_support_pass_conditional_only",
           "assumptions":"Exactly 10pct of observed bee visits replaced 1:1 by same-plant fly visits; stable pollen per visit and flowers; no other interactions; observational only",
           "fly_to_bee_same_plant_virgin_stigma_mean_ratio":ratio,
           "hypothetical_change_in_bee_fly_contributed_pollen_fraction":
               .10*(ratio-1),
           "not_moth_or_ALAN_pollinator_replacement":True,
        }
    else:
        function["hypothetical_10pct_bee_visits_replaced_with_fly_visits"]={
            "status":"source_fly_within_plant_replication_STOP",
            "minimum_orig_fly_visits":10,
            "minimum_source_fly_days":4,
            "observed_orig_fly_visits":nfly,
            "observed_orig_fly_days":nflydays,
            "not_moth_or_ALAN_pollinator_replacement":True
        }
    return function

def day_cluster_bootstrap(records,plant):
    source=[x for x in records if x["plant"]==plant]
    days=defaultdict(list)
    for r in source:days[r["date"]].append(r)
    keys=list(days)
    if len(keys)<5:
        return {"status":"fewer_than_five_observed_source_days_STOP",
                "original_dates":len(keys)}
    rng=random.Random(SEED+sum(ord(x) for x in plant)+424)
    mismatch=0
    top_func_sums=0.
    top_visit_sums=0.
    rare_count=0
    for _ in range(BOOT):
        s=[r for day in rng.choices(keys,k=len(keys)) for r in days[day]]
        summary=summarize_plant(s,plant)
        if summary is None:continue
        mismatch+=not summary["source_same_taxon_tops_both"]
        top_func_sums+=summary["source_pollen_miss_if_top_pollen_taxon_missing"]
        top_visit_sums+=summary["source_pollen_miss_if_top_visitation_taxon_missing"]
        rare_count+=len(summary["original_rare_under5pct_visits_over10pct_pollen_taxa"])>0
    return {
      "status":"descriptive_source_day_cluster_bootstrap",
      "original_source_days":len(keys),
      "source_bootstrap_resamples":BOOT,
      "fraction_resamples_where_top_visitation_vs_pollen_taxon_differ":mismatch/BOOT,
      "mean_bootstrap_most_functional_taxon_pollen_share":top_func_sums/BOOT,
      "mean_bootstrap_top_visitation_taxon_pollen_share":top_visit_sums/BOOT,
      "fraction_bootstraps_with_rare_taxon_under5pct_visits_above10pct_pollen":rare_count/BOOT,
      "this_is_bootstrap_stability_not_a_causal_light_intervention_test":True
    }

def execute():
    source=normalize(load())
    assert len(source)==1338
    results={}
    for plant in PLANTS:
        direct=summarize_plant(source,plant)
        direct["source_day_cluster_bootstrap_stability"]=day_cluster_bootstrap(source,plant)
        results[plant]=direct
    mismatched=[p for p,x in results.items() if not x["source_same_taxon_tops_both"]]
    disproportionate=[p for p,x in results.items()
                      if x["original_rare_under5pct_visits_over10pct_pollen_taxa"]]
    return {
        "schema":SCHEMA,
        "original_source_Dryad_DOI":"10.5061/dryad.17pp3",
        "original_source_article_DOI":"10.1098/rspb.2015.1130",
        "original_source_exact_SHA256":"6994304ac377540c87b57560ecec50f72e8514494f9edd034eae2b5ed18f61c1",
        "original_visits_to_virgin_stigmas":1338,
        "original_independent_field_communities":1,
        "original_five_heath_plant_species":5,
        "main_alan_moisture_affinity_8p1_candidate_E_tetralix":True,
        "source_effectiveness_network":results,
        "original_plants_with_top_visitation_not_top_functional_contributor":mismatched,
        "original_plant_n_with_mismatch":len(mismatched),
        "plant_species_with_less5pct_visitor_but_over10pct_functional_delivery":disproportionate,
        "n_plants_rare_but_functional":len(disproportionate),
        "no_ALAN_treatment":True,
        "no_moth_visit_observed":True,
        "no_fertilization_or_seeds_measured":True,
        "do_not_count_n_plant_species_as_independent_field_replicates":True,
        "source_visitation_weighted_network_and_actual_pollen_network_not_equivalent":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    data=execute()
    dst=Path(a.out);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(data,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("ALAN_HEATH_FIVEPLANTS_FUNCTIONAL_POLLEN_NETWORK_RANK",json.dumps({
       "species_n":data["original_five_heath_plant_species"],
       "mismatch_n":data["original_plant_n_with_mismatch"],
       "mismatch_species":data["original_plants_with_top_visitation_not_top_functional_contributor"],
       "rare_function_n":data["n_plants_rare_but_functional"],
       "rare_function_species":data["plant_species_with_less5pct_visitor_but_over10pct_functional_delivery"],
       "species":{p:{
         "top_visitor":s["top_source_visitation_taxon"],
         "top_pollen":s["top_source_functional_pollen_depositor"],
         "top_visit_pollen_share":s["top_visitor_functional_pollen_fraction"],
         "top_pollen_contributor_pollen_share":s["top_functional_taxon_pollen_fraction"],
         "rare_important":s["original_rare_under5pct_visits_over10pct_pollen_taxa"],
         "fly_10percent_scenario":s["hypothetical_10pct_bee_visits_replaced_with_fly_visits"],
         "bootstrap_mismatch_frequency":
           s["source_day_cluster_bootstrap_stability"].get(
               "fraction_resamples_where_top_visitation_vs_pollen_taxon_differ")
         } for p,s in data["source_effectiveness_network"].items()}
    },sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
