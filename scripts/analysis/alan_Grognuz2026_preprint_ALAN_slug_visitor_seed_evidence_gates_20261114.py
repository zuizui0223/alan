#!/usr/bin/env python3
"""Auditable 2026 OPEN PREPRINT light×herbivory×daytime-visitor×seed path gate.

The published preprint (Grognuz & Knop, Mar 2026) describes an actual
2×2 ALAN × slug factorial, directly measured foliar/floral damage, daytime
insect visitor frequencies and seed set. Its raw Mendeley files were NOT
downloaded or re-fitted by this project. These source-*reported* coefficients
are not independent experimental replications and may have different
covariate scales/conditional interpretations.

This source audit answers a specific macro-mechanism question:
 Do individually reported significant ALAN→herbivore and herbivore→visitor
 effects form the SAME damage pathway and independently supported mediated
 net reproductive effect? NO: light changed folivory while the stronger
 herbivory×ALAN visitor interaction concerns florivory, which showed no
 detected direct light response. Neither damage×light seed interaction was
 significant in this original preprint.

Preserve all positive/non-positive source claims and never multiply or add
different-model betas as an estimated indirect effect.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

FILE="config/alan_public_Grognuz2026_herbivory_daypollination_seed_preprint_20261114.json"
SCHEMA="ALAN_public_2026_factorial_herbivore_light_pollination_seed_primary_manuscript_ledger_v1"

def require(cond,msg):
    if not bool(cond):raise ValueError("ALAN_2026_DIRECT_AND_HERBIVORY_PATH: "+str(msg))

def tested_sign(x):
    if "CI95" not in x:return "UNRESOLVED_OR_UNREPORTED"
    ci=x["CI95"]
    require(isinstance(ci,list) and len(ci)==2 and
            all(isinstance(a,(int,float)) and math.isfinite(a) for a in ci) and
            ci[0]<ci[1],"Original preprint source confidence interval malformed")
    beta=x.get("beta")
    require(beta is not None and ci[0]<beta<ci[1],
            "Original preprint summary beta outside reported confidence interval")
    if ci[0]>0:return "POSITIVE_REPORTED_CI_EXCLUDES_ZERO"
    if ci[1]<0:return "NEGATIVE_REPORTED_CI_EXCLUDES_ZERO"
    return "INCONCLUSIVE_REPORTED_CI_CONTAINS_ZERO"

def path(s):
    require(s["schema"]==SCHEMA,"Original source evidence contract schema mismatch")
    require(s["primary_source"]["status"]=="PUBLIC_PREPRINT_NOT_CONFIRMED_PEER_REVIEWED",
            "Do not promote unreviewed source to a published confirmed article")
    require(s["primary_source"]["raw_data_bytes_reacquired_in_this_ALAN_analysis"] is False,
            "Source raw Mendeley reanalysis has not been executed")
    d=s["design"];m=s["published_reported_coefficients"]
    require(d["plant_species"]==2 and d["four_treatment_combinations"] and
            d["light_exposure_nights"]==3 and
            d["total_visit_events_reported"]==819 and
            d["herbivory_treated_plant_source_count"]==276 and
            d["original_inflorescence_assessed_for_florivory"]==355,
            "Expected true preprint treatment design/samples changed")
    require(abs(m["ALAN_without_herbivores_day_pollinator_visits"]["beta"]+.38)<1e-9 and
            abs(m["ALAN_by_florivory_on_day_visitor_response"]["beta"]-4.6)<1e-9 and
            abs(m["ALAN_on_seed_set"]["beta"]+.40)<1e-9 and
            abs(m["ALAN_slug_damaged_leaf_odds_ratio"]["odds_ratio"]-.21)<1e-9,
            "Original preprint reported quantitative point values changed")
    arms={
      "ALAN_direct_day_visits":tested_sign(m["ALAN_without_herbivores_day_pollinator_visits"]),
      "ALAN_direct_Hymenoptera_day_visits":tested_sign(m["ALAN_day_hymenoptera_visits"]),
      "ALAN_direct_Diptera_day_visits":tested_sign(m["ALAN_day_diptera_visits"]),
      "ALAN_by_florivory_VISITS_interaction":tested_sign(m["ALAN_by_florivory_on_day_visitor_response"]),
      "ALAN_by_folivory_VISITS_interaction":tested_sign(m["ALAN_by_folivory_on_day_visitor_response"]),
      "florivory_on_VISITS":tested_sign(m["florivory_on_day_visits"]),
      "folivory_on_VISITS":tested_sign(m["folivory_on_day_visits"]),
      "ALAN_on_SEED_SET":tested_sign(m["ALAN_on_seed_set"])
    }
    leaf=m["ALAN_slug_damaged_leaf_odds_ratio"]
    require(leaf["CI95"][0]<leaf["odds_ratio"]<leaf["CI95"][1] and
            0<leaf["CI95"][0]<leaf["CI95"][1]<1,
            "Reported leaf herbivory odds cannot be called negative")
    require(arms["ALAN_by_florivory_VISITS_interaction"].startswith("POSITIVE") and
            arms["ALAN_by_folivory_VISITS_interaction"].startswith("INCONCLUSIVE") and
            arms["florivory_on_VISITS"].startswith("NEGATIVE") and
            arms["folivory_on_VISITS"].startswith("INCONCLUSIVE") and
            arms["ALAN_on_SEED_SET"].startswith("NEGATIVE"),
            "Original source sign gate must be documented honestly")
    require(m["ALAN_on_florivory_incidence"]["inferential_status"]==
            "NO_DETECTED_DIRECT_ALAN_TREATMENT_EFFECT_ON_FLORIVORY" and
            m["ALAN_by_herbivory_seed_set"]["inferential_status"]==
            "NOT_RESOLVED_STATISTICALLY_IN_ORIGINAL_PREPRINT",
            "The original preprint lacks a verified completed light→florivory→seed chain")
    leaf_chain={
      "original_published_light_to_LEAF_DAMAGE_effect":{
          "reported_OR":leaf["odds_ratio"],"reported_OR95CI":leaf["CI95"],
          "report_odds_reduction_percent":100*(1-leaf["odds_ratio"]),
          "not_a_direct_predator_kill_effect":True,
          "reported_CI_below_one":True
      },
      "original_folivory_to_day_visits_CI_status":arms["folivory_on_VISITS"],
      "original_ALAN_by_folivory_visitor_effect_CI_status":
            arms["ALAN_by_folivory_VISITS_interaction"],
      "complete_positive_pollination_mediation_established":False,
      "original_ALAN_by_herbivory_seed_set_confirmed":False
    }
    flower_chain={
      "original_light_to_FLOWER_DAMAGE_detected":False,
      "original_florivory_to_day_visits_CI_status":arms["florivory_on_VISITS"],
      "original_light_by_florivory_visitor_effect_CI_status":
          arms["ALAN_by_florivory_VISITS_interaction"],
      "complete_positive_pollination_mediation_established":False,
      "original_ALAN_by_herbivory_seed_set_confirmed":False
    }
    return {
      "schema":"ALAN_20261114_Grognuz2026_direct_and_herbivory_interaction_path_closure_v1",
      "status":"ACTUAL_PREPRINT_REPORTED_RESULT_PATH_AUDIT_NO_RAW_DATA_REANALYSIS",
      "source_version_and_inference_tier":{
         "preprint_DOI":s["primary_source"]["preprint_DOI"],
         "raw_dataset_DOI":s["primary_source"]["raw_data_DOI"],
         "independent_publications":1,
         "independent_regions":1,
         "species":2,
         "visits_original_preprint_reported":d["total_visit_events_reported"],
         "original_factorial_2x2":True,
         "raw_experiment_records_recomputed":False,
         "source_awaits_peer_review":True,
         "no_join_to_original_Giavi_source_2016_replicates":True
      },
      "original_published_effect_direction_CI_checks":arms,
      "source_foliar_damage_mediation_requirements":leaf_chain,
      "source_floral_damage_mediation_requirements":flower_chain,
      "seed_set":{
        "original_ALAN_seed_beta":m["ALAN_on_seed_set"]["beta"],
        "original_ALAN_seed_beta_CI95":m["ALAN_on_seed_set"]["CI95"],
        "original_ALAN_seed_p":m["ALAN_on_seed_set"]["p"],
        "reported_ALAN_by_herbivory_on_SEED_SET_status":
            m["ALAN_by_herbivory_seed_set"]["inferential_status"],
        "NO_actual_seed_assurance_buffering_by_herbivore_state_demonstrated":True
      },
      "source_quality_flags":{
        "different_damage_types_must_not_be_merged_in_a_single_mediator":True,
        "observed_interaction_does_not_prove_indirect_product_path":True,
        "primary_preprint_seed_effect_not_extrapolated_to_island_or_nocturnal_moth":True,
        "slug_feeding_reduced_not_slug_predator_pollination":True,
        "preprint_beta_scale_not_converted_to_percent":True,
        "raw_data_access_not_presumed":True,
        "no_invented_indirect_effect_standard_errors_or_p_values":True,
        "no_global_ALAN_pollen_service_meta_estimate":True
      }
    }

def main(repo,out):
    out=Path(out)
    require(not out.exists(),"Refuse overwrite prior source audit")
    s=json.loads((Path(repo)/FILE).read_text(encoding="utf-8"))
    result=path(s)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                   encoding="utf-8")
    print("ALAN_PREPRINT_HERBIVORY_LIGHT_VISIT_SEED_FULL_CAUSAL_CHAIN_GATES",
          json.dumps(result,ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo",type=Path,default=Path("."))
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    main(a.repo,a.out)
