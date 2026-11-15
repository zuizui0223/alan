#!/usr/bin/env python3
"""Source-locked Li et al (2026) ALAN reproduction stage accounting.

2026-07-10 preprint DOI 10.64898/2026.07.10.737511.
Real 79 author plant records in Dongboli/experimental-data git blob
b616aca5beeba1457549841477d61ca4f796a147. 3 source hayfield
blocks, 3 moth-associated species, C/P/L lighting treatment combinations.

Novel accounting distinguishes:
  (1) probability Seeds>0,
  (2) number of Seeds CONDITIONAL on (1),
  (3) fruit production / INITIAL flower-bud opportunity and
  (4) seeds per fruit CONDITIONAL on fruit production.
Arithmetic identities, not causal mediation or unobserved pollen service.
All original individual rows stay on the ephemeral CI runner, only anonymous
field-level aggregates in the artifact. Never impute missing plant observations.
"""
from __future__ import annotations
import argparse, math, json
from pathlib import Path
import numpy as np
from scripts.analysis.alan_Li2026_real_geometry_seed_hayfield_source_20261114 import (
    source_bytes,parse,require,GIT_BLOB_SHA1,TREATMENTS,SPECIES
)

def records(arm):
    known=[r for r in arm if r["pollinated"]!="NA"]
    if not known:return None
    require(all(r["Starting"]!="NA" for r in known),"initial flower buds absent")
    starting=sum(int(r["Starting"]) for r in known)
    require(starting>0,"expected original positive buds")
    pos=[r for r in known if int(r["pollinated"])>0]
    seeds=sum(int(r["Seeds"]) for r in known)
    fruits=sum(int(r["Fruit"]) for r in known)
    require(len(pos)>0 and fruits>0,"positive reproductive events required")
    require(all(int(r["Starting"])>=int(r["Fruit"]) for r in known),
            "source fruits exceed initial buds; exposure denominator invalid")
    require(all((int(r["Seeds"])>0)==(int(r["Fruit"])>0) for r in known),
            "fruit and seed presence disagree in original study")
    largest=max(known,key=lambda x:int(x["Seeds"]))
    n=len(known)
    return {
      "source_plant_rows":len(arm),"known_seed_outcome_plants":n,
      "missing_reproductive_outcomes":len(arm)-n,
      "seed_positive_plants":len(pos),
      "probability_any_seeds":len(pos)/n,
      "total_seeds_in_known_plants":seeds,
      "total_fruits_in_known_plants":fruits,
      "total_initial_buds_in_known_plants":starting,
      "mean_seeds_per_plant_all_known":seeds/n,
      "mean_fruits_per_plant_all_known":fruits/n,
      "mean_seeds_conditional_on_any_seeds":seeds/len(pos),
      "mean_seeds_per_initial_flower_bud":seeds/starting,
      "mean_fruits_per_initial_flower_bud":fruits/starting,
      "mean_seeds_per_fruit":seeds/fruits,
      "largest_seed_count_in_one_plant":int(largest["Seeds"]),
      "largest_plant_taxon_code":largest["species"],
      "hurdle_identity_absolute_error":
          abs((len(pos)/n)*(seeds/len(pos))-seeds/n),
      "fruit_identity_absolute_error":
          abs((fruits/starting)*(seeds/fruits)-seeds/starting),
      "effective_moth_visits_or_stigma_pollen_measured":False,
      "initial_buds_are_not_all_eventual_pollination_opportunities":True
    }

def ratio(x,y,k):
    require(x[k]>0 and y[k]>0,"positive source denominator needed: "+k)
    return x[k]/y[k]

def two_stage(a,b):
    """Exact identities on source observed means, non-mediation factorization."""
    probabilistic=ratio(a,b,"probability_any_seeds")
    conditional=ratio(a,b,"mean_seeds_conditional_on_any_seeds")
    perplant=ratio(a,b,"mean_seeds_per_plant_all_known")
    fruit=ratio(a,b,"mean_fruits_per_initial_flower_bud")
    seedperfruit=ratio(a,b,"mean_seeds_per_fruit")
    perbud=ratio(a,b,"mean_seeds_per_initial_flower_bud")
    require(math.isclose(probabilistic*conditional,perplant,rel_tol=1e-12),
            "source hurdle probability x conditional seed identity failed")
    require(math.isclose(fruit*seedperfruit,perbud,rel_tol=1e-12),
            "source fruits/bud x seeds/fruit identity failed")
    def logshare(z,whole):
        return math.log(z)/math.log(whole) if whole!=1 else None
    return {
      "seed_positive_probability_ratio":probabilistic,
      "conditional_seeds_given_positive_ratio":conditional,
      "observed_seeds_per_known_plant_ratio":perplant,
      "probability_and_conditional_seed_identity_error":probabilistic*conditional-perplant,
      "observed_fruits_per_initial_bud_ratio":fruit,
      "observed_seeds_per_fruit_ratio":seedperfruit,
      "observed_seeds_per_initial_bud_ratio":perbud,
      "fruit_and_seeds_per_fruit_identity_error":fruit*seedperfruit-perbud,
      "log_seed_per_bud_change_share_arithmetic_fruit_stage":logshare(fruit,perbud),
      "log_seed_per_bud_change_share_arithmetic_seed_per_fruit":logshare(seedperfruit,perbud),
      "noncausal_arithmetic_not_statistical_mediation":True
    }

def analyze(rows):
    known=[r for r in rows if r["pollinated"]!="NA"]
    require(len(rows)==79 and len(known)==78,"Original 79 author source rows and 78 known outcomes")
    full={t:records([r for r in rows if r["Treatment"]==t]) for t in ("C","P","L")}
    require([(t,full[t]["known_seed_outcome_plants"],full[t]["seed_positive_plants"])
             for t in ("C","P","L")]==[("C",26,8),("P",25,12),("L",27,20)],
            "Original source factorial arm sample changed")
    require([(t,full[t]["total_seeds_in_known_plants"],
              full[t]["total_fruits_in_known_plants"])
             for t in ("C","P","L")]==[("C",672,32),("P",967,33),("L",3729,149)],
            "Original seeds/fruits source source drift")
    require([(t,full[t]["total_initial_buds_in_known_plants"])
             for t in ("C","P","L")]==[("C",369),("P",408),("L",483)],
            "Original initial budding opportunity drift")
    fields={}
    for field in ("1","2","3"):
        q=[r for r in rows if r["replicate"]==field]
        fs={t:records([r for r in q if r["Treatment"]==t]) for t in ("C","P","L")}
        comp=two_stage(fs["L"],fs["C"])
        fields[field]={"original_plant_reproductive_stage_by_treatment":fs,
                      "linear_vs_dark_reproductive_stage":comp,
                      "linear_vs_dark_seed_presence_fraction_difference":
                         fs["L"]["probability_any_seeds"]-fs["C"]["probability_any_seeds"]}
    require(all(fields[f]["linear_vs_dark_reproductive_stage"]["observed_seeds_per_initial_bud_ratio"]>1
                for f in ("1","2","3")),
            "Positive L/C seed:initialbud real source not present in all 3 fields")
    require(all(fields[f]["linear_vs_dark_reproductive_stage"]["observed_fruits_per_initial_bud_ratio"]>1
                for f in ("1","2","3")),
            "Positive L/C fruit:initialbud real source not present in all 3 fields")

    # The L sample extreme plant is the highest seed yield of a true
    # originally observed author row, no duplicate-ID-imputation.
    array=[r for r in known if r["Treatment"]=="L"]
    maximum=max(array,key=lambda r:int(r["Seeds"]))
    without=[r for r in array if r is not maximum]
    robust=records(without)
    compare=two_stage(full["L"],full["C"])
    sample_outlier={
      "source_top_plant_seeds":int(maximum["Seeds"]),
      "source_top_plant_species_code":maximum["species"],
      "source_top_plant_original_field_label_anonymized":True,
      "top_plant_share_of_linear_total_seeds":
         int(maximum["Seeds"])/full["L"]["total_seeds_in_known_plants"],
      "original_linear_seeds_per_bud":full["L"]["mean_seeds_per_initial_flower_bud"],
      "drop_largest_linear_plant_seeds_per_bud":robust["mean_seeds_per_initial_flower_bud"],
      "drop_largest_linear_plant_seed_bud_ratio_to_unlit":
         robust["mean_seeds_per_initial_flower_bud"]/full["C"]["mean_seeds_per_initial_flower_bud"],
      "leave_one_field_out_ratios_seed_per_bud":{},
      "leave_one_field_out_ratios_fruits_per_bud":{}
    }
    for excluded in ("1","2","3"):
        a={t:records([r for r in known if r["replicate"]!=excluded and r["Treatment"]==t])
           for t in ("C","L")}
        c=two_stage(a["L"],a["C"])
        sample_outlier["leave_one_field_out_ratios_seed_per_bud"][excluded]=c["observed_seeds_per_initial_bud_ratio"]
        sample_outlier["leave_one_field_out_ratios_fruits_per_bud"][excluded]=c["observed_fruits_per_initial_bud_ratio"]
    # Multi-stage seed output source was not an hour-matched observation of
    # moth behavior or a true fruit-set causal mediating effect.
    return {
      "schema":"ALAN_Li2026_authentic_original_three_field_reproductive_hurdle_fruit_stage_v1",
      "original_author_data_git_blob_SHA1":GIT_BLOB_SHA1,
      "article_preprint_DOI":"10.64898/2026.07.10.737511",
      "original_experiment":{"independent_hayfields":3,"plant_species":list(SPECIES.values()),
       "lighting_treatments":TREATMENTS,
       "actual_source_plant_rows":len(rows),"observed_outcome_plants":len(known),
       "one_point_light_plant_seed_outcome_missing_not_imputed":True,
       "published_protocol_possible_deployment_slots":90,
       "no_evidence_11_potential_plants_were_actually_measured":True},
      "author_source_arm_reproductive_stages":full,
      "three_true_hayfield_reproductive_stages":fields,
      "linear_vs_dark_observed_composite":compare,
      "point_vs_dark_observed_composite":two_stage(full["P"],full["C"]),
      "largest_original_plant_and_leave_one_field_out_sensitivity":sample_outlier,
      "strict_biology_boundaries":{
         "seed_presence_equals_actual_measured_stigma_pollen":False,
         "original_nocturnal_visitors_or_pollen_transfer_measured":False,
         "same_light_brightness_for_point_and_array":False,
         "seed_count_is_conditional_on_recorded_experimental_plants":True,
         "seeds_per_initial_bud_not_seed_set_probability_per_flower":True,
         "positive_seed_samples_selected_postoutcome_for_conditional_mean":True,
         "true_independent_light_randomizations_not_78_plants":True,
         "published_preprint_not_peer_review_certified":True,
         "ecological_compensation_and_moth_corridor_unproven":True,
         "global_ALAN_moth_syndrome_coefficient_not_claimed":True}
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    opt=ap.parse_args()
    out=opt.out
    if out.exists():raise FileExistsError("No overwriting checked author data receipts")
    rows=parse(source_bytes())
    result=analyze(rows)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("LI2026_TRUE_SEED_REPRO_HURDLE_SOURCE",json.dumps({
      "arms":{k:{fld:v[fld] for fld in ("known_seed_outcome_plants","seed_positive_plants",
            "mean_seeds_conditional_on_any_seeds","mean_fruits_per_initial_flower_bud",
            "mean_seeds_per_fruit","mean_seeds_per_initial_flower_bud")}
            for k,v in result["author_source_arm_reproductive_stages"].items()},
      "linear_over_dark":result["linear_vs_dark_observed_composite"],
      "three_field_ratio":{k:z["linear_vs_dark_reproductive_stage"]["observed_seeds_per_initial_bud_ratio"]
          for k,z in result["three_true_hayfield_reproductive_stages"].items()},
      "outlier":result["largest_original_plant_and_leave_one_field_out_sensitivity"]},ensure_ascii=False))
if __name__=="__main__":main()
