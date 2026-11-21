#!/usr/bin/env python3
"""ALAN 2023 review: exact, source-authentic 17-primary-paper ENDPOINT eligibility.

Do not pool effect sizes: "herbivory" in the review includes foraging animals,
marine biofilm snails, and simulated nectary experiments. The actual plant
reproductive outcomes are from only three original publications in this cohort.
"""
import argparse
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path

from scripts.analysis.alan_seymoure2023_true_workbook_coding_audit_20261119 import (
    SOURCE_SHA256,SOURCE_URL,SOURCE_FILE_ID,SOURCE_NAME,get,source_rows
)
from scripts.acquisition.alan_seymoure2023_plant17_primary_bibliography_20261119 import source_bibliography

DEFAULT_REGISTRY=Path("results/all_species_forest/biotic_interactions/seymoure2023_plant17_primary_endpoint_tiers_20261119.json")
SCHEMA="alan_seymoure2023_source_17plant_guild_primary_endpoint_tiers_20261119_v1"
STAGES={
"flowering_plant_reproduction":3,
"plant_resource_and_herbivore_demography":2,
"vertebrate_artificial_food_foraging":4,
"moth_body_pollen_transport":1,
"caterpillar_biomass_proxy":1,
"gastropod_surface_food_collecting":2,
"marine_microbiofilm_grazing":1,
"moth_community_abundance_only":1,
"modeled_sensory_vision_only":1,
"hawkmoth_artificial_flower_learning":1,
}

def classify(registry,publications,strict=True):
    if registry.get("schema")!=SCHEMA:raise ValueError("Primary endpoint scheme changed")
    if registry.get("source_original_review_SHA256")!=SOURCE_SHA256:raise ValueError("Review provenance drift")
    rows=registry.get("rows",[])
    if len(rows)!=len(publications):raise ValueError("Full source bibliography not covered")
    source={r["doi"]:r for r in publications}
    classified={r["doi"]:r for r in rows}
    if len(source)!=len(publications) or len(classified)!=len(rows):raise ValueError("DOI double counted")
    if set(source)!=set(classified):raise ValueError("Missing or unverified original source DOI")
    nsource=sum(r["original_2023_review_source_interaction_count"] for r in publications)
    if strict and (len(source)!=17 or nsource!=18):raise ValueError("Original 17 papers/18 interaction membership changed")
    counts=Counter()
    byreview=defaultdict(Counter)
    publication_stage=[]
    actual_outcomes=[]
    for doi,src in sorted(source.items()):
        r=classified[doi]
        stage=r["stage"]
        if stage not in STAGES:raise ValueError("Unrecognized primary biological endpoint class")
        if not r["primary_article"].startswith("https://"):raise ValueError("Missing public primary source")
        if not r.get("actually_measured_endpoint") or not r.get("evidence_limit"):
            raise ValueError("Missing measured endpoint or source design limitation")
        if r["measured_living_flowering_plant_reproductive_fitness"] != (stage=="flowering_plant_reproduction"):
            raise ValueError("Plant seed/fruit endpoint mismapped")
        if r["measured_direct_plant_resource_growth"] != (stage=="plant_resource_and_herbivore_demography"):
            raise ValueError("Plant growth and reproduction incorrectly pooled")
        if r.get("has_verified_independent_light_effect_SE_for_same_comparable_plant_seed_or_fruit_endpoint") is not False:
            raise ValueError("No cross-paper common randomized plant effect SE is source verified")
        types=src["original_2023_review_source_types"]
        if not types or any(t not in ("pollination","herbivory") for t in types):
            raise ValueError("Review interaction type drift")
        counts[stage]+=1
        for t in types:byreview[t][stage]+=1
        publication_stage.append({"DOI":doi,"review_source_types":types,
                                  "actual_measured_endpoint_stage":stage,
                                  "primary_paper_source":r["primary_article"],
                                  "study_design_as_primary_article_described":r["source_design"],
                                  "original_primary_outcome_and_limit":r["actually_measured_endpoint"],
                                  "source_uncertainty":r["evidence_limit"],
                                  "plant_reproductive_endpoint_measured":r["measured_living_flowering_plant_reproductive_fitness"]})
        if r["measured_living_flowering_plant_reproductive_fitness"]:
            actual_outcomes.append({"doi":doi,"endpoint":r["actually_measured_endpoint"],"source":r["primary_article"]})
    if strict and dict(counts)!=STAGES:raise ValueError("Primary endpoint tiers require evidence audit; class count drift")
    if strict and sum(byreview["pollination"].values())!=6:raise ValueError("Review pollination denominator drift")
    if strict and sum(byreview["herbivory"].values())!=12:raise ValueError("Review herbivory denominator drift")
    # Source "pollination" should not be interpreted as 6 measured plant seed/fruit studies.
    if strict and byreview["pollination"]["flowering_plant_reproduction"]!=3:
        raise ValueError("Original source pollination stage count drift")
    # The shared plant reproduction study also was coded in review herbivory.
    if strict and byreview["herbivory"]["flowering_plant_reproduction"]!=1:
        raise ValueError("Dual-coded original study missing")
    return {
       "schema":"alan_source_17_plant_primary_endpoint_stage_eligibility_results_v1",
       "genuine_Seymoure2023_original_source_SHA256":SOURCE_SHA256,
       "source_review_paper_n":79,"source_review_interaction_n":81,
       "priority_original_paper_n":len(source),
       "priority_original_review_interaction_n":nsource,
       "mutually_exclusive_primary_measured_endpoints_by_paper":dict(sorted(counts.items())),
       "review_guild_vs_real_primary_outcome_matrix":{t:dict(sorted(v.items())) for t,v in sorted(byreview.items())},
       "actual_flowering_plant_fruit_seed_papers_n":len(actual_outcomes),
       "actual_flowering_plant_fruit_seed_primary_papers":actual_outcomes,
       "plant_resources_growth_without_reproductive_fitness_papers_n":counts["plant_resource_and_herbivore_demography"],
       "other_stage_papers_n":len(source)-len(actual_outcomes)-counts["plant_resource_and_herbivore_demography"],
       "source_published_paper_stage_identifiers":publication_stage,
       "effect_size_meta_gate":{
           "plant_reproductive_primary_papers_eligible_for_paper_review":len(actual_outcomes),
           "verified_common_same_endpoint_effects_with_original_SE":0,
           "reproductive_endpoint_one_kind_and_species_across_all_three":False,
           "source_study_lighting_independent_clusters_verified_by_2023_review":False,
           "one_source_paper_counts_as_two_independent_original_studies":False,
           "no_meta_pooled_mean_or_confidence_interval_fitted":True,
           "plant_vs_microbe_or_bait_foraging_cannot_be_pooled":True,
           "observational_skyglow_gradient_not_randomized_lighting":True,
       }
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--registry",type=Path,default=DEFAULT_REGISTRY)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    registry=json.loads(args.registry.read_text())
    meta=json.loads(get(SOURCE_URL))
    x=[v for v in meta.get("files",[]) if v.get("id")==SOURCE_FILE_ID and v.get("name")==SOURCE_NAME]
    if len(x)!=1:raise ValueError("Source file identity changed")
    b=get(x[0]["download_url"])
    if hashlib.sha256(b).hexdigest()!=SOURCE_SHA256:raise ValueError("Original source hash failed")
    hdr,rows=source_rows(b)
    if len(hdr)!=20:raise ValueError("Original source column set changed")
    bibliography=source_bibliography(rows,strict=True)
    result=classify(registry,bibliography,strict=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_TRUE_PRIMARY_ENDPOINT_STAGE_COUNTS",
          result["mutually_exclusive_primary_measured_endpoints_by_paper"])
    print("ALAN_REVIEW_POLLINATION_HERBIVORY_VS_ACTUAL_MEASUREMENT",
          result["review_guild_vs_real_primary_outcome_matrix"])
    print("ALAN_PRIMARY_STUDY_COMMON_FITNESS_META_GATE",
          result["effect_size_meta_gate"])

if __name__=="__main__":
    main()
