#!/usr/bin/env python3
"""ALAN 3-stage systematic-review SEED readiness census and PRISMA guard.

The records were chosen from prior analyses, NOT obtained from independent
WoS/Scopus searches. Their counts must never populate PRISMA study inclusion
or produce rma.mv results. Keep review type, original paper, physical lighting
programme and endpoint stage distinct. Source conditions are tested now,
without manufacturing a 'k' of eligible quantitative primary light contrasts.
"""
import argparse
import json
from collections import Counter
from pathlib import Path

SOURCE=Path("config/alan_prisma_three_stage_seed_census_20261119.json")
SCHEMA="alan_prisma_three_stage_source_seed_v1"
ALLOW_STAGES={"S1","S2"}
REVIEW_DOIS={"10.1038/s41559-020-01322-x","10.1002/ece3.9990","10.1098/rstb.2022.0356"}

def audit(inp):
    if inp.get("schema")!=SCHEMA:
        raise ValueError("Three-stage source seed schema drift")
    if inp.get("review_status")!="DRAFT_PROTOCOL_NOT_PROSPECTIVE_PRE_REGISTRATION_AND_NOT_SYSTEMATICALLY_SEARCHED":
        raise ValueError("Prior literature discoveries cannot be called newly preregistered")
    if not inp.get("no_PRISMA_record_counts_claimed") or not inp.get("effect_observed_pre_protocol"):
        raise ValueError("Post-hoc existing results need an explicit PRISMA and temporal disclosure")
    if not inp.get("programme_vs_paper_vs_effect_vs_field_unit_distinguished"):
        raise ValueError("Original unit of independent ALAN illumination changed")
    baseline=inp["baseline_reviews"]
    if len(baseline)!=3 or {x["doi"] for x in baseline}!=REVIEW_DOIS:
        raise ValueError("Three pre-existing reviews are not original research experiments")
    published=next(r for r in baseline if r["doi"]=="10.1038/s41559-020-01322-x")
    if published["published_126_papers"]!=126 or published["published_1304_effects"]!=1304:
        raise ValueError("Original Sanders source size changed")
    if published["source_metric"]!="Hedges d; NOT automatically lnRR":
        raise ValueError("Sanders' original standardized effects are not automatically lnRR")
    original_g=next(r for r in baseline if r["doi"]=="10.1002/ece3.9990")
    if (original_g["source_light_subpanel_articles"]!=22 or
        original_g["source_UVB_light_articles"]!=8 or
        original_g["source_ALAN_nonUVB_light_candidate_articles"]!=14):
        raise ValueError("Published 2023 review UV-B contamination audit drift")
    original_s=next(r for r in baseline if r["doi"]=="10.1098/rstb.2022.0356")
    if original_s["published_papers"]!=79 or original_s["published_interaction_investigations"]!=81:
        raise ValueError("Original interaction review 81 coded endpoints are not 81 independent ALAN effects")
    overlap=inp["source_reference_id_concordance"]
    if not all(overlap.get(k) is True for k in
       ("penryn_2015_2018_2022_same_54_troughs",
        "lizard_2017_2018_same_30_light_plots",
        "Lu2026_JPE_and_EE_overlap_unverified",
        "westhavelland_2017_2025_same_1lit_1dark_sitepair",
        "Seymoure_review_and_Guenat_review_original_papers_must_deduplicate")):
        raise ValueError("Physical source independence assumptions missing")
    rows=inp.get("seed_candidates",[])
    if any(r.get("stage") not in ALLOW_STAGES for r in rows):
        raise ValueError("Only two light-outcome literature branches are first-screen source seeds")
    if len({(r["stage"],r["doi"]) for r in rows})!=len(rows):
        raise ValueError("Duplicate stage×DOI in manually seeded original study list")
    if len({(r["stage"],r["study_programme_id"]) for r in rows})!=len(rows):
        raise ValueError("Nonunique original programme key within stage")
    if any(r["doi"] in REVIEW_DOIS for r in rows):
        raise ValueError("A review paper cannot become a primary lighting experiment")
    for r in rows:
        if r["not_PRISMA_included_yet"] is not True or r["screening_status"]!="SEED_REQUIRES_FULLTEXT_AND_EFFECT_VARIANCE":
            raise ValueError("Known prior literature is not systematically included data")
        if not r.get("independence_note") or not r.get("endpoint_candidates") or not r.get("design"):
            raise ValueError("Primary study source or endpoint not identified")
    counts=Counter(r["stage"] for r in rows)
    if counts!={"S1":11,"S2":11}:
        raise ValueError("Current evidence map 11 A/11 B post-hoc source paper seeds changed")
    if (inp["counts"].get("S1_initial_seed_at_protocol_creation_20261119")!=9 or
        inp["counts"].get("S1_added_after_initial_seed_before_PRISMA_search")!=2):
        raise ValueError("Original nine-paper baseline plus two later known primary papers not preserved")
    if (inp["counts"]["S1_purposeful_candidate_paper_records_not_search_inclusions"]!=counts["S1"] or
        inp["counts"]["S2_purposeful_candidate_paper_records_not_search_inclusions"]!=counts["S2"]):
        raise ValueError("Seed counts must reflect known candidates, not PRISMA inclusions")
    reported=inp["counts"]
    for key in ("S1_confirmed_readily_extractable_lnRR_with_correct_variance",
                "S2_confirmed_readily_extractable_lnRR_with_correct_variance"):
        if reported[key]!="UNKNOWN_NEEDS_ENDPOINT_RECHECK":
            raise ValueError("Meta-eligible independent original effects remain unassessed")
    if not any(r["study_programme_id"]=="Manfrin2025" and r["habitat"]=="aquatic_riparian" for r in rows):
        raise ValueError("2025 wetland-subsidy counterexample excluded")
    west=[r for r in rows if r.get("physical_exposure_programme_id")=="WesthavellandDrainageDitch2012"]
    if len(west)!=2 or {r["doi"] for r in west}!={"10.3389/fenvs.2017.00061","10.1007/s00027-025-01161-7"}:
        raise ValueError("Two riparian papers must share one original Westhavelland field-site lighting programme")
    if not all("one" in r["independence_note"].lower() or "ONE" in r["independence_note"] for r in west):
        raise ValueError("Never mistake within-site streetlights or traps for independent experiment sites")
    if not any(r["study_programme_id"]=="Macgregor2019" for r in rows):
        raise ValueError("Known positive reproduction counterexample excluded")
    if not any(r["study_programme_id"]=="Farnworth2026" for r in rows):
        raise ValueError("Known forest opposite-signed predator counterexample excluded")
    return {
      "schema":"alan_20261119_three_stage_PRISMA_source_seed_readiness_guard_v1",
      "source_status":"POST_HOC_INITIAL_SEED_CENSUS_NOT_A_SYSTEMATIC_SEARCH",
      "S1_known_intentionally_seeded_primary_paper_records":counts["S1"],
      "S2_known_intentionally_seeded_primary_paper_records":counts["S2"],
      "published_review_counts_for_source_discovery_only":{
        "Sanders_126_papers_effects":1304,
        "Guenat_nightlight_nonUVB_candidates":14,
        "Seymoure_original_papers":79,
        "Seymoure_original_interaction_records":81,
      },
      "original_source_2018_2022_54trough_three_papers_counted_as_one_lighting_programme":True,
      "original_source_2017_2018_30plot_two_papers_counted_as_one_lighting_programme":True,
      "two_Lu2026_publications_original_plot_overlap_unverified":True,
      "Westhavelland_2017_2025_share_one_lit_one_dark_field_pair":True,
      "source_A_11_seed_papers_not_11_independent_light_programmes":True,
      "S1_primary_META_ready_k":None,
      "S2_primary_META_ready_k":None,
      "Stage3_plant_specific_conspecific_pollen_compensation_ready_k":None,
      "WoS_Scopus_original_searches_executed":False,
      "PRISMA_records_identified_from_database":None,
      "PRISMA_records_title_screened":None,
      "PRISMA_fulltext_screened":None,
      "PRISMA_S1_included":None,
      "PRISMA_S2_included":None,
      "PRISMA_S3_singlevisit_efficiency_included":None,
      "zero_is_not_used_for_unknown_PRISMA_counts":True,
      "no_lnRR_meta_or_causal_mediation_coefficient_estimated":True,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,default=SOURCE)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    source=json.loads(args.source.read_text())
    out=audit(source)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_PRISMA_STAGE1_KNOWN_NON_SYSTEMATIC_SEED",out["S1_known_intentionally_seeded_primary_paper_records"])
    print("ALAN_PRISMA_STAGE2_KNOWN_NON_SYSTEMATIC_SEED",out["S2_known_intentionally_seeded_primary_paper_records"])
    print("ALAN_PRISMA_UNKNOWN_ELIGIBLE_AND_SEARCH_COUNTS",out["S1_primary_META_ready_k"],out["S2_primary_META_ready_k"],out["PRISMA_records_identified_from_database"])

if __name__=="__main__":
    main()
