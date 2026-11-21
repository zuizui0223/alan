#!/usr/bin/env python3
"""Study-preserving source-coded ALAN interaction contrast and observation audit.

This is a *descriptive source-codebook analysis*, not a random-effects
effect-size meta-analysis, nor an estimate of causal light effect on plants.
Original source: Seymoure et al. 2023 Figshare SHA-fixed xlsx, 81
interactions nested in 79 papers.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random

from scripts.analysis.alan_seymoure2023_true_workbook_coding_audit_20261119 import (
    SOURCE_SHA256, SOURCE_URL, SOURCE_FILE_ID, SOURCE_NAME, get, source_rows
)
from scripts.analysis.alan_seymoure2023_original_direction_study_units_20261119 import (
    _source_scientific_rows, study_key, summarize_source
)

SEED=20261119
N_BOOT=7999


def fisher_two_sided(a,b,c,d):
    """Exact fixed-margin 2x2 Fisher, probability-ordering two-sided."""
    nums=[int(z) for z in (a,b,c,d)]
    if min(nums)<0:
        raise ValueError("Negative table cell")
    a,b,c,d=nums
    N=a+b+c+d
    if not N:raise ValueError("Empty 2x2")
    R=a+b
    C=a+c
    def comb(n,k):return math.comb(n,k) if 0<=k<=n else 0
    denominator=comb(N,C)
    prob=lambda x:comb(R,x)*comb(N-R,C-x)/denominator
    pobs=prob(a)
    low=max(0,C-(N-R));high=min(R,C)
    return min(1.0,sum(prob(x) for x in range(low,high+1)
                      if prob(x)<=pobs*(1+1e-12)))


def paired_binomial_two_sided(a,b):
    """McNemar exact conditional on discordant coded source records."""
    if a<0 or b<0:raise ValueError("Negative counts")
    n=a+b
    if n==0:return 1.
    k=min(a,b)
    return min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n)


def quantiles(xs):
    xs=sorted(xs)
    if not xs:raise ValueError("No bootstrap replications")
    return [xs[math.floor((len(xs)-1)*0.025)],
            xs[math.ceil((len(xs)-1)*0.975)]]


def source_role_audit(original_rows,strict=True):
    rows,footer=_source_scientific_rows(original_rows,strict=strict)
    # Import source validation and original 79 DOI clusters.
    base=summarize_source(original_rows,strict=strict)
    papers=defaultdict(list)
    for r in rows:
        if r["interaction strength (IS) reduction"] not in ("0","1"):
            raise ValueError("Reduction source code changed")
        if r["IS enhancement"] not in ("0","1"):
            raise ValueError("Enhancement source code changed")
        if r["temporal mismatch induced"] not in ("0","1"):
            raise ValueError("Temporal code changed")
        if r["spatial mismatch induced"] not in ("0","1"):
            raise ValueError("Spatial code changed")
        papers[study_key(r)].append({
            "type":r["interaction_type"],
            "red":int(r["interaction strength (IS) reduction"]),
            "enh":int(r["IS enhancement"]),
            "tem":int(r["temporal mismatch induced"]),
            "spa":int(r["spatial mismatch induced"]),
            "method_is_recorded":bool(r.get("method","")),
            "year":r.get("Publication Year",""),
        })
    if strict and (len(papers)!=79 or len(rows)!=81):
        raise ValueError("Article / investigation grain drift")
    groups={}
    for name in ("predator-prey","herbivory","pollination","detritivory"):
        ps={k:[x for x in v if x["type"]==name] for k,v in papers.items()}
        ps={k:v for k,v in ps.items() if v}
        stage={
            "publication_n":len(ps),
            "investigation_n":sum(map(len,ps.values())),
            "source_binary_enhancement":sum(z["enh"] for v in ps.values() for z in v),
            "source_binary_reduction":sum(z["red"] for v in ps.values() for z in v),
            "source_mixed_sign_investigations":sum(z["red"] and z["enh"] for v in ps.values() for z in v),
            "source_method_recorded_investigations":sum(z["method_is_recorded"] for v in ps.values() for z in v)
        }
        groups[name]=stage
    predator={k:v for k,v in papers.items() if any(x["type"]=="predator-prey" for x in v)}
    herbivore={k:v for k,v in papers.items() if any(x["type"]=="herbivory" for x in v)}
    # Prevent non-independent same-publication comparison when type varies.
    overlapping=set(predator)&set(herbivore)
    if overlapping:
        raise ValueError("Cross-type comparisons share publications; requires joint cluster bootstrap")

    def select_paper_row(selected,typ):
        result=[]
        for k,v in selected.items():
            subset=[x for x in v if x["type"]==typ]
            if len(subset)!=1:
                raise ValueError("More than one same-type observation per publication")
            result.append(subset[0])
        return result
    pp=select_paper_row(predator,"predator-prey")
    hh=select_paper_row(herbivore,"herbivory")
    if strict and (len(pp)!=55 or len(hh)!=12):
        raise ValueError("Original guild denominators changed")

    def directional_record(v):
        inc=sum(x["enh"]==1 and x["red"]==0 for x in v)
        dec=sum(x["red"]==1 and x["enh"]==0 for x in v)
        both=sum(x["red"]==1 and x["enh"]==1 for x in v)
        neither=sum(x["red"]==0 and x["enh"]==0 for x in v)
        return {"enhancement_only":inc,"reduction_only":dec,
                "both":both,"neither":neither,"n":len(v),
                "source_signed_score_per_study":(inc-dec)/len(v)}
    pr=directional_record(pp); hr=directional_record(hh)
    odds=(pr["enhancement_only"]*hr["reduction_only"]/
          (pr["reduction_only"]*hr["enhancement_only"])
          if pr["reduction_only"]*hr["enhancement_only"] else None)
    p=fisher_two_sided(pr["enhancement_only"],pr["reduction_only"],
                       hr["enhancement_only"],hr["reduction_only"])
    expected_or=(19*7)/(13*2)
    if strict and (pr["enhancement_only"],pr["reduction_only"],pr["both"],
                   hr["enhancement_only"],hr["reduction_only"])!=(19,13,2,2,7):
        raise ValueError("Original exclusive direction 2x2 source changed")
    if strict and abs(odds-expected_or)>1e-10:
        raise ValueError("Exclusive directions OR drift")

    rng=random.Random(SEED)
    boots=[]
    for _ in range(N_BOOT):
        pboot=[pp[rng.randrange(len(pp))] for _ in pp]
        hboot=[hh[rng.randrange(len(hh))] for _ in hh]
        boots.append(directional_record(pboot)["source_signed_score_per_study"]-
                     directional_record(hboot)["source_signed_score_per_study"])

    allflags=[]
    for v in papers.values():
        # Preserve 79-paper source unit; combine flags only, never sum duplicates.
        allflags.append({"tem":int(any(z["tem"] for z in v)),
                         "spa":int(any(z["spa"] for z in v))})
    spatial_only=sum(x["spa"] and not x["tem"] for x in allflags)
    temporal_only=sum(x["tem"] and not x["spa"] for x in allflags)
    both=sum(x["tem"] and x["spa"] for x in allflags)
    neither=sum(not x["tem"] and not x["spa"] for x in allflags)
    original_s=base["publication_level_DOI_aggregated_counts"]
    if strict and not (spatial_only==24 and temporal_only==4 and both==3 and neither==48):
        raise ValueError("Review source time-space discordances changed")
    if strict and (original_s["coded_any_spatial"],original_s["coded_any_temporal"])!=(27,7):
        raise ValueError("Original paper-time and space coding changed")

    result={
        "schema":"alan_original_2023_seymoure_79paper_guild_time_space_descriptive_associations_v1",
        "source_SHA256":SOURCE_SHA256,
        "independent_publications_in_review":len(papers),
        "original_source_interactions":len(rows),
        "original_programme":"Seymoure et al. 2023 source review through 2022",
        "by_source_interaction_group":groups,
        "predator_prey_versus_herbivory_coding":{
            "predator_prey":pr,"herbivory":hr,
            "predator_prey_publication_overlap_herbivory":len(overlapping),
            "exclusive_enhancement_reduction_2x2_order":[
                [pr["enhancement_only"],pr["reduction_only"]],
                [hr["enhancement_only"],hr["reduction_only"]]],
            "odds_ratio_for_enhancement_vs_reduction_in_reported_exclusive_signs":odds,
            "exploratory_fisher_two_sided_fixed_margins_p":p,
            "all_coded_studies_signed_score_mean_difference":(
                pr["source_signed_score_per_study"]-
                hr["source_signed_score_per_study"]),
            "whole_publication_stratified_bootstrap_signed_score_difference_95pct_descriptive_interval":quantiles(boots),
            "bootstrap_n":N_BOOT,
            "bootstrap_seed":SEED,
            "p_is_exploratory_not_pre_registered":True,
            "not_estimate_of_impact_on_actual_species_interaction_strength":True,
            "interaction_type_methods_taxon_imbalanced":True
        },
        "within_publication_original_spatiotemporal_coding":{
            "spatial_only":spatial_only,
            "temporal_only":temporal_only,
            "both":both,
            "neither":neither,
            "spatial_total":spatial_only+both,
            "temporal_total":temporal_only+both,
            "source_coded_spatial_minus_temporal_proportion":(spatial_only-temporal_only)/len(papers),
            "exploratory_exact_McNemar_p":paired_binomial_two_sided(spatial_only,temporal_only),
            "not_a_test_of_biological_mechanism_size":True,
            "undetected_or_not_evaluated_not_equal_true_zero":True
        },
        "inference_gate":{
            "original_lighting_installation_counts_available":False,
            "original_effect_size_standard_errors_available":False,
            "source_method_missing_for_63_of_79_papers":True,
            "reviews_are_not_probabilistic_global_ecological_samples":True,
            "not_a_random_effects_effect_size_meta_analysis":True,
            "no_predator_mediation_of_pollen_or_seed_production_fitted":True,
            "multiple_hypotheses_unadjusted_exploratory_screen":True,
            "no_original_study_row_or_DOI_export":True
        }
    }
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    meta=json.loads(get(SOURCE_URL))
    xs=[f for f in meta.get("files",[]) if f.get("id")==SOURCE_FILE_ID and f.get("name")==SOURCE_NAME]
    if len(xs)!=1:raise ValueError("Original public Excel source missing")
    original=get(xs[0]["download_url"])
    if hashlib.sha256(original).hexdigest()!=SOURCE_SHA256:
        raise ValueError("Original source exact SHA mismatch")
    header, rows=source_rows(original)
    if len(header)!=20:raise ValueError("Unexpected original source header")
    ans=source_role_audit(rows,strict=True)
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(ans,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_SOURCE_REVIEW_GUILD_DIRECTION",ans["predator_prey_versus_herbivory_coding"])
    print("ALAN_SOURCE_REVIEW_TIME_SPACE",ans["within_publication_original_spatiotemporal_coding"])


if __name__=="__main__":
    main()
