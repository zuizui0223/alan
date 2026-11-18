#!/usr/bin/env python3
"""Source-authentic ALAN biotic mechanism: hidden pollination-function bottlenecks.

Compare visitation rates v_j versus pollen delivery proxies p_j =
(mean conspecific stigma pollen grains per visit q_j) * v_j across each
original Page et al. 2021 study×plant×survey-period×source metric cell.

Original Figshare #16160265 623 source taxon entries, pinned SHA256:
e28ebe99e0dfb27b73fbe6638e568bf98fa95b92c920e1d3d981e8eaea27519a

Population frozen by source-gates in the previously audited repo:
72 comparable plant-network cells, 33 source studies, 62 named plants.
This is not 72 independent lighting interventions, a causal extinction
model, a complete effective pollination network, or a direct ALAN trial.

Primary descriptive diagnostic:
D_i = 0.5 * sum_j |V_ij / sum_k V_ik - (q_ij V_ij) / sum_k (q_ik V_ik)|
     = total variation (TV) between the visitor-frequency and the
       conspecific-pollen-delivery-weighted source networks.

This is an observational mismatch of compositional network weights.
It cannot say that pollination q remains fixed after a species loss,
that a predator specifically intercepts the effective pollinator,
or that an ALAN treatment caused that loss. A post-removal q*V proxy
loss is a hypothetical frozen-network vulnerability only.

Supplementary source-matched diagnostics:
  A: pollen effective-taxa Hill2 versus visitation effective taxa Hill2
  B: top frequent taxon's |pollen share - visitor share|
  C: how much pollen importance is missed when conserving just the
     most frequent versus the actually top pollen supplier
  D: signs/magnitude stratified by 1 vs ≥2 original visitor guilds,
     2, 3+, and 5+ taxa, and no recorded q==0.
All multi-plant source studies are block-resampled at study, not taxon,
to form descriptive percentile bootstrap intervals. All endpoints were
defined here before first data-source rerun and reported including nulls.

No direct moth↔fly q on 58 selected wet ALAN plant taxa is produced.
"""
from __future__ import annotations
import argparse, collections, json, math, random, statistics
from pathlib import Path
from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import (
    get_source, num, VKEY
)
from scripts.analysis.alan_page2021_original_macro_functional_network_mismatch import (
    candidate_groups
)
SCHEMA="alan_page2021_33study_72plant_network_functional_tv_and_loss_v1"
SHA="e28ebe99e0dfb27b73fbe6638e568bf98fa95b92c920e1d3d981e8eaea27519a"
SEED=20261118
NBOOT=4999

def checked(ok,msg):
    if not bool(ok):raise ValueError("ALAN_BIOTIC_FUNCTIONAL_TV_SOURCE_STOP: "+str(msg))

def measures(source_taxa):
    """Pollinator taxon tuples of (visits, mean conspecific pollen) in same plant."""
    checked(len(source_taxa)>=2,"requires at least two observed taxa")
    V=[float(t[0]) for t in source_taxa]
    Q=[float(t[1]) for t in source_taxa]
    checked(all(math.isfinite(v) and v>0 for v in V),"invalid source visitor rate")
    checked(all(math.isfinite(q) and q>=0 for q in Q),"invalid source pollen mean")
    p=[v*q for v,q in zip(V,Q)]
    checked(sum(V)>0 and sum(p)>0,"source zero pollen flux")
    vs=[v/sum(V) for v in V]
    ps=[x/sum(p) for x in p]
    tv=.5*sum(abs(v-f) for v,f in zip(vs,ps))
    checked(-1e-9<=tv<=1+1e-9,"invalid compositional total variation")
    visit_top=max(range(len(V)),key=lambda i:V[i])
    func_top=max(range(len(V)),key=lambda i:p[i])
    n2=lambda proportions:1/sum(t*t for t in proportions)
    return {
        "source_taxon_visitation_composition":vs,
        "source_taxon_pollen_delivery_composition":ps,
        "functional_compositional_total_variation":tv,
        "visitation_effective_taxa_Hill2":n2(vs),
        "pollen_delivery_effective_taxa_Hill2":n2(ps),
        "top_visitor_taxon_index":visit_top,
        "top_pollen_taxon_index":func_top,
        "top_visitor_visitation_share":vs[visit_top],
        "top_visitor_true_pollen_share":ps[visit_top],
        "top_pollen_taxon_true_pollen_share":ps[func_top],
        "top_visit_proxy_absolute_pollen_share_error":
            abs(ps[visit_top]-vs[visit_top]),
        "top_pollen_contributor_extra_share_over_top_visitor":
            ps[func_top]-ps[visit_top],
        "functional_network_narrower":n2(ps)<n2(vs),
        "source_zero_pollen_means":sum(q==0 for q in Q),
        "all_outcomes_are_observational_not_ALAN_experimental":True,
    }

def audit_networks(raw,taxon_revision=False):
    accepted,excluded,source_cells=candidate_groups(raw,taxon_revision)
    checked(len(accepted)==72 and len({x["source_study"] for x in accepted})==33,
            "original 72 plant network /33 source studies eligibility changed")
    original=collections.defaultdict(list)
    for r in raw:
        if taxon_revision and r["pollinator.group"].strip()=="moth" and (
            r["pollinator.taxon"].strip()=="Vanessa carye"):
            continue
        original[tuple(r[k].strip() for k in VKEY)].append(r)
    audited=[]
    for a in accepted:
        key=(a["source_study"],a["source_plant"],a["source_study_year_period"],
             a["original_SVE_measure"],"mean",a["original_visitation_rate_unit"])
        data=original.get(key,[])
        d=[]
        guilds=set()
        for r in data:
            q=num(r["effectiveness.average"])
            rate=num(r["visit.value"])
            if q is not None and rate is not None and q>=0 and rate>0:
                d.append((rate,q))
                guilds.add(r["pollinator.group"].strip())
        checked(len(d)==a["original_source_pollinator_taxa"],
                "source matching does not reproduce source-pinned taxa count")
        stats=measures(d)
        checked(abs(stats["visitation_effective_taxa_Hill2"]-
                    a["Hill2_effective_taxon_number_by_visit_rate"])<1e-9,
                "expected visit Hill2 source mismatch")
        checked(abs(stats["pollen_delivery_effective_taxa_Hill2"]-
                    a["Hill2_effective_taxon_number_by_pollen_delivery"])<1e-9,
                "expected pollen Hill2 source mismatch")
        checked(abs(stats["top_visit_proxy_absolute_pollen_share_error"]-
                    abs(a["top_visitor_functional_pollen_share"]-
                        a["top_visitor_visit_rate_share"]))<1e-9,
                "expected dominant-visitor proxy source mismatch")
        audited.append({
            "source_study":a["source_study"],
            "source_plant":a["source_plant"],
            "source_year_period":a["source_study_year_period"],
            "original_visitor_taxa":len(d),
            "original_visitor_guilds":len(guilds),
            "source_zero_pollen_means":stats["source_zero_pollen_means"],
            "functional_total_variation":stats["functional_compositional_total_variation"],
            "visitor_Hill2":stats["visitation_effective_taxa_Hill2"],
            "functional_Hill2":stats["pollen_delivery_effective_taxa_Hill2"],
            "functional_Hill2_ratio_to_visitor_Hill2":
                 stats["pollen_delivery_effective_taxa_Hill2"]/
                 stats["visitation_effective_taxa_Hill2"],
            "functionally_narrower":stats["functional_network_narrower"],
            "frequent_visitor_proxy_abs_error":
                 stats["top_visit_proxy_absolute_pollen_share_error"],
            "frequent_visitor_proxy_signed_error":
                 stats["top_visitor_true_pollen_share"]-
                 stats["top_visitor_visitation_share"],
            "functional_top_taxon_pollen_share_gap":
                 stats["top_pollen_contributor_extra_share_over_top_visitor"],
            "most_visited_differs_most_functional":
                 not a["top_visitor_equals_top_functional_pollen"],
            "original_rare_functional_taxon_present":a["rare_functional_taxa_n"]>0,
            "original_primary_endpoint_is_pollen_not_seed_set":True
        })
    return audited,excluded,source_cells

def smedian(x):
    return statistics.median(x) if x else None

def average(x):
    return statistics.fmean(x) if x else None

def aggregate(records):
    if not records:return {"status":"source_support_0_STOP"}
    g=lambda key:[x[key] for x in records]
    bystudy=collections.defaultdict(list)
    for r in records:bystudy[r["source_study"]].append(r)
    stats={
      "status":"source_descriptive",
      "source_plant_network_cells":len(records),
      "source_studies":len(bystudy),
      "distinct_source_plant_labels":len({r["source_plant"] for r in records}),
      "functional_weight_vs_visitation_weight_TV_mean":average(g("functional_total_variation")),
      "functional_weight_vs_visitation_weight_TV_median":smedian(g("functional_total_variation")),
      "TV_fraction_ge0p2":sum(x>=.2 for x in g("functional_total_variation"))/len(records),
      "functional_Hill2_mean":average(g("functional_Hill2")),
      "visitor_Hill2_mean":average(g("visitor_Hill2")),
      "within_plant_Hill2_difference_function_minus_visitor_mean":average([
          r["functional_Hill2"]-r["visitor_Hill2"] for r in records]),
      "fraction_functional_Hill2_less_than_visitor_Hill2":average([
          float(x) for x in g("functionally_narrower")]),
      "frequent_visitor_proxy_abs_pollen_share_error_mean":
          average(g("frequent_visitor_proxy_abs_error")),
      "frequent_visitor_proxy_abs_error_fraction_ge0p2":
          sum(x>=.2 for x in g("frequent_visitor_proxy_abs_error"))/len(records),
      "top_functional_contributor_over_top_visitor_pollen_share_gap_mean":
          average(g("functional_top_taxon_pollen_share_gap")),
      "top_functional_share_gap_among_rank_mismatch_only":
          average([r["functional_top_taxon_pollen_share_gap"] for r in records
                   if r["most_visited_differs_most_functional"]]),
      "rank_mismatch_fraction":average([float(x) for x in g(
          "most_visited_differs_most_functional")]),
      "functional_rare_taxon_group_fraction":average([float(x)
          for x in g("original_rare_functional_taxon_present")]),
      "source_study_balanced_TV_mean":average([
          average([r["functional_total_variation"] for r in group])
          for group in bystudy.values()]),
      "source_study_balanced_functionally_narrow_fraction":average([
          average([float(r["functionally_narrower"]) for r in group])
          for group in bystudy.values()]),
      "original_source_study_may_contain_multiple_plant_networks":True
    }
    return stats

METRICS=(
 "functional_weight_vs_visitation_weight_TV_mean",
 "fraction_functional_Hill2_less_than_visitor_Hill2",
 "within_plant_Hill2_difference_function_minus_visitor_mean",
 "frequent_visitor_proxy_abs_pollen_share_error_mean",
 "top_functional_contributor_over_top_visitor_pollen_share_gap_mean",
 "rank_mismatch_fraction"
)
def study_bootstrap(records,reps=NBOOT):
    studies=collections.defaultdict(list)
    for r in records:studies[r["source_study"]].append(r)
    ids=list(studies)
    if len(ids)<5:return {"status":"fewer_than_five_source_studies_STOP",
                           "source_studies":len(ids)}
    rng=random.Random(SEED)
    draws={key:[] for key in METRICS}
    for _ in range(reps):
        batch=[r for sid in rng.choices(ids,k=len(ids)) for r in studies[sid]]
        q=aggregate(batch)
        for key in METRICS:draws[key].append(q[key])
    intervals={}
    for key,vals in draws.items():
        vals.sort()
        intervals[key]=[vals[int(.025*(reps-1))],vals[int(.975*(reps-1))]]
    return {
      "status":"descriptive_bootstrap_resampling_original_study_ID",
      "original_independent_source_studies":len(ids),
      "bootstrap_replicates":reps,
      "nominal_descriptive_95intervals":intervals,
      "not_an_ALAN_or_pollinator_loss_randomized_test":True
    }

def execute(raw):
    versions={}
    for label,rev in (("as_published",False),
                      ("excluding_verified_Vanessa_butterfly",True)):
        audited,ex,source_cells=audit_networks(raw,rev)
        subsets={
          "all_source_comparable_networks":lambda x:True,
          "no_zero_single_visit_pollen_means":lambda x:x["source_zero_pollen_means"]==0,
          "source_visitor_richness_at_least_3":lambda x:x["original_visitor_taxa"]>=3,
          "source_visitor_richness_at_least_5":lambda x:x["original_visitor_taxa"]>=5,
          "two_or_more_source_visitor_guilds":lambda x:x["original_visitor_guilds"]>=2,
          "only_one_source_visitor_guild":lambda x:x["original_visitor_guilds"]==1,
          "no_zero_SVE_and_multiple_visitor_guilds":
              lambda x:x["source_zero_pollen_means"]==0 and x["original_visitor_guilds"]>=2
        }
        sub={}
        for name,flt in subsets.items():
            entries=[x for x in audited if flt(x)]
            sub[name]={"source_summary":aggregate(entries),
                       "study_bootstrap":study_bootstrap(entries,1999 if name!=
                           "all_source_comparable_networks" else NBOOT)}
        # The multi- vs single-guild contrast is descriptive only; it
        # compares nonrandomly assembled ecological study populations.
        main=sub["all_source_comparable_networks"]["source_summary"]
        multi=sub["two_or_more_source_visitor_guilds"]["source_summary"]
        single=sub["only_one_source_visitor_guild"]["source_summary"]
        differential=(
            multi["functional_weight_vs_visitation_weight_TV_mean"]-
            single["functional_weight_vs_visitation_weight_TV_mean"]
            if multi["status"]=="source_descriptive" and single["status"]=="source_descriptive"
            else None)
        versions[label]={
           "schema":SCHEMA,"source_unmodified_original_623_SHA256":SHA,
           "source_matched_eligible_plant_networks":len(audited),
           "source_matched_original_studies":len({r["source_study"] for r in audited}),
           "source_quality_excluded_original_network_cells":ex,
           "subsets_and_source_study_bootstrap":sub,
           "multi_guild_minus_single_guild_average_TV_descriptive":differential,
           "example_high_functional_composition_mismatch":sorted(audited,
                key=lambda z:-z["functional_total_variation"])[:12],
           "example_high_top_visitor_prediction_error":sorted(audited,
                key=lambda z:-z["frequent_visitor_proxy_abs_error"])[:12],
           "observed_study_plant_networks_are_not_independent_ALAN_treatments":True
        }
        checked(len(audited)==72 and main["source_plant_network_cells"]==72,
                "verified primary population changed")
    return {
        "schema":SCHEMA,
        "original_page_paper_DOI":"10.1002/ajb2.1764",
        "original_taxon_source_Figshare_DOI":"10.6084/m9.figshare.16160265.v1",
        "original_csv_sha256":SHA,
        "source_original_taxon_rows":len(raw),
        "source_sensitivity_versions":versions,
        "estimand":"compositional total variation between observed visitor rate fractions and source q*V conspecific-stigma-pollen delivery fractions per same study/plant/period",
        "no_experimental_ALAN_light_treatment":True,
        "no_targeted_pollinator_removal_or_predator_intervention":True,
        "cannot_infer_stigma_to_seed_fitness":True,
        "no_moth_fly_effectiveness_transferred_to_wet_ALAN_flora":True,
        "method":"True source 33-study cluster descriptive bootstrap; q=0 reports are measured zero; no pseudo-replication across taxa or plant networks"
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    result=execute(get_source())
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("PAGE2021_ORIGINAL_FUNCTIONAL_BOTTLENECK_STUDY_CLUSTER",json.dumps({
       version:{
          "all":d["subsets_and_source_study_bootstrap"]["all_source_comparable_networks"],
          "nozero":d["subsets_and_source_study_bootstrap"]["no_zero_single_visit_pollen_means"]["source_summary"],
          "richness_ge5":d["subsets_and_source_study_bootstrap"]["source_visitor_richness_at_least_5"]["source_summary"],
          "multi_guild":d["subsets_and_source_study_bootstrap"]["two_or_more_source_visitor_guilds"]["source_summary"],
          "one_guild":d["subsets_and_source_study_bootstrap"]["only_one_source_visitor_guild"]["source_summary"],
          "multi_minus_single_TV":d["multi_guild_minus_single_guild_average_TV_descriptive"]
        } for version,d in result["source_sensitivity_versions"].items()
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
