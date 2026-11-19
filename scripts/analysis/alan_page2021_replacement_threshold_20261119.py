#!/usr/bin/env python3
"""Source-grounded potential for *numerical* visitor replacement to restore pollen.

Uses independently published Page et al. 2021 original Figshare taxon-level
visitation and mean single-visit conspecific-pollen deposition data.
Observational static scenario only; NOT an ALAN experiment or observed
compensatory behavior, and NOT a wet-plant trait association.

Scenario: remove the originally most frequent visitor taxon from a source
study x plant x time x original metric cell. Distribute EXACTLY the same
number of lost visits among remaining taxa proportionally to their ORIGINAL
visit rates, holding each observed mean pollen-delivery q fixed.

Required visit multiplier kappa = q_lost / weighted_mean(q_remaining).
When kappa > 1, replacing 100% of lost visits cannot restore total pollen.
When others have zero observed pollen, replacement is unidentifiable/infinite.
Outcome is only potential conspecific-pollen delivery, not seed/fitness.

Supports original study bootstrap and zeros/taxon richness sensitivity;
the output never treats 72 source networks as 72 ALAN interventions.
"""
from __future__ import annotations
import argparse, collections, hashlib, io, json, math, random, statistics
from pathlib import Path
from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import get_source, num, VKEY, SHA
from scripts.analysis.alan_page2021_original_macro_functional_network_mismatch import candidate_groups
from scripts.analysis.alan_page2021_source_functional_bottlenecks_20261118 import audit_networks

SCHEMA="alan_original_page33study_72network_visit_replacement_threshold_v1"
SEED=20261119
TOL=1e-9

def require(ok,why):
    if not ok:
        raise RuntimeError("PAGE_ORIGINAL_FUNCTIONAL_REPLACEMENT_SOURCE_STOP: "+why)

def calculate(visits, efficiency):
    """Counterfactual without fitting effects from lighting or pollen fitness."""
    require(len(visits)==len(efficiency) and len(visits)>=2, "2+ taxa required")
    v=[float(x) for x in visits]
    q=[float(x) for x in efficiency]
    require(all(math.isfinite(x) and x>0 for x in v), "visits positive")
    require(all(math.isfinite(x) and x>=0 for x in q), "q nonnegative")
    total_v=sum(v)
    total_p=sum(a*b for a,b in zip(v,q))
    require(total_p>0, "positive original pollen")
    j=max(range(len(v)),key=lambda i:v[i])
    other_v=total_v-v[j]
    other_p=total_p-v[j]*q[j]
    require(other_v>0 and other_p>=-1e-8, "valid other taxa")
    other_p=max(0.,other_p)
    other_q=other_p/other_v
    # To restore pollen lost by removal, extra visit requirement:
    # V_new = V_lost*q_lost / qbar_remaining.
    impossible=(other_q==0 and q[j]>0)
    if impossible:
        multiplier=None
        retained=0.
    else:
        multiplier=q[j]/other_q if other_q>0 else 0.
        extra_p=v[j]*other_q
        retained=(other_p+extra_p)/total_p
    if impossible:
        # Equal-count replacement cannot deliver any pollen via remaining taxa.
        retained=0.
    # Best-case upper envelope: every lost visit goes to the single remaining
    # taxon with highest measured q. This is a bound, not observed choice.
    best_q=max(q[i] for i in range(len(q)) if i!=j)
    best_retained=(other_p+v[j]*best_q)/total_p
    best_case_failure=q[j]>best_q+TOL*max(q[j],best_q,1.)
    require(best_retained>=-TOL and math.isfinite(best_retained),
            "invalid best-case bound")
    require(retained<=best_retained+1e-8,
            "proportional replacement exceeded optimal replacement bound")
    require(retained>=-TOL and math.isfinite(retained), "invalid scenario retained")
    return {
      "taxa":len(v),
      "original_top_visitor_visit_share":v[j]/total_v,
      "original_top_visitor_pollen_share":v[j]*q[j]/total_p,
      "source_zero_q_taxa":sum(x==0 for x in q),
      "requires_more_than_equal_visits":impossible or multiplier>1+TOL,
      "needs_twice_lost_visits":impossible or multiplier>2+TOL,
      "needs_tenfold_lost_visits":impossible or multiplier>10+TOL,
      "requires_infinite_or_unavailable_visits":impossible,
      "required_new_visits_per_lost_visit":multiplier,
      "equal_visits_restored_pollen_ratio":retained,
      "even_best_remaining_taxon_fails_equal_visit_compensation":best_case_failure,
      "best_case_restored_pollen_ratio":best_retained,
      "best_case_functional_shortfall":max(0.,1-best_retained),
      "proportional_fails_but_best_case_could_succeed":
          (retained<1-TOL and best_retained>=1-TOL),
      "equal_visits_functional_shortfall":max(0.,1-retained),
      "equal_visits_functional_gain":max(0.,retained-1),
      "conditional_half_top_visits_shift_pollen_ratio":(1+retained)/2,
      "new_visits_not_actually_observed":True,
    }

def networks(source, revised):
    accepted, _, _=candidate_groups(source, revised)
    audited, _, _=audit_networks(source, revised)
    require(len(accepted)==72 and len(audited)==72, "72 source networks changed")
    observed=collections.defaultdict(list)
    for row in source:
        if revised and row["pollinator.group"].strip()=="moth" and row["pollinator.taxon"].strip()=="Vanessa carye":
            continue
        observed[tuple(row[k].strip() for k in VKEY)].append(row)
    result=[]
    for item in accepted:
        key=(item["source_study"],item["source_plant"],item["source_study_year_period"],
             item["original_SVE_measure"],"mean",item["original_visitation_rate_unit"])
        rows=observed.get(key,[])
        valid=[r for r in rows if num(r["visit.value"]) is not None
               and num(r["visit.value"])>0
               and num(r["effectiveness.average"]) is not None
               and num(r["effectiveness.average"])>=0]
        require(len(valid)==item["original_source_pollinator_taxa"],
                "original visitor taxa not aligned")
        q=[num(r["effectiveness.average"]) for r in valid]
        v=[num(r["visit.value"]) for r in valid]
        m=calculate(v,q)
        m.update({
          "source_study":item["source_study"],
          "source_plant":item["source_plant"],
          "source_year_period":item["source_study_year_period"],
          "source_visitor_guilds":len(set(r["pollinator.group"].strip() for r in valid)),
          "source_top_visitor_group":valid[max(range(len(valid)),key=lambda z:v[z])]["pollinator.group"].strip(),
        })
        result.append(m)
    require(len({x["source_study"] for x in result})==33, "33 source studies changed")
    return result

def mean(a):return statistics.fmean(a) if a else None
def median(a):return statistics.median(a) if a else None

def summary(rows):
    if not rows:return {"status":"NO_SOURCE_SUPPORT","networks":0}
    numeric=[x["required_new_visits_per_lost_visit"] for x in rows
             if x["required_new_visits_per_lost_visit"] is not None]
    clusters=collections.defaultdict(list)
    for r in rows:clusters[r["source_study"]].append(r)
    return {
       "status":"static_observational_counterfactual",
       "networks":len(rows),
       "independent_original_studies":len(clusters),
       "source_distinct_plant_names":len({r["source_plant"] for r in rows}),
       "most_frequent_visitor_requires_more_than_equal_visit_replacement_fraction":
          mean([int(x["requires_more_than_equal_visits"]) for x in rows]),
       "requires_more_than_equal_visit_replacement_count":
          sum(x["requires_more_than_equal_visits"] for x in rows),
       "requires_twice_visit_replacement_fraction":
          mean([int(x["needs_twice_lost_visits"]) for x in rows]),
       "requires_tenfold_visit_replacement_fraction":
          mean([int(x["needs_tenfold_lost_visits"]) for x in rows]),
       "unavailable_positive_q_remaining_networks":
          sum(x["requires_infinite_or_unavailable_visits"] for x in rows),
       "even_optimal_equal_count_replacement_fails_fraction":
          mean([int(x["even_best_remaining_taxon_fails_equal_visit_compensation"])
                for x in rows]),
       "even_optimal_equal_count_replacement_fails_count":
          sum(x["even_best_remaining_taxon_fails_equal_visit_compensation"]
              for x in rows),
       "proportional_fails_but_optimal_could_succeed_fraction":
          mean([int(x["proportional_fails_but_best_case_could_succeed"])
                for x in rows]),
       "best_case_pollen_shortfall_mean":
          mean([x["best_case_functional_shortfall"] for x in rows]),
       "required_visit_multiplier_median_when_finite":median(numeric),
       "required_visit_multiplier_75pct_when_finite":(
          float(sorted(numeric)[int(.75*(len(numeric)-1))]) if numeric else None),
       "relative_pollen_shortfall_after_full_numerical_replacement_mean":
          mean([x["equal_visits_functional_shortfall"] for x in rows]),
       "relative_pollen_shortfall_after_full_numerical_replacement_median":
          median([x["equal_visits_functional_shortfall"] for x in rows]),
       "relative_pollen_retained_after_full_numerical_replacement_median":
          median([x["equal_visits_restored_pollen_ratio"] for x in rows]),
       "study_balanced_fraction_requires_more_visits":mean([
          mean([int(x["requires_more_than_equal_visits"]) for x in z])
          for z in clusters.values()
       ]),
    }

BOOT_KEYS=("most_frequent_visitor_requires_more_than_equal_visit_replacement_fraction",
           "even_optimal_equal_count_replacement_fails_fraction",
           "proportional_fails_but_optimal_could_succeed_fraction",
           "best_case_pollen_shortfall_mean",
           "requires_twice_visit_replacement_fraction",
           "requires_tenfold_visit_replacement_fraction",
           "relative_pollen_shortfall_after_full_numerical_replacement_mean",
           "study_balanced_fraction_requires_more_visits")
def bootstrap(rows,reps=2999):
    clusters=collections.defaultdict(list)
    for x in rows:clusters[x["source_study"]].append(x)
    if len(clusters)<5:
        return {"status":"UNDER_FIVE_INDEPENDENT_STUDIES_STOP","studies":len(clusters)}
    keys=list(clusters)
    rng=random.Random(SEED)
    vals={k:[] for k in BOOT_KEYS}
    for _ in range(reps):
        batch=[x for key in rng.choices(keys,k=len(keys)) for x in clusters[key]]
        # Study-balanced proportions computed at draw-level by *original study*
        # would collapse repeats; use drawn cluster means separately.
        z=summary(batch)
        for k in BOOT_KEYS[:-1]:vals[k].append(z[k])
        chosen=rng # no extra sampling: use same draws via individual statistic below
        # The source study-balanced fraction is recomputed in summary(batch);
        # diagnostic CI is descriptive and counts selected study at most once.
        vals[BOOT_KEYS[-1]].append(z[BOOT_KEYS[-1]])
    out={}
    for k,x in vals.items():
        x.sort()
        out[k]=[x[int(.025*(reps-1))],x[int(.975*(reps-1))]]
    return {"status":"descriptive_original_study_cluster_bootstrap",
            "original_source_studies":len(keys),"replicates":reps,
            "descriptive_95pct":out}

def execute(source):
    require(len(source)==623, "original Page source rows changed")
    results={}
    for label,exclude_butterfly in [
        ("as_published",False),("verified_butterfly_taxon_excluded",True)]:
        net=networks(source,exclude_butterfly)
        subsets={
          "all":net,
          "no_source_zero_q":[a for a in net if a["source_zero_q_taxa"]==0],
          "at_least_five_visitors":[a for a in net if a["taxa"]>=5],
          "at_least_three_visitors":[a for a in net if a["taxa"]>=3],
          "at_least_two_guilds":[a for a in net if a["source_visitor_guilds"]>=2],
          "moth_as_most_frequent":[a for a in net if a["source_top_visitor_group"]=="moth"],
          "fly_as_most_frequent":[a for a in net if a["source_top_visitor_group"]=="fly"],
        }
        results[label]={
          "subsets":{k:{"summary":summary(v),
                        "bootstrap":bootstrap(v,2999 if k=="all" else 999)}
                     for k,v in subsets.items()},
          "studies_count":len({x["source_study"] for x in net}),
          "source_networks":len(net),
          "per_network_source_provenance_and_static_scenario":net,
        }
    return {
      "schema":SCHEMA,
      "original_data_doi":"10.6084/m9.figshare.16160265.v1",
      "original_paper_doi":"10.1002/ajb2.1764",
      "source_sha256":SHA,
      "source_rows":len(source),
      "scenarios":results,
      "lost_visitor":"most_frequent_original_visitor_taxon_within_each_study_plant_period",
      "restoration_rule":"reassign_same_number_of_visits_across_original_remaining_taxa_proportional_to_their_existing_visit_rates",
      "counterfactual_changes_only_visitor_allocation_not_per_visit_pollen_delivery_q":True,
      "no_alan_experiment_joined":True,
      "no_effect_of_alan_on_plant_fitness":True,
      "no_observed_behavioural_compensation":True,
      "cannot_generalize_to_wet_plant_ALAN_decline":True,
      "variance_of_source_q_unavailable":True,
      "scope":"static pollen-delivery potential, not causal mediation or pollination services after actual insect loss",
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    r=execute(get_source())
    where=Path(args.out);where.parent.mkdir(parents=True,exist_ok=True)
    where.write_text(json.dumps(r,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    compact={}
    for version,v in r["scenarios"].items():
        compact[version]={k:{
          "summary":q["summary"],
          "intervals":q["bootstrap"].get("descriptive_95pct"),
          "bootstrap_status":q["bootstrap"]["status"]
        } for k,q in v["subsets"].items()}
    print("PAGE33_STUDY_72_NETWORK_VISIT_REPLACEMENT_AUDIT "+
          json.dumps(compact,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
