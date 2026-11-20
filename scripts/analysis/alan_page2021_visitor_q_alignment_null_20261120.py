#!/usr/bin/env python3
"""Does observed frequency-effectiveness alignment amplify functional replacement failures?

Source: original Page et al 2021 Figshare 623 taxa, 33 original studies,
72 plant×study×period networks, original SHA pinned by existing get_source().
This is a conditional within-network permutation diagnostic of natural q/V
alignment, NOT a randomized ALAN or observed pollinator loss experiment.

Freeze each network's original visitation-rate vector, species names and
pollen-efficiency multiset. Permute the assignment of measured efficiency q
among pollinator taxa, leaving q and V distributions unchanged, and compare:
(a) the highest-visit taxon's q exceeds all remaining q (optimally irrecoverable)
(b) q_top exceeds remaining visit-weighted average q (proportional replacement
    insufficient).
Within-original-study balanced effects and study-cluster intervals are reported.
Additional within-guild taxon permutations diagnose whether alignment is
explained by broad visitor guild assignments.
"""
from __future__ import annotations
import argparse
import collections
import json
import math
import random
from pathlib import Path

from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import get_source, VKEY, num, SHA
from scripts.analysis.alan_page2021_original_macro_functional_network_mismatch import candidate_groups
from scripts.analysis.alan_page2021_replacement_threshold_20261119 import calculate, require

SEED=20261120
SCHEMA="alan_page2021_source33study_visitor_efficiency_alignment_null_v1"

def source_networks(raw, correction=True):
    accepted, _, _ = candidate_groups(raw, correction)
    require(len(accepted)==72,"original 72 networks changed")
    bykey=collections.defaultdict(list)
    for r in raw:
        if correction and r["pollinator.taxon"].strip()=="Vanessa carye" and r["pollinator.group"].strip()=="moth":
            continue
        bykey[tuple(r[k].strip() for k in VKEY)].append(r)
    networks=[]
    for a in accepted:
        key=(a["source_study"],a["source_plant"],a["source_study_year_period"],
             a["original_SVE_measure"],"mean",a["original_visitation_rate_unit"])
        eligible=[r for r in bykey[key]
                  if num(r["visit.value"]) is not None and num(r["visit.value"])>0
                  and num(r["effectiveness.average"]) is not None
                  and num(r["effectiveness.average"])>=0]
        require(len(eligible)==a["original_source_pollinator_taxa"],"source taxon mismatch")
        networks.append({
            "source_study":a["source_study"],
            "source_plant":a["source_plant"],
            "visits":[num(r["visit.value"]) for r in eligible],
            "q":[num(r["effectiveness.average"]) for r in eligible],
            "guilds":[r["pollinator.group"].strip() for r in eligible],
        })
    require(len({z["source_study"] for z in networks})==33,"study count changed")
    return networks

def thresholds(visits, q):
    """Exact two threshold indicators; target is first source maximum V."""
    require(len(visits)==len(q) and len(q)>=2,"2+ source visitors")
    idx=max(range(len(visits)),key=lambda i:visits[i])
    others=[j for j in range(len(q)) if j!=idx]
    qbar=sum(visits[j]*q[j] for j in others)/sum(visits[j] for j in others)
    return (int(q[idx]>max(q[j] for j in others)+1e-9*max(1.,max(q))),
            int(q[idx]>qbar+1e-9*max(1.,max(q))))

def permute(q, guilds, rng, restricted=False):
    if not restricted:
        out=list(q);rng.shuffle(out);return out
    byguild=collections.defaultdict(list)
    for i,g in enumerate(guilds):byguild[g].append(i)
    out=list(q)
    for indices in byguild.values():
        vals=[q[i] for i in indices]
        rng.shuffle(vals)
        for i,v in zip(indices,vals):out[i]=v
    return out

def network_null(z, rng, reps=1999, restricted=False):
    v,q,g=z["visits"],z["q"],z["guilds"]
    obs_best,obs_prop=thresholds(v,q)
    b,p=0,0
    changed=0
    # Independent conditional null replicates within this network.
    for _ in range(reps):
        shuffled=permute(q,g,rng,restricted)
        changed+=int(shuffled!=q)
        bb,pp=thresholds(v,shuffled)
        b+=bb;p+=pp
    return {
        "source_study":z["source_study"],
        "source_plant":z["source_plant"],
        "taxa":len(v),
        "original_guilds":len(set(g)),
        "guild_permutable":sum(len([x for x in g if x==t])>1 for t in set(g))>0,
        "observed_best_failure":obs_best,
        "observed_proportional_failure":obs_prop,
        "null_probability_best_failure":b/reps,
        "null_probability_proportional_failure":p/reps,
        "null_replicates_changed_q_assignment":changed,
        "conditional_best_failure_null_exact_if_unique_max_q":
            1./len(q) if sum(x==max(q) for x in q)==1 else 0.,
    }

def summarize(paired):
    bystudy=collections.defaultdict(list)
    for x in paired:bystudy[x["source_study"]].append(x)
    study_metrics=[]
    for study, group in bystudy.items():
        vals={}
        for typ in ("best","proportional"):
            observed=sum(r[f"observed_{typ}_failure"] for r in group)/len(group)
            exp=sum(r[f"null_probability_{typ}_failure"] for r in group)/len(group)
            vals[f"{typ}_observed"]=observed
            vals[f"{typ}_null"]=exp
            vals[f"{typ}_observed_minus_null"]=observed-exp
        study_metrics.append(vals)
    def avg(vals):return sum(vals)/len(vals)
    ans={"networks":len(paired),"source_studies":len(bystudy),
         "source_study_balanced":{},"source_network_balanced":{}}
    for typ in ("best","proportional"):
        obs=[a[f"observed_{typ}_failure"] for a in paired]
        exp=[a[f"null_probability_{typ}_failure"] for a in paired]
        ans["source_network_balanced"][typ]={
            "observed_fraction":avg(obs),
            "null_expected_fraction":avg(exp),
            "observed_minus_null_fraction":avg(obs)-avg(exp),
            "observed_count":sum(obs),
            "null_expected_count":sum(exp),
        }
        ans["source_study_balanced"][typ]={
            "observed_fraction":avg([x[f"{typ}_observed"] for x in study_metrics]),
            "null_expected_fraction":avg([x[f"{typ}_null"] for x in study_metrics]),
            "observed_minus_null_fraction":
                avg([x[f"{typ}_observed_minus_null"] for x in study_metrics]),
        }
    return ans,study_metrics

def study_bootstrap(study_metrics, reps=3999):
    rng=random.Random(SEED+987)
    names=["best_observed_minus_null","proportional_observed_minus_null"]
    sims={k:[] for k in names}
    for _ in range(reps):
        chosen=rng.choices(study_metrics,k=len(study_metrics))
        for key in names:
            sims[key].append(sum(z[key] for z in chosen)/len(chosen))
    return {k:[sorted(v)[int((reps-1)*.025)],
               sorted(v)[int((reps-1)*.975)]] for k,v in sims.items()}

def run(raw,permutations=1999, boot=3999):
    require(len(raw)==623,"source changed")
    source=source_networks(raw,True)
    outputs={}
    for label,restricted in [
        ("all_taxa_random_q_assignment",False),
        ("only_within_same_source_pollinator_guild",True)]:
        rng=random.Random(SEED + (100 if restricted else 0))
        paired=[network_null(x,rng,permutations,restricted) for x in source]
        aggregate,study= summarize(paired)
        aggregate["original_source_study_bootstrap_of_alignment_excess"]=study_bootstrap(study,boot)
        aggregate["source_count_with_permutable_taxa"]=sum(x["null_replicates_changed_q_assignment"]>0 for x in paired)
        aggregate["source_original_per_network_test_statistics"]=paired
        outputs[label]=aggregate
    return {
      "schema":SCHEMA,
      "source_page_paper_doi":"10.1002/ajb2.1764",
      "source_Figshare_original_csv_sha256":SHA,
      "source_taxon_rows":len(raw),
      "null_replications_per_network":permutations,
      "study_bootstrap_reps":boot,
      "results":outputs,
      "no_alan_treatment_in_source":True,
      "source_study_cluster_95_intervals_are_descriptive":True,
      "null_exchanges_q_among_taxa_not_a_physical_mechanism":True,
      "primary_question":"Is the most frequent visitor unusually efficient compared with randomly assigning the same empirical q values among visitor taxa?",
      "not_an_observed_plant_reproduction_or_compensation_result":True,
      "source_taxa_are_selected_honeybee_SVE_literature_not_global_random_sample":True,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    ap.add_argument("--permutations",type=int,default=1999)
    args=ap.parse_args()
    require(args.permutations>=99,"minimum null replication")
    r=run(get_source(),args.permutations)
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("PAGE33_Q_V_NULL_ALIGNMENT_AUDIT "+json.dumps({
      k:{g:v[g] for g in ("networks","source_studies",
                          "source_count_with_permutable_taxa","source_study_balanced",
                          "source_network_balanced",
                          "original_source_study_bootstrap_of_alignment_excess")}
      for k,v in r["results"].items()
    },sort_keys=True),flush=True)
if __name__=="__main__":
    main()
