#!/usr/bin/env python3
"""ALAN 2026 Lu original: does same-flower fly gain *target* plants losing night moths?

Observed ALAN treatments: FIVE spatial lit/dark pairs, original full-network
1,605 source rows from the 1,606-row supplemental "plant–pollinator network"
sheet. This is a conditional *cross-plant within-block* overlap test; it does
not infer realized pollen delivery, plant reproduction, or a wet-plant effect.

For 33 plant × spatial-block cells with both light arms and both diel phases,
freeze the exact 5 spatial blocks and group-wise marginal frequencies of
night moth decline and night/day fly increase. Under the conditional null,
reallocate fly-gain plant labels uniformly within the SAME block without
changing the number of fly-gain plants in that block. Each block's overlap
has a hypergeometric distribution. Convolve five block distributions to
get exact one/two-sided probabilities. Keep separate "day", "night", "either"
channels without adding incomparable visit counts across sampling periods.

Evidence gate: observed plant-block pair signs are NOT independent ALAN
experiments; causal claim about visitors selecting flowers after moth loss
would require longitudinal, within-flower mediation and effort normalization.
"""
from __future__ import annotations
import argparse, collections, json, math
from pathlib import Path

from scripts.acquisition.alan_lu2026_original_workbook_schema import source_open, sheet_data, ORIGINAL_WORKBOOK_SHA
from scripts.analysis.alan_lu2026_original_fullnetwork_sameflower_20261118 import parse_source
from scripts.analysis.alan_lu2026_original_sameplant_guild_replacement import match_plants

SCHEMA="alan_lu2026_five_spatial_block_conditional_fly_moth_overlap_v1"
SOURCE_SHA=ORIGINAL_WORKBOOK_SHA

def chk(cond,why):
    if not cond:raise RuntimeError("LU_ORIGINAL_BLOCK_CONDITIONAL_SOURCE_STOP: "+str(why))

def cells(source):
    paired=match_plants(source,[])["same_plant_block_any_round"]
    grouping={}
    for p in paired:
        key=(p["plant"],p["block"])
        t=grouping.setdefault(key,{})
        chk(p["phase"] not in t,"duplicate source plant block diel record")
        t[p["phase"]]=p
    output=[]
    for (plant,block),record in grouping.items():
        if not {"day","night"}.issubset(record):
            continue
        night,day=record["night"],record["day"]
        output.append({
            "plant":plant,"block":block,
            "moth_loss":night["moth_delta"]<0.,
            "night_fly_gain":night["fly_delta"]>0.,
            "day_fly_gain":day["fly_delta"]>0.,
            "either_fly_gain":night["fly_delta"]>0. or day["fly_delta"]>0.,
            "both_fly_gain":night["fly_delta"]>0. and day["fly_delta"]>0.,
            "original_night_moth_delta":night["moth_delta"],
            "original_night_fly_delta":night["fly_delta"],
            "original_day_fly_delta":day["fly_delta"]
        })
    chk(len(output)==33,"original 33 comparable same-plant same-block diel pairs changed")
    chk(len({x["block"] for x in output})==5,"original five independent light pairs changed")
    chk(sum(x["moth_loss"] for x in output)==21,"original moth-loss plant block cohort changed")
    chk(sum(x["moth_loss"] and x["either_fly_gain"] for x in output)==15,
        "previous source either-diel fly-increase count mismatch")
    return output

def hypergeom_pmf(n,losses,gains):
    chk(0<=losses<=n and 0<=gains<=n,"invalid group margins")
    lower=max(0,losses+gains-n)
    upper=min(losses,gains)
    d=math.comb(n,gains)
    return {k:math.comb(losses,k)*math.comb(n-losses,gains-k)/d
            for k in range(lower,upper+1)}

def convolve(a,b):
    out=collections.defaultdict(float)
    for ka,pa in a.items():
        for kb,pb in b.items():
            out[ka+kb]+=pa*pb
    return dict(out)

def exact_overlap(original,channel):
    b=collections.defaultdict(list)
    for x in original:b[x["block"]].append(x)
    dist={0:1.}
    observed=0
    total_loss=0
    total_gain=0
    byblock=[]
    for block,sub in sorted(b.items()):
        n=len(sub)
        m=sum(z["moth_loss"] for z in sub)
        k=sum(z[channel] for z in sub)
        intersection=sum(z["moth_loss"] and z[channel] for z in sub)
        observed+=intersection
        total_loss+=m
        total_gain+=k
        pmf=hypergeom_pmf(n,m,k)
        dist=convolve(dist,pmf)
        byblock.append({"block":block,"comparable_plant_cells":n,
            "night_moth_loss_plant_cells":m,
            "fly_gain_plant_cells":k,
            "both_observed":intersection,
            "null_expected_both":m*k/n})
    chk(abs(sum(dist.values())-1)<1e-9,"conditional null not normalized")
    exp=sum(k*v for k,v in dist.items())
    vv=sum((k-exp)**2*v for k,v in dist.items())
    tail=sum(p for k,p in dist.items() if k>=observed)
    lower=sum(p for k,p in dist.items() if k<=observed)
    # Probability-ordering exact two-sided test (not 2*one-sided).
    p_obs=dist.get(observed,0.)
    two=sum(p for k,p in dist.items() if p<=p_obs+1e-12)
    sorted_pmf=sorted(dist.items())
    def qtile(t):
        cs=0.
        for k,p in sorted_pmf:
            cs+=p
            if cs+1e-12>=t:return k
        return sorted_pmf[-1][0]
    return {
        "channel":channel,
        "comparable_plant_block_cells":len(original),
        "independent_original_spatial_light_pairs":len(b),
        "night_moth_loss_cells":total_loss,
        "fly_gain_cells":total_gain,
        "observed_moth_loss_and_fly_gain_cells":observed,
        "null_expected_intersection":exp,
        "observed_minus_conditional_null":observed-exp,
        "conditional_null_std":math.sqrt(vv),
        "conditional_null_95pct_integer_range":[qtile(.025),qtile(.975)],
        "exact_one_sided_enrichment_probability":tail,
        "exact_one_sided_depletion_probability":lower,
        "exact_probability_ordering_two_sided_p":min(1.,two),
        "original_source_block_margins":byblock,
        "null_is_exchangeability_of_plant_labels_within_spatial_block":True,
        "not_randomization_of_original_field_treatment_assignments":True,
    }

def execute():
    x,code,meta=source_open()
    chk(meta["source_workbook_sha256"]==SOURCE_SHA,"Lu Excel byte SHA mismatch")
    full,validation=parse_source(sheet_data(x)["plant–pollinator network"])
    same=cells(full)
    outcomes={
      k:exact_overlap(same,k)
      for k in ("night_fly_gain","day_fly_gain","either_fly_gain","both_fly_gain")
    }
    chk(outcomes["night_fly_gain"]["observed_moth_loss_and_fly_gain_cells"]==11,
        "night overlap mismatch")
    chk(outcomes["day_fly_gain"]["observed_moth_loss_and_fly_gain_cells"]==12,
        "day overlap mismatch")
    chk(outcomes["both_fly_gain"]["observed_moth_loss_and_fly_gain_cells"]==8,
        "both overlap mismatch")
    return {
      "schema":SCHEMA,
      "original_source_doi":"10.1002/ece3.74368",
      "original_author_excel_sha256":SOURCE_SHA,
      "source_validation":validation,
      "experimental_independent_light_pairs":5,
      "original_same_plant_both_diel_paired_cells":len(same),
      "original_focal_plant_names":len({x["plant"] for x in same}),
      "conditional_results":outcomes,
      "source_unseen_plant_visits_never_imputed_as_zeros":True,
      "source_diel_visitation_denominators_not_summed":True,
      "no_direct_functional_pollen_or_seed_delivery_tested":True,
      "no_causal_within_flower_fly_substitution_inferred":True,
      "cannot_explain_wet_plant_macro_bias_with_this_experiment":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    r=execute()
    z=Path(args.out);z.parent.mkdir(parents=True,exist_ok=True)
    z.write_text(json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("LU2026_SOURCE_5BLOCK_EXACT_FLY_MOTH_CONDITIONAL_OVERLAP "+
          json.dumps({k:{z:v[z] for z in (
          "comparable_plant_block_cells","night_moth_loss_cells",
          "fly_gain_cells","observed_moth_loss_and_fly_gain_cells",
          "null_expected_intersection","observed_minus_conditional_null",
          "conditional_null_95pct_integer_range",
          "exact_one_sided_enrichment_probability",
          "exact_probability_ordering_two_sided_p")}
          for k,v in r["conditional_results"].items()},sort_keys=True),flush=True)
if __name__=="__main__":
    main()
