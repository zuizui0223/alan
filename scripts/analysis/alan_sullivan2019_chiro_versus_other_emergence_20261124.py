#!/usr/bin/env python3
"""Exploratory taxon-composition diagnostic within authenticated Sullivan2013 wetland.

Contrast author-labelled Chironomidae versus every *other counted label* in
the same AquaticEmergent traps. Other = non-Chironomidae catch, NOT verified
aquatic-origin insects or a coherent functional guild. Never count 4 lamps
and 4 repeated sampling dates as independent light-allocations.

For each of three physical wetland blocks g, pool the original four lighting
treatment catches with fixed equal weights, and calculate the cross-taxon
ratio-of-response-ratios (RRR):
 [sum_c Chiro_lit,g,c / Chiro_dark,g] /
 [sum_c Other_lit,g,c / Other_dark,g].
Factors 1/4 cancel (both classes sampled using same four source lamps).
The block-level 3 RRRs yield a DESCRIPTIVE t(df=2) interval and a sign
reference, never a causal mediation or source programme meta replicate.

Original DOI 10.1002/eap.1821 ; Dryad DOI 10.5061/dryad.gb9534p ;
source authentic TSV SHA256 385444ab318ec4b15edd3161c740aa7f4b2d9b973fe4c2ddf6c29cc56cdec484.
"""
from __future__ import annotations
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path
import statistics

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"scripts/analysis/alan_sullivan2019_chironomid_source_gate_20261124.py"
S=importlib.util.spec_from_file_location("source_sullivan",SOURCE)
mod=importlib.util.module_from_spec(S)
S.loader.exec_module(mod)

T2=4.302652729911275
LIGHTS=("LEDI","LEDH","HPSI","HPSH")

def evaluate(source):
    if (source.get("status")!="SOURCE_LABEL_RETRIEVED_NO_META_EFFECT" or
        source.get("original_2013_author_TSV_sha256")!=mod.DIGEST):
        raise ValueError("STOP original Sullivan2013 source identity missing")
    if source.get("chironomid_source_labels")!=["Chironomidae"]:
        raise ValueError("STOP Chironomidae original taxon name not matched")
    if source.get("independent_habitat_blocks")!=3 or source.get("original_stations")!=15:
        raise ValueError("STOP original physical experimental replication changed")
    expected={"DRK":356,"LEDI":1162,"LEDH":1295,"HPSI":1457,"HPSH":690}
    if source.get("chironomid_count_by_treatment")!=expected:
        raise ValueError("STOP source family treatment total changed")
    events={}
    for row in source["source_block_treatment_counts"]:
        g,t=str(row["group"]),row["light"]
        key=(g,t)
        if g not in ("1","2","3") or t not in ("DRK",)+LIGHTS or key in events:
            raise ValueError("STOP group/treatment key collision")
        total=int(row["all_aquatic_capture"])
        chiro=int(row["explicit_chironomid_recorded_capture"])
        if int(row["four_original_dates"])!=4 or total<0 or not 0<=chiro<=total:
            raise ValueError("STOP source group/treatment biological count inconsistent")
        events[key]={"chiro":chiro,"other":total-chiro,"total":total}
    if len(events)!=15:raise ValueError("STOP missing original 3x5 treatments")
    if sum(row["total"] for row in events.values())!=5372:
        raise ValueError("STOP source AquaticEmergent raw counts no longer 5372")
    overall={t:{k:sum(events[g,t][k] for g in ("1","2","3")) for k in ("chiro","other","total")}
             for t in ("DRK",)+LIGHTS}
    if sum(overall[t]["chiro"] for t in overall)!=4960:
        raise ValueError("STOP original aquatic emergence family total no longer 4960")

    logs=[];physical=[]
    for g in ("1","2","3"):
        dark=events[g,"DRK"]
        pooled={k:sum(events[g,t][k] for t in LIGHTS) for k in ("chiro","other","total")}
        # never make a pseudo zero count; stop rather than invent a correction.
        if any(x<=0 for x in (dark["chiro"],dark["other"],pooled["chiro"],pooled["other"])):
            raise ValueError(f"STOP: a physical block has zero family/complement count: {g}")
        rrr=(pooled["chiro"]/dark["chiro"])/(pooled["other"]/dark["other"])
        logs.append(math.log(rrr))
        physical.append({"block":g,"dark":dark,"four_lamp_aggregate":pooled,
                         "conditional_composition_RRR":rrr,
                         "log_conditional_composition_RRR":math.log(rrr)})
    avg=statistics.mean(logs)
    se=statistics.stdev(logs)/math.sqrt(3)
    ci=[math.exp(avg-T2*se),math.exp(avg+T2*se)]
    sign=0
    for choices in itertools.product((-1,1),repeat=3):
        q=abs(sum(s*y for s,y in zip(choices,logs))/3)
        if q>=abs(avg)-1e-12:sign+=1
    dark=overall["DRK"]
    illuminated={k:sum(overall[t][k] for t in LIGHTS) for k in ("chiro","other","total")}
    return {
        "doi":"10.1002/eap.1821",
        "data_doi":"10.5061/dryad.gb9534p",
        "source_original_sha256":mod.DIGEST,
        "source_and_taxon_matched":True,
        "true_independent_habitat_block_n":3,
        "source_all_2013_aquatic_emergence_count":5372,
        "source_chironomid_count":4960,
        "dark_capture_source_counts":dark,
        "four_lamp_aggregate_source_counts":illuminated,
        "dark_chironomid_proportion":dark["chiro"]/dark["total"],
        "four_lamp_chironomid_proportion":illuminated["chiro"]/illuminated["total"],
        "arithmetic_proportion_point_difference":illuminated["chiro"]/illuminated["total"]-dark["chiro"]/dark["total"],
        "pooled_count_ratio_of_RRs_descriptive":(illuminated["chiro"]/dark["chiro"])/(illuminated["other"]/dark["other"]),
        "three_true_physical_block_diagnostics":physical,
        "block_equal_weight_logRRR":avg,
        "block_equal_weight_RRR":math.exp(avg),
        "block_equal_weight_SE":se,
        "block_t2_descriptive_RRR_95":ci,
        "three_block_signflip_reference_p":sign/8,
        "source_site_level_sampling":"3 wetland groups x one dark + four original fixed light stations x four repeated catch dates",
        "never_equate_other_trap_catches_to_verified_other_aquatic_emergence":True,
        "composition_diagnostic_is_a_different_estimand_from_2013_emergence_meta":True,
        "existing_source_meta_result_not_overwritten":True,
        "new_meta_independent_k_added":0,
        "status":"SOURCE_PINNED_EXPLORATORY_COMPOSITION_NOT_REPRODUCTIVE_OR_PREDATION_MEDIATION",
        "limits":[
          "Other includes all source-recorded non-Chironomidae families in aquatic emergence traps, not an independently validated aquatic functional guild.",
          "The author had different spectra and intensity lamps sharing a dark station; four lamps do not supply independent meta programmes.",
          "The within-trap taxon contrast uses exact raw counts but wetland blocks n=3; t2 interval is descriptive.",
          "The light signal may reflect taxon attraction/retention rather than increased aquatic production.",
          "Separately fitting aquatic chironomid catch and other taxa does not establish insect demographic effects or plant reproductive mediation.",
        ]
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path,default=Path("results/alan_sullivan2019_chiro_vs_other_emergence_20261124.json"))
    args=a.parse_args()
    raw=mod.analyze(mod.get_original())
    out=evaluate(raw)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("SULLIVAN_CHIRO_V_OTHER_SOURCE",json.dumps({
      "source_SHA256_verified":True,"blockwise":out["three_true_physical_block_diagnostics"],
      "pooled":out["pooled_count_ratio_of_RRs_descriptive"],
      "equal_block_RRR":out["block_equal_weight_RRR"],
      "t2_interval":out["block_t2_descriptive_RRR_95"],
      "signflip":out["three_block_signflip_reference_p"],
      "meta_added":0
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
