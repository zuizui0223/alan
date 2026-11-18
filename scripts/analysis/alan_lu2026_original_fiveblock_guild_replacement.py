#!/usr/bin/env python3
"""Author-origin 2026 Lu paired ALAN pollination source: five-block guild trade-off.

Reads byte-exact official "Date.xlsx" (the original author workbook filename)
from Europe PMC PMC13645961 with a minimal stdlib OOXML reader. The sheet
'visitation' has the *original* block, light, survey and diel columns. Five
lit/unlit spatial pairs (NOT visitor-level 2160/3381 independent sites) are
independent contrasts. Source-wide treatment totals and block-specific paired
effect sizes are audit outputs. Exact 2^5 within-block treatment label sign
flips are descriptive and have min two-sided p 0.0625, so NOT significant at
0.05 for any uniformly signed 5-block association.

Species visitation source is used only to estimate common-support coverage
from truly same plant-in-both-arms-in-same-block observations. Positive flower
visits do NOT measure plant fecundity or stigma-compatible pollen deposition.
No per-species p-value search. The authors' JPE seed-set paper may share plots.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import itertools
import json
import math
import statistics
from pathlib import Path

from scripts.acquisition.alan_lu2026_original_workbook_schema import source_open,sheet_data,ORIGINAL_WORKBOOK_SHA

SCHEMA="alan_lu2026_2024_original_five_block_diel_pollinator_guild_tradeoff_v1"
PRIMARY={"total":"abundance","Lepidoptera":"L_abundance",
         "Diptera":"D_abundance"}
DAYTIME={"Hymenoptera":"H_abundance","Coleoptera":"C_abundance",
          "Other":"O_abundance"}
GUILDS={"Lepidoptera":"Lepidoptera","Diptera":"Diptera","all":"all"}

def need(z,why):
    if not bool(z):raise ValueError("ALAN_LU2026_FIVE_BLOCK_SOURCE_STOP: "+why)

def rows(s):
    h=[str(z).strip() for z in s[0]]
    need(len(set(h))==len(h),"duplicate Excel field names")
    return [dict(zip(h, r+[""]*(len(h)-len(r)))) for r in s[1:]]

def num(v):
    try:
        if str(v).strip()=="":return None
        x=float(v)
        if not math.isfinite(x):return None
        return x
    except (TypeError,ValueError):return None

def normalize_time(s):
    v=str(s).strip().lower()
    if v in ("nocturnal","night","n","nighttime"):return "night"
    if v in ("diurnal","day","d","daytime"):return "day"
    return "unknown:"+v

def normalize_light(s):
    v=str(s).strip().lower()
    if v in ("lit","alan","light","lighted","illuminated"):return "lit"
    if v in ("unlit","dark","no alan","control","unilluminated"):return "dark"
    return "unknown:"+v

def validate_diel(table):
    need(len(table)==400,"original 400 visitation source rows changed")
    need(all({"time","sampling_time","treatment","block","repeat",
              "abundance","L_abundance","D_abundance",
              "H_abundance","C_abundance","O_abundance"} <= set(row) for row in table),
         "missing original author sheets input")
    original=Counter()
    parsed=[]
    for r in table:
        day=normalize_time(r["sampling_time"]);light=normalize_light(r["treatment"])
        need(day in ("day","night") and light in ("lit","dark"),
             "source diel or original treatment coding unknown")
        block=str(r["block"]).strip()
        survey_round=str(r["time"]).strip()
        rep=str(r["repeat"]).strip()
        need(block and survey_round and rep,"original block/round/replicate key absent")
        v={k:num(r[col]) for k,col in {**PRIMARY,**DAYTIME}.items()}
        need(all(x is not None and x>=0 and abs(x-round(x))<1e-8
                 for x in v.values()),"visits must be original nonnegative integer counts")
        original[(day,light,block)]+=1
        parsed.append({"day":day,"light":light,"block":block,
                       "round":survey_round,"repeat":rep,**v})
    blocks=sorted({v["block"] for v in parsed})
    need(len(blocks)==5,"exactly five independent original spatial blocks required")
    # In the author's original Data.xlsx, combinations of diel/treatment/
    # block/survey-round/repeat may occur more than once. Preserve ALL original
    # rows, but do not count those repeated keys as additional experimental
    # light installations. Report this issue as source-granularity QA.
    source_keys=Counter((x["day"],x["light"],x["block"],x["round"],x["repeat"])
                         for x in parsed)
    repeated=sum(v-1 for v in source_keys.values())
    whole_rows=Counter((x["day"],x["light"],x["block"],x["round"],x["repeat"],
                        *(x[k] for k in (*PRIMARY,*DAYTIME))) for x in parsed)
    exact_duplicate_rows=sum(v-1 for v in whole_rows.values())
    need(all(("day",arm,block) in original and ("night",arm,block) in original
             for block in blocks for arm in ("lit","dark")),"not all 5 treatment pairs share day and night")
    return parsed,{"n_raw_visits_frame_rows":len(table),
       "original_independent_block_pairs":len(blocks),
       "source_number_of_unique_diel_arm_block_round_repeat_keys":len(source_keys),
       "source_rows_with_reused_observation_keys":repeated,
       "source_completely_identical_full_rows_beyond_first":exact_duplicate_rows,
       "source_max_repeated_key_multiplicity":max(source_keys.values()),
       "repeated_keys_retained_as_author_original_not_extra_light_treatments":True,
       "balanced_source_group_record_count":len(set(original.values()))==1,
       "per_diel_arm_block_record_counts":sorted(set(original.values())),
       "independent_treatment_installations":10,
       "source_visitors_are_NOT_independent_treatments":True}

def exact_five_signflip(ds):
    need(len(ds)==5,"only five independent original site pairs")
    obs=statistics.mean(ds)
    sims=[statistics.mean([d*sgn for d,sgn in zip(ds,signs)])
          for signs in itertools.product((-1,1),repeat=5)]
    p=sum(abs(x)>=abs(obs)-1e-12 for x in sims)/len(sims)
    return {"paired_5_block_mean":obs,
            "descriptive_exact_32_paired_label_signflips_two_sided_p":p,
            "number_of_source_randomization_signflips":32,
            "not_a_large_n_significant_test":p>=.0625-1e-12,
            "independent_block_positive_effects":sum(x>0 for x in ds),
            "independent_block_negative_effects":sum(x<0 for x in ds),
            "minimum_two_sided_p_is_0p0625":True}

def fit_five_blocks(data):
    value=defaultdict(float)
    block_obs=defaultdict(int)
    for r in data:
        q=(r["day"],r["light"],r["block"])
        block_obs[q]+=1
        for guild in (*PRIMARY,*DAYTIME):
            value[(q,guild)]+=r[guild]
    blocks=sorted({x["block"] for x in data})
    out={"source_5_pair_diel_guilds":{}}
    for t in ("night","day"):
        groups={}
        for g in (*PRIMARY,*DAYTIME):
            lit=[value[((t,"lit",b),g)] for b in blocks]
            dark=[value[((t,"dark",b),g)] for b in blocks]
            # Geometric treatment ratio of independent paired *source counts*
            # +0.5 pseudocount for blocks with no visitors in a guild.
            delta=[math.log((l+.5)/(d+.5)) for l,d in zip(lit,dark)]
            x=exact_five_signflip(delta)
            x.update({
                "source_total_visits_lit":sum(lit),
                "source_total_visits_dark":sum(dark),
                "lit_to_dark_total_visits_ratio":
                    sum(lit)/sum(dark) if sum(dark)>0 else None,
                "geometric_mean_paired_block_lit_dark_ratio":math.exp(x["paired_5_block_mean"]),
                "original_pair_delta_log1pseudocount_min":min(delta),
                "original_pair_delta_log1pseudocount_max":max(delta)
            })
            groups[g]=x
        # Guild tradeoff, not simply a total animal visit reduction.
        moth=[math.log((value[((t,"lit",b),"Lepidoptera")]+.5)/
                       (value[((t,"dark",b),"Lepidoptera")]+.5)) for b in blocks]
        fly=[math.log((value[((t,"lit",b),"Diptera")]+.5)/
                       (value[((t,"dark",b),"Diptera")]+.5)) for b in blocks]
        relative=[m-f for m,f in zip(moth,fly)]
        joint=exact_five_signflip(relative)
        joint.update({
            "original_pairwise_lit_dark_relative_moth_to_fly_geomean_ratio":
               math.exp(joint["paired_5_block_mean"]),
            "source_definition":"log((lit moth+.5)/(dark moth+.5)) - log((lit fly+.5)/(dark fly+.5)) on each of 5 pairs",
            "pseudocount_sensitivity_needed_for_sparse_blocks":True,
            "replacement_in_effective_pollen_delivery_NOT_measured":True
        })
        out["source_5_pair_diel_guilds"][t]={"guilds":groups,
                                  "paired_moth_relative_to_fly":joint}
    out["original_plant_pollinator_visitor_total_vs_paper"]={
      "source_night_total_visits":sum(r["total"] for r in data if r["day"]=="night"),
      "source_day_total_visits":sum(r["total"] for r in data if r["day"]=="day"),
      "published_night_network_event_rows_2160":2160,
      "published_day_network_event_rows_3381":3381,
      "different_source_network_events_and_visitation_summary_possible":True
    }
    return out

def plant_source_support(s):
    data=rows(s)
    need(len(data)==1038,"source original 1038 species visitation rows changed")
    need({"plant","time","sampling_time","treatment","block","repeat",
          "all","Lepidoptera","Diptera"}.issubset(set(data[0])),
         "original plant-level species visit schema changed")
    specimen=defaultdict(lambda:{"lit":0,"dark":0})
    plant_guild=defaultdict(lambda:{"Lepidoptera":0,"Diptera":0,"all":0})
    species=set()
    for r in data:
        name=str(r["plant"]).strip()
        day=normalize_time(r["sampling_time"])
        arm=normalize_light(r["treatment"])
        block=str(r["block"]).strip()
        need(day in ("night","day") and arm in ("lit","dark"),
             "species-visit source missing original treatment/diel classification")
        need(name and block,"source plant name/block missing")
        need(all(num(r[k]) is not None for k in GUILDS),"source species visit count unknown")
        species.add(name)
        specimen[(name,day,block)][arm]+=1
        for guild in GUILDS:
            plant_guild[(name,day,block)][guild]+=num(r[guild])
    common=[(p,t,b) for (p,t,b),a in specimen.items()
            if a["lit"]>0 and a["dark"]>0]
    by_stage={}
    for t in ("night","day"):
        p={a for a,d,b in common if d==t}
        count=Counter(p for p,d,b in common if d==t)
        by_stage[t]={
          "species_with_same_block_both_light_arms":len(p),
          "species_with_at_least_two_independent_treatment_pairs":
             sum(n>=2 for n in count.values()),
          "species_with_all_five_block_treatment_pairs":sum(n==5 for n in count.values())
        }
    return {"original_source_species_visitation_rows":len(data),
            "distinct_original_plant_taxon_labels":len(species),
            "original_joint_source_species_block_diel_light_treatment_support":by_stage,
            "species_with_missing_light_arm_is_NOT_a_species_without_visits":True,
            "species_specific_fecundity_cannot_be_inferred":True}

def calculate():
    wb,rcode,src=source_open()
    need(src["source_workbook_sha256"]==ORIGINAL_WORKBOOK_SHA,"author original XLSX hash changed")
    sheets=sheet_data(wb)
    parsed,coverage=validate_diel(rows(sheets["visitation"]))
    out={"schema":SCHEMA,"source_original_DOI":"10.1002/ece3.74368",
         "source_original_PMC":"PMC13645961",
         "source_workbook_sha256":ORIGINAL_WORKBOOK_SHA,
         "source_paper_independent_treatment_control_spatial_pairs":5,
         "source_original_2024_field_experiment":True,
         "source_design":coverage,
         "source_fit":fit_five_blocks(parsed),
         "species_cross_arm_common_support":plant_source_support(sheets["species visitation"]),
         "not_identified":[
            "pollen delivery from a Dipteran visit compared with a Lepidopteran visit",
            "plant seed/fruit success from five paired visitation blocks",
            "independent experiment overlap with JPE paper DOI 10.1093/jpe/rtag160",
            "wet-habitat (EIV MOIST>7) moderation of five-block Tibetan alpine plant pollination"
         ]}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    result=calculate()
    f=Path(a.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("ALAN_LU2026_ORIGINAL_FIVE_SPATIAL_BLOCK_GUILD_REPLACEMENT",json.dumps({
      "source":result["source_design"],
      "night_guilds":{k:{kk:x[kk] for kk in (
           "source_total_visits_lit","source_total_visits_dark",
           "geometric_mean_paired_block_lit_dark_ratio",
           "descriptive_exact_32_paired_label_signflips_two_sided_p")}
           for k,x in result["source_fit"]["source_5_pair_diel_guilds"]["night"]["guilds"].items()},
      "night_moth_minus_fly":result["source_fit"]["source_5_pair_diel_guilds"]["night"]["paired_moth_relative_to_fly"],
      "source_totals":result["source_fit"]["original_plant_pollinator_visitor_total_vs_paper"],
      "species_common_support":result["species_cross_arm_common_support"]
    },sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
