#!/usr/bin/env python3
"""Source-authentic Lu 2026 plant-specific moth vs fly redistribution analysis.

Read publication-original Data S1 at fixed XLSX checksum. The
author's focal 11-plant 'species visitation' original table has 1,038 source
records; it is not the COMPLETE 29 night/40 day plant community. The separate
'plant–pollinator network' worksheet contains non-tabular sparse/blank
identifier fields and is not silently imputed from previous row values.
For a plant to enter same-flower contrasts, both original lit/dark arms must have source records
in the *same spatial block and sampling period*, rather than infer unlisted
plant species as zeros or local absences. Compare (i) plant×block×diel
source support, (ii) plant×block×diel×survey-round/repeat strict support,
(iii) net moth losses vs fly gains received by the SAME flowering species.

This is a descriptive network rewiring audit. Flower pollen deposition,
quality of daytime vs nocturnal pollinator replacement and plant seed set
cannot be inferred from these source visits. Five spatial blocks are the
only original treatment replication; plant-block entries are not five
independent experimental studies.
"""
from __future__ import annotations
import argparse
import collections
import json
import math
import re
from pathlib import Path

from scripts.acquisition.alan_lu2026_original_workbook_schema import source_open,sheet_data,ORIGINAL_WORKBOOK_SHA
from scripts.analysis.alan_lu2026_original_fiveblock_guild_replacement import (
    rows,num,normalize_time,normalize_light,need,exact_five_signflip
)

SCHEMA="alan_lu2026_original_same_plant_guild_replacement_20261118_v1"
G=("all","Lepidoptera","Diptera","Hymenoptera","Coleoptera","Others")
RAW_SHEET="species visitation"
N_SOURCE=1038

def indexed_original(s):
    src=rows(s)
    need(len(src)==N_SOURCE,"original full plant×pollinator network source row count drift")
    need({"time","sampling_time","treatment","block","repeat","plant",
          *G}.issubset(src[0]),"original plant-visit matrix columns drift")
    parsed=[]
    for r in src:
        phase=normalize_time(r["sampling_time"])
        arm=normalize_light(r["treatment"])
        block=str(r["block"]).strip()
        survey_round=str(r["time"]).strip()
        repeat=str(r["repeat"]).strip()
        plant=str(r["plant"]).strip()
        need(phase in ("day","night") and arm in ("lit","dark"),
             "unknown source diel/treatment coding: sampling_time="
             +repr(str(r["sampling_time"])[:60])+", treatment="
             +repr(str(r["treatment"])[:60]))
        need(block and plant and survey_round and repeat,
             "original source plant, block or observation period missing")
        v={g:num(r[g]) for g in G}
        need(all(a is not None and a>=0 and a.is_integer() for a in v.values()),
             "original guild visitor counts must be nonnegative integer")
        need(v["all"]==sum(v[g] for g in G if g!="all"),
             "original plant×visitor summed guilds != original all")
        parsed.append({
            "phase":phase,"arm":arm,"block":block,
            "round":survey_round,"repeat":repeat,"plant":plant,**v})
    blocks=sorted(set(x["block"] for x in parsed))
    need(len(blocks)==5,"original five independently treated light/dark blocks changed")
    return parsed,blocks

def summarize(source,blocks,visit_rows):
    # Verify the original complete matrix agrees with the independently
    # supplied all-taxa visitation totals by diel stage. No forced rescaling.
    network=collections.Counter()
    main=collections.Counter()
    for r in source:
        for g in G:network[(r["phase"],g)]+=int(r[g])
    for r in visit_rows:
        t=normalize_time(r["sampling_time"])
        need(t in ("day","night"),"visitation sheet diel key missing")
        mapping={"all":"abundance","Lepidoptera":"L_abundance",
                 "Diptera":"D_abundance","Hymenoptera":"H_abundance",
                 "Coleoptera":"C_abundance","Others":"O_abundance"}
        for g,col in mapping.items():
            k=num(r[col])
            need(k is not None,"original visitation supporting source count missing")
            main[(t,g)]+=int(k)
    audit={"source_selected_focal_plant_visitation_rows":len(source),
        "source_focal_11_plants_NOT_complete_29_night_40_day_network":True,
        "source_original_full_network_1606_rows_not_used_due_sparse_identifiers":True,
        "original_visitation_rows":len(visit_rows),
        "original_unique_plant_labels":len({r["plant"] for r in source}),
        "original_plants_by_diel":{t:len({r["plant"] for r in source if r["phase"]==t})
                                   for t in ("day","night")},
        "source_focal_plant_counts_equal_whole_community_visitation":
            all(network[k]==main[k] for k in main),
        "source_focal_plants_are_subset_of_full_community":True,
        "per_diel_focal_plant_visits_vs_total_community":{
             t:{g:{"plant_network":network[(t,g)],"visitation_sheet":main[(t,g)]}
                for g in G} for t in ("night","day")},
        "network_rows_with_zero_total_visitors":
            sum(r["all"]==0 for r in source),
        "zero_visitor_source_row_not_evidence_of_local_plant_absence":True
    }
    return audit

def match_plants(source,blocks):
    groups=collections.defaultdict(lambda:{"lit":[],"dark":[]})
    groups_time=collections.defaultdict(lambda:{"lit":[],"dark":[]})
    for r in source:
        key=(r["plant"],r["phase"],r["block"])
        groups[key][r["arm"]].append(r)
        key2=(r["plant"],r["phase"],r["block"],r["round"],r["repeat"])
        groups_time[key2][r["arm"]].append(r)
    # The latter is stricter; duplicated published source keys are preserved
    # as extra source data rows, not separate treatment blocks.
    by_mode={}
    for mode,entries in (("same_plant_block_any_round",groups),
                         ("same_plant_block_round_repeat",groups_time)):
        matched=[]
        for key,arm in entries.items():
            if not arm["lit"] or not arm["dark"]:continue
            l={g:sum(r[g] for r in arm["lit"]) for g in G}
            d={g:sum(r[g] for r in arm["dark"]) for g in G}
            matched.append({
                "plant":key[0],"phase":key[1],"block":key[2],
                "round":key[3] if len(key)>=4 else None,
                "rep":key[4] if len(key)>=5 else None,
                "source_record_lit":len(arm["lit"]),
                "source_record_dark":len(arm["dark"]),
                "lit":l,"dark":d,
                "moth_delta":l["Lepidoptera"]-d["Lepidoptera"],
                "fly_delta":l["Diptera"]-d["Diptera"],
                "total_delta":l["all"]-d["all"]
            })
        by_mode[mode]=matched
    return by_mode

def source_results(rows0,blocks,author_taxon_map=None):
    if author_taxon_map is None:author_taxon_map={}
    from collections import defaultdict
    matched=match_plants(rows0,blocks)
    out={}
    for mode,rows in matched.items():
        result={"original_source_matched_row_keys":len(rows),
                "original_distinct_matched_plants":len({r["plant"] for r in rows}),
                "unique_source_plant_block_diel_cells":len({
                    (r["plant"],r["block"],r["phase"]) for r in rows}),
                "source_plants_same_block_treatment_pairs_by_diel":{
                    phase:len({(r["plant"],r["block"]) for r in rows if r["phase"]==phase})
                    for phase in ("night","day")},
                "source_not_listed_plant_never_called_zero":True,
                "counts_treatments_source_observations_not_independent_sites":True
                }
        n=[r for r in rows if r["phase"]=="night"]
        strict_positive=[r for r in n if r["moth_delta"]<0 and r["fly_delta"]>0]
        result["night_same_plant_moth_down_fly_up"]={
            "source_plant_block_round_repeat_rows":len(strict_positive),
            "source_unique_plant_block_pairs":len({
                (r["plant"],r["block"]) for r in strict_positive}),
            "source_distinct_plant_taxa":len({r["plant"] for r in strict_positive}),
            "source_independent_light_blocks_with_at_least_one_crossover":
                len({r["block"] for r in strict_positive}),
            "total_moth_decrease_on_these_same_plants":-sum(r["moth_delta"] for r in strict_positive),
            "total_fly_increase_on_these_same_plants":sum(r["fly_delta"] for r in strict_positive),
            "source_pairs_with_both_moth_and_fly_source_visits":
                sum(r["lit"]["Lepidoptera"]+r["dark"]["Lepidoptera"]>0 and
                    r["lit"]["Diptera"]+r["dark"]["Diptera"]>0 for r in n),
            "non_causal_replacement":True
        }
        # Not species-level replicated p-values. Catalogue the biggest
        # original plant contributors among those with >=2 paired blocks.
        plant=defaultdict(lambda:{"block_ids":set(),"moth_delta":0,
                                  "fly_delta":0,"total_delta":0,"source_rows":0})
        for r in n:
            x=plant[r["plant"]]
            x["block_ids"].add(r["block"])
            x["moth_delta"]+=r["moth_delta"]
            x["fly_delta"]+=r["fly_delta"]
            x["total_delta"]+=r["total_delta"]
            x["source_rows"]+=r["source_record_lit"]+r["source_record_dark"]
        contributors=[]
        for s,x in plant.items():
            if len(x["block_ids"])<2:continue
            contributors.append({
              "original_taxon_code":s,
              "source_author_R_abbreviated_latin_name":author_taxon_map.get(s),
              "source_spatial_block_pairs":len(x["block_ids"]),
              "moth_visits_lit_minus_dark":int(x["moth_delta"]),
              "fly_visits_lit_minus_dark":int(x["fly_delta"]),
              "all_visits_lit_minus_dark":int(x["total_delta"]),
              "original_source_plant_rows":x["source_rows"],
              "same_plant_moth_decrease_and_fly_increase":
                 x["moth_delta"]<0 and x["fly_delta"]>0})
        contributors.sort(key=lambda x:-(abs(x["moth_visits_lit_minus_dark"])+
                                        abs(x["fly_visits_lit_minus_dark"])))
        result["night_source_species_with_at_least_two_independent_block_pairs"]=len(contributors)
        result["night_largest_documented_species_changes_descriptive"]=contributors[:16]
        result["all_plant_rows_not_independent_experiments"]=True
        out[mode]=result
    return out


def read_author_taxon_abbreviations(r_source, plant_codes):
    """Use author R mapping only when explicitly present in exact supplement."""
    mapped={}
    for key,abbr in re.findall(
        r'"([a-z]{2,8})"\s*=\s*"([A-Z]\.\s*[a-z][a-z-]+)"', r_source):
        if key in mapped and mapped[key]!=abbr:
            raise ValueError("Conflicting original Lu author code to taxon mapping")
        mapped[key]=abbr
    return {x:mapped.get(x) for x in sorted(plant_codes)}

def calculate():
    wb,rcode,metadata=source_open()
    need(metadata["source_workbook_sha256"]==ORIGINAL_WORKBOOK_SHA,
         "source original 2026 supplemental XLSX changed")
    sheets=sheet_data(wb)
    source,blocks=indexed_original(sheets[RAW_SHEET])
    visitation=rows(sheets["visitation"])
    verify=summarize(source,blocks,visitation)
    plantcodes={x["plant"] for x in source}
    mapping=read_author_taxon_abbreviations(rcode,plantcodes)
    byplant=source_results(source,blocks,mapping)
    return {
        "schema":SCHEMA,
        "author_R_source_verified_plant_codes_to_latin_abbreviations":mapping,
        "plant_code_mapping_fraction_verified":
            sum(v is not None for v in mapping.values())/len(mapping),
        "original_source_workbook_SHA256":ORIGINAL_WORKBOOK_SHA,
        "DOI":"10.1002/ece3.74368",
        "original_independently_lit_control_spatial_pairs":5,
        "original_source_verification":verify,
        "original_same_plant_within_block_subsets":byplant,
        "source_moth_to_fly_functional_pollen_substitution_measured":False,
        "plant_specific_seed_set_measured":False,
        "original_JPE_2026_light_installation_independence_verified":False,
        "limitations":[
            "Focal 11 plant species visitation sheet is not a partition of the complete plant community visitation sheet",
            "Only plants observed in both source lighting arms of same spatial block enter individual flower contrasts; unmatched plant not assigned zero",
            "Strict matched time and survey replicate compared separately from any-round block-level support",
            "The original five independent blocks cannot be replaced by plant-by-survey pseudo-experimental replicates",
            "Plant visits do not identify pollinator effectiveness, compatible pollen deposition, or seed set",
            "Species labels remain author codes unless exactly mapped by original source R code"
        ]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    result=calculate()
    dest=Path(a.out)
    dest.parent.mkdir(exist_ok=True,parents=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("LU2026_ORIGINAL_SAME_FLOWERS_GUILD_REPLACEMENT",json.dumps({
        "source":result["original_source_verification"],
        "source_taxon_map":result["author_R_source_verified_plant_codes_to_latin_abbreviations"],
        "support":{k:{
            "matched_source_rowkeys":v["original_source_matched_row_keys"],
            "night_block_plant_pairs":
                v["source_plants_same_block_treatment_pairs_by_diel"]["night"],
            "night_mothdown_flyup":v["night_same_plant_moth_down_fly_up"],
            "strongest_species":v["night_largest_documented_species_changes_descriptive"][:10]
        } for k,v in result["original_same_plant_within_block_subsets"].items()}
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
