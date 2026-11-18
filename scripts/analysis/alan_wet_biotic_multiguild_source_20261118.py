#!/usr/bin/env python3
"""Source-frozen macro inventory: wet flora => adult visitors => predator guilds.

This measures graph *documentary support*, not local abundance of insects,
light-specific attack rates, pollination of any plant, or causal mediation.
Original wet plant cohort came from alan/main: 58 selected decline candidates
and 305 other high-MOIST species. Source name-joined two-step GloBI/DoPI
adult prey links were built separately; 272 chain types are NOT independent
predator trials and many are observational iNaturalist associations.

Read the existing published summary files at their true Git blob hashes.
Never equate no interaction record with a known absence.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from collections import Counter,defaultdict

SCHEMA="alan_wet_flora_adult_visitor_predator_multiguild_source_audit_20261118_v1"
FILES={
  "wet_shore/shared_adult_chains.csv":"449a2b9b7af8be7a0d3ba54b2b42758e2212ba6e",
  "wet_partner_composition/species.csv":"bdca428269614c51e2c592f8757384d60f78438e",
  "wet_reproduction/habitat_labelled_moth_spiders.csv":"edee32f4cb8afc19bad638b53dcceb4e888fcaab"
}
EXPECTED_GUILDS={
    "spiders":91,
    "other_insects":73,
    "birds":47,
    "mantises":34,
    "dragonflies_damselflies":21,
    "reptiles":3,
    "amphibians":2,
    "other_animals":1
}

def require(cond,reason):
    if not bool(cond):raise ValueError("ALAN_WET_BIOTIC_MACRO_SOURCE_STOP: "+reason)
def gitsha(path):
    raw=Path(path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def load(root,relative):
    path=root/relative
    require(path.exists() and gitsha(path)==FILES[relative],
            "Original Git blob identity changed "+relative)
    with path.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))
def audit(root):
    chains=load(root,"wet_shore/shared_adult_chains.csv")
    plants=load(root,"wet_partner_composition/species.csv")
    moths=load(root,"wet_reproduction/habitat_labelled_moth_spiders.csv")
    require(len(chains)==272 and len(plants)==363 and len(moths)==8,
            "Original 272 two-step chain types/363 wet plants/8 labelled moth-spider edges changed")
    spids=[r["sp"] for r in plants]
    require(len(set(spids))==len(spids),"Original botanical species ID duplicated")
    require(sum(r["candidate"]=="True" for r in plants)==58,"Case cohort not frozen")
    require(sum(r["candidate"]=="False" for r in plants)==305,"Reference cohort not frozen")
    for r in chains:
        require(int(r["candidates"])>=1 and int(r["references"])>=1,
                "Each original shared visitor-predator pair must touch both plant groups")
        require(int(r["common_source_pairs"])>=1,
                "Source-matched visitor/diet evidence missing")
        require(int(r["source_pairs_without_inat"])<=int(r["common_source_pairs"]),
                "non-iNaturalist sources exceed total")
    classes=Counter(r["guild"] for r in chains)
    require(dict(classes)==EXPECTED_GUILDS,"Original guilds changed; do not fit emergent pathways")
    eligible=[r for r in chains if int(r["source_pairs_without_inat"])>0]
    require(len(eligible)==2,"Original source-independent of iNat support changed")
    mothplant={r["plant"] for r in moths}
    groups={r["accepted"]:r["candidate"] for r in plants
            if not any(x["accepted"]==r["accepted"] and x["candidate"]!=r["candidate"] for x in plants)}
    require(all(r["wet_site"]=="True" for r in moths),
            "Moth data not strictly source-labelled as riparian")
    require(all(r["guild"]=="spiders" for r in moths),
            "All 8 source-labelled adult chains used spiders; guild source changed")
    source_observations={
      "source_wet_flora_species":len(plants),
      "source_candidate_plants":58,
      "source_other_wet_plants":305,
      "documented_source_matched_adult_visitor_predator_pairs":len(chains),
      "distinct_named_adult_visitor_species":len({r["visitor"] for r in chains}),
      "distinct_named_predator_species":len({r["predator"] for r in chains}),
      "chains_with_at_least_one_non_iNaturalist_source_pair":len(eligible),
      "chains_with_only_iNaturalist_containing_source_pairs":len(chains)-len(eligible),
      "non_iNaturalist_source_pairs_present":{
           r["visitor"]+" / "+r["predator"]:{
              "predator_guild":r["guild"],
              "case_plant_species":int(r["candidates"]),
              "other_wet_plant_species":int(r["references"]),
              "independent_paired_sources_without_inaturalist":int(r["source_pairs_without_inat"])
           } for r in eligible
      },
      "all_source_pairs_are_evidence_of_global_links_not_same_field_plot":True
    }
    by_guild={}
    for g in sorted(classes):
        x=[r for r in chains if r["guild"]==g]
        by_guild[g]={
           "shared_visitor_predator_pair_types":len(x),
           "distinct_recorded_visitors":len({r["visitor"] for r in x}),
           "distinct_recorded_predators":len({r["predator"] for r in x}),
           "pair_types_with_non_iNaturalist_source_pair":
              sum(int(r["source_pairs_without_inat"])>0 for r in x),
           "cross_group_comparisons_with_nonoverlapping_reported_habitat_optima":
              sum(int(r["cross_group_disjoint_habitat_pairs"]) for r in x),
           "not_independent_study_replicates":True
        }
    riparian={
       "original_moth_to_spider_source_rows":len(moths),
       "source_distinct_botanical_species":len(mothplant),
       "source_candidate_botanical_species":sum(groups.get(n)=="True" for n in mothplant),
       "source_other_botanical_species":sum(groups.get(n)=="False" for n in mothplant),
       "source_distinct_moth_species":len({r["visitor"] for r in moths}),
       "source_distinct_spider_species":len({r["predator"] for r in moths}),
       "original_local_plant_spider_pollen_transfer_same_place_year_verified":False,
       "source_diet_edges_reliant_on_iNaturalist":
           sum("inaturalist" in r["dataset"].lower() for r in moths),
       "cannot_infer_functional_pollination_from_flower_visiting_or_moth_body_pollen":True
    }
    return {"schema":SCHEMA,"alan_main_plant_source_files":FILES,
       "plant_source_moisture_gate":">7 European EIV MOIST, not measured local water table",
       "two_step_biotic_graph_support":source_observations,
       "predator_guilds":by_guild,
       "documented_riparian_moth_adult_spider_chains":riparian,
       "mechanism_candidates":{
         "aquatic_subsidy_to_web_spiders":
             "source-supported by independent Parkinson experimental mesocosm; plant-pollination step missing",
         "non_web_flower_predators":
             "Misumena crab-spider flower visitor effects meta is light-unmanipulated; do not combine with Tetragnatha mesocosm treatment",
         "odonate_predators":
             "aquatic larval guild documented as adult predator in source graph, not verified same site-year",
         "reptiles_or_amphibians":
             "documented adult prey in a few source graph types; no proof of gecko pollination role at wet plant sites",
         "alternative_daytime_pollination":
             "Pladias insect-pollination classes are not species-level day/night functional redundancy"
       },
       "strict_not_estimated":[
          "ALAN x adult visitor trophic capture in the same field wet plant community",
          "ALAN x flower pollen delivery or reproductive success",
          "case-specific GloBI predator selection independent of the selection cohort",
          "interacting guild-level causal mediation or a cross-study pollination meta-coefficient"
       ]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--base",default="results/all_species_forest")
    p.add_argument("--out",required=True)
    args=p.parse_args()
    data=audit(Path(args.base))
    out=Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("ALAN_WET_BIOTIC_MULTIGUILD_SOURCE_INVENTORY",
          json.dumps({"support":data["two_step_biotic_graph_support"],
              "guilds":data["predator_guilds"],
              "riparian":data["documented_riparian_moth_adult_spider_chains"]},
              ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
