#!/usr/bin/env python3
"""Compare Page original 623 single-visit records to frozen 363 wet ALAN flora.

This is COVERAGE / source transportability, not a plant-decline ALAN model.
Use literal exact accepted binomials, no synonyms/fuzzy nearest names.
Original main wet 363 (selected 58; reference 305), Git blob SHA frozen.
Source q*V pair comparisons refer to five *other* original plants,
and cannot be transferred to Lu seven alpine plants without direct traits.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from collections import Counter,defaultdict
from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import (
 get_source,compare,VKEY
)
SCHEMA="alan_page2021_original_sve_source_versus_363_wet_plant_coverage_v1"
MAIN_FILE=Path("results/all_species_forest/wet_partner_composition/species.csv")
MAIN_BLOB="bdca428269614c51e2c592f8757384d60f78438e"
LU7={"Pedicularis kansuensis","Carum carvi",
     "Leontopodium leontopodioides","Anaphalis hancockii",
     "Dasiphora fruticosa","Cirsium setosum","Geum aleppicum"}

def source():
    raw=MAIN_FILE.read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if actual!=MAIN_BLOB:raise RuntimeError("Original alan/main wet flora Git blob changed")
    with MAIN_FILE.open(newline="",encoding="utf-8") as f:rows=list(csv.DictReader(f))
    if len(rows)!=363 or sum(r["candidate"]=="True" for r in rows)!=58:
        raise RuntimeError("Frozen original 58/305 wet-case source mutated")
    # The 363 original source rows contain 2 same-binomial wet references
    # measured in distinct site regions. A dictionary must NOT silently
    # erase their source row membership or invent 363 different species.
    by=defaultdict(list)
    for r in rows:by[r["accepted"]].append(r)
    duplicates={k:v for k,v in by.items() if len(v)>1}
    for k,values in duplicates.items():
        assert len({v["candidate"] for v in values})==1,(
            "same-binomial candidate/ref source group conflict")
        assert len({v["indicator"] for v in values})==1,(
            "same-binomial source MOIST conflict")
    return ({k:v[0] for k,v in by.items()},{
       "source_original_rows":len(rows),
       "source_distinct_accepted_species":len(by),
       "source_duplicated_accepted_binomial_names":sorted(duplicates),
       "source_distinct_ALAN_candidates":
           len({r["accepted"] for r in rows if r["candidate"]=="True"}),
       "source_distinct_other_wet_species":
           len({r["accepted"] for r in rows if r["candidate"]=="False"}),
       "source_original_case_rows":58,
       "source_original_reference_rows":305
    })

def execute():
    wet,original_meta=source()
    rows=get_source()
    normalize=lambda r:r.replace("_"," ").strip()
    lookup=defaultdict(list)
    for r in rows:lookup[normalize(r["plant.species"])].append(r)
    page_plants=set(lookup)
    matched=sorted(page_plants&set(wet))
    matchedcases=[x for x in matched if wet[x]["candidate"]=="True"]
    matchedrefs=[x for x in matched if wet[x]["candidate"]=="False"]
    five=compare(rows,VKEY,True)["source_comparisons"]
    paired={normalize(x["source_plant"]) for x in five}
    return {
     "schema":SCHEMA,
     "main_wet_flora_source_git_blob_SHA1":MAIN_BLOB,
     "source_original_Page_taxon_table_sha256":
         "e28ebe99e0dfb27b73fbe6638e568bf98fa95b92c920e1d3d981e8eaea27519a",
     "original_main_wet_species":original_meta["source_original_rows"],
     "original_main_distinct_accepted_plant_binomials":
         original_meta["source_distinct_accepted_species"],
     "original_main_source_duplicates_not_imputed_as_unique_species":
         original_meta["source_duplicated_accepted_binomial_names"],
     "original_main_distinct_selected_candidate_species":
         original_meta["source_distinct_ALAN_candidates"],
     "original_main_distinct_reference_species":
         original_meta["source_distinct_other_wet_species"],
     "original_main_selected_ALAN_recording_rate_decline_cases":58,
     "original_main_other_high_MOIST_species":305,
     "Page_source_distinct_study_plant_labels":len(page_plants),
     "Page_original_plant_labels_exact_match_main_wet":len(matched),
     "Page_source_wet_candidate_cases_matched":len(matchedcases),
     "Page_source_wet_reference_species_matched":len(matchedrefs),
     "original_candidate_case_species_names_with_SVE":[
        {"species":n,"source_MOIST":float(wet[n]["indicator"]),
         "original_Page_pollinator_guilds":sorted(set(
              r["pollinator.group"] for r in lookup[n])),
         "original_Page_taxonlevel_effect_rows":len(lookup[n]),
         "original_Page_study_ids":sorted(set(
              r["studyID"] for r in lookup[n])),
         "original_Page_SVE_effectiveness_metrics":sorted(set(
              r["effectiveness.metric"] for r in lookup[n]))}
        for n in matchedcases],
     "source_wet_reference_species_names_with_SVE":[
        {"species":n,"MOIST":float(wet[n]["indicator"]),
         "Page_guilds":sorted(set(
              r["pollinator.group"] for r in lookup[n])),
         "Page_taxonlevel_effect_rows":len(lookup[n]),
         "Page_original_study_ids":sorted(set(
              r["studyID"] for r in lookup[n]))}
        for n in matchedrefs],
     "source_same_study_plant_moth_and_fly_pairs":len(five),
     "source_Pagemothfly_pair_plant_names":sorted(paired),
     "source_five_mothfly_plants_in_wetflora":sorted(paired&set(wet)),
     "source_five_mothfly_plants_in_Lu_alpine_sameflower_seven":
          sorted(paired&LU7),
     "source_wet_candidate_coverage_fraction":len(matchedcases)/58,
     "source_other_wet_reference_coverage_fraction":len(matchedrefs)/305,
     "source_unique_reference_binomial_coverage_fraction":
          len(matchedrefs)/original_meta["source_distinct_other_wet_species"],
     "no_ecological_absence_inferred_from_missing_single_visit_data":True,
     "no_species_level_ALAN_effect_estimated":True,
     "no_pollen_efficiency_transferred_between_nonidentical_species":True,
     "source_is_plant_trait_data_not_illumination_experiment":True
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    args=p.parse_args()
    d=execute()
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(d,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("PAGE_ORIGINAL_623_SVE_VS_WET_FLORA_CASE_COVERAGE",
          json.dumps(d,ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
