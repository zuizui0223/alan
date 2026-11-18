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
    return {r["accepted"]:r for r in rows}

def execute():
    wet=source()
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
     "original_main_wet_species":len(wet),
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
         "original_Page_taxonlevel_effect_rows":len(lookup[n])}
        for n in matchedcases],
     "source_wet_reference_species_names_with_SVE":[
        {"species":n,"MOIST":float(wet[n]["indicator"]),
         "Page_guilds":sorted(set(
              r["pollinator.group"] for r in lookup[n])),
         "Page_taxonlevel_effect_rows":len(lookup[n])}
        for n in matchedrefs],
     "source_same_study_plant_moth_and_fly_pairs":len(five),
     "source_Pagemothfly_pair_plant_names":sorted(paired),
     "source_five_mothfly_plants_in_wetflora":sorted(paired&set(wet)),
     "source_five_mothfly_plants_in_Lu_alpine_sameflower_seven":
          sorted(paired&LU7),
     "source_wet_candidate_coverage_fraction":len(matchedcases)/58,
     "source_other_wet_reference_coverage_fraction":len(matchedrefs)/305,
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
