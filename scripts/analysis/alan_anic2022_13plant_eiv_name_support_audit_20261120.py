#!/usr/bin/env python3
"""Inspect official European EIV taxonomic ambiguities for Anic 13 species.

This is source-label audit ONLY. Unmatched EIV species are not counted dry.
Suffix/authorship/aggregate, genus and fuzzy-neighbour candidates are logged
for external expert/GBIF/WCVP taxonomic evidence; NO fuzzy value transfer.
"""
from __future__ import annotations
import argparse,collections,difflib,hashlib,json,re,unicodedata
from pathlib import Path
import pandas as pd

SCHEMA="alan_anic2022_13species_eiv_unmatched_taxonomic_candidate_audit_v1"
SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
def audit(source_path,trait_path):
    source=json.loads(Path(source_path).read_text())
    if source["schema"]!="alan_Anic2022_13plant_exact_EIV_MOIST_x_unmapped_light_codes_v1":
        raise ValueError("Anic original 13 reported source result absent")
    if hashlib.sha256(Path(trait_path).read_bytes()).hexdigest()!=SHA:
        raise ValueError("European EIV original source byte hash drift")
    d=pd.read_csv(trait_path,encoding="latin1",index_col=0,low_memory=False)
    d["MOIST"]=pd.to_numeric(d["MOIST"],errors="coerce")
    candidate=collections.defaultdict(set)
    for name,group in d.groupby(level=0).MOIST:
        if not isinstance(name,str):continue
        if group.dropna().nunique()<=1:
            value=group.dropna().iloc[0] if len(group.dropna()) else None
            candidate[name].add(None if value is None else float(value))
        else:candidate[name].add("CONFLICTING_VALUES")
    all_names=list(candidate)
    output=[]
    for source_row in source["original_species_measurements"]:
        name=source_row["species"]
        genus=name.split()[0]
        near=[x for x in all_names if x.casefold().startswith(name.casefold()+" ")]
        genus_sim=[x for x in all_names if x.split()[0].casefold()==genus.casefold()]
        close=sorted(genus_sim,key=lambda x:-difflib.SequenceMatcher(None,name.casefold(),x.casefold()).ratio())[:10]
        fuzzy=[{"eiv_source_name":x,"similarity":round(difflib.SequenceMatcher(None,name.casefold(),x.casefold()).ratio(),4),
                "source_trait_measure":sorted(v,key=str)[:4]}
               for x in close for v in [candidate[x]]
               if difflib.SequenceMatcher(None,name.casefold(),x.casefold()).ratio()>=.65]
        output.append({
          "original_source_plant":name,
          "original_source_exact_MOIST":source_row["MOIST_exact"],
          "exact_name_source_valid":source_row["MOIST_exact"] is not None,
          "EIV_source_rows_with_exact_binomial_plus_suffix":near[:15],
          "same_genus_EIV_rows":len(genus_sim),
          "same_genus_high_similarity_source_labels":fuzzy[:10],
          "no_suggested_mapping_automatically_accepted":True,
        })
    return {
     "schema":SCHEMA,"source_exact_EIV_SHA256":SHA,
     "source_Aniс_13_plants_verified":True,
     "original_13_plant_names":len(output),
     "original_exact_EIV_plant_names":sum(z["exact_name_source_valid"] for z in output),
     "source_ambiguous_or_unmatched":sum(not z["exact_name_source_valid"] for z in output),
     "candidate_label_audit_only":output,
     "no_trait_values_inferred_from_genus_or_fuzzy_match":True,
     "no_biological_moisture_effect_estimated":True,
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True);p.add_argument("--eiv",required=True);p.add_argument("--out",required=True)
    a=p.parse_args()
    result=audit(a.source,a.eiv)
    q=Path(a.out);q.parent.mkdir(exist_ok=True,parents=True)
    q.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("ANIC2022_EIV_13PLANT_UNMATCHED_TAXON_CANDIDATES "+json.dumps(result,sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
