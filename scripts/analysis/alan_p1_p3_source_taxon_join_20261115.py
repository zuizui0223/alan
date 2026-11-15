#!/usr/bin/env python3
"""Source-locked P1 (diel exclusion) × P3 (ALAN) species crosswalk.

P1 original Kendall/Nicholson 1652 rows and P3 true Li preprint 79 plants.
Do not use population-level night dependence from a different study as a
within-population ALAN moderator. A true effect interaction needs BOTH
exposures on common ecological support, not merely the same taxon name.
"""
from __future__ import annotations
import argparse, csv,hashlib,io,json,math
from collections import Counter
from pathlib import Path
import pandas as pd
from scripts.analysis.alan_functional_dependence_20261115 import (
    SOURCE_SHA, match_panels, endpoint_subtype
)

LI_SOURCE_GIT_BLOB="b616aca5beeba1457549841477d61ca4f796a147"
LI_SOURCE_DOI="10.64898/2026.07.10.737511"
P1_SOURCE_DOI="10.1111/ele.70036"
LI_CODES={"ML":"Matthiola longipetala","SN":"Silene nutans",
          "ZC":"Zaluzianskya capensis"}
REPRO_METRICS={"fruit set","seed set"}

def require(ok,why):
    if not ok:raise ValueError("ALAN_P1_P3_REAL_SOURCE: "+why)

def li_blob_verify(raw):
    digest=hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+b"\x00"+raw).hexdigest()
    require(digest==LI_SOURCE_GIT_BLOB,"author Li Git blob mismatch "+digest)
    return digest

def li_original(raw):
    li_blob_verify(raw)
    data=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    require(len(data)==79,"Li true source changed from 79 original rows")
    require({r["species"] for r in data}==set(LI_CODES),"Li plant codes changed")
    require({r["Treatment"] for r in data}=={"L","P","C"},"Li exposure codes changed")
    require({r["replicate"] for r in data}=={"1","2","3"},"Li 3 field blocks missing")
    for r in data:
        if r["pollinated"]!="NA":
            require(r["pollinated"] in {"0","1"},"seed-presence flag")
            require(r["Seeds"]!="NA","known seed presence but missing seed count")
            require(int(r["pollinated"])==int(float(r["Seeds"])>0),"seed flag != Seeds>0")
    return data

def li_silene(data):
    rows=[r for r in data if r["species"]=="SN"]
    require(len(rows)==18,"source Silene nutans rows changed")
    arm={}
    for a in ("C","P","L"):
        q=[r for r in rows if r["Treatment"]==a]
        valid=[r for r in q if r["pollinated"]!="NA"]
        arm[a]={"rows":len(q),"valid_outcomes":len(valid),
                "seed_present":sum(int(r["pollinated"]) for r in valid),
                "total_seeds":sum(int(float(r["Seeds"])) for r in valid),
                "mean_seeds_per_known_plant":sum(int(float(r["Seeds"])) for r in valid)/len(valid)}
        require((arm[a]["rows"],arm[a]["valid_outcomes"])==(6,6),"SN six known outcomes per light arm")
    require([(arm[a]["seed_present"],arm[a]["total_seeds"]) for a in ("C","P","L")]==[
        (2,237),(5,481),(5,1636)],"Li original Silene nutans outcomes changed")
    site=[]
    for field in ("1","2","3"):
        t=[r for r in rows if r["replicate"]==field]
        n={}
        for a in ("C","P","L"):
            g=[r for r in t if r["Treatment"]==a]
            require(len(g)==2,"each Li SN field/treatment needs two original plants")
            n[a]={"known":2,"seed_present":sum(int(r["pollinated"]) for r in g),
                  "seeds":sum(int(r["Seeds"]) for r in g)}
        site.append({"field":int(field),"arms":n,
                     "linear_minus_dark_seed_presence_fraction":n["L"]["seed_present"]/2-n["C"]["seed_present"]/2,
                     "point_minus_dark_seed_presence_fraction":n["P"]["seed_present"]/2-n["C"]["seed_present"]/2,
                     "linear_minus_dark_seed_count":n["L"]["seeds"]-n["C"]["seeds"]})
    return {
      "source":"original author git SHA-verified 79-plant seedset.csv",
      "species":"Silene nutans","p3_preprint_DOI":LI_SOURCE_DOI,
      "n_independent_hayfield_blocks":3,"arms":arm,
      "linear_vs_dark_seed_presence_risk_difference":arm["L"]["seed_present"]/6-arm["C"]["seed_present"]/6,
      "linear_vs_dark_raw_seed_presence_ratio":arm["L"]["seed_present"]/arm["C"]["seed_present"],
      "linear_vs_dark_seed_presence_fraction_ratio":
        (arm["L"]["seed_present"]/6)/(arm["C"]["seed_present"]/6),
      "point_vs_dark_seed_presence_risk_difference":arm["P"]["seed_present"]/6-arm["C"]["seed_present"]/6,
      "linear_vs_point_seed_presence_risk_difference":arm["L"]["seed_present"]/6-arm["P"]["seed_present"]/6,
      "field_level_descriptive":site,
      "n_fields_linear_positive":sum(v["linear_minus_dark_seed_presence_fraction"]>0 for v in site),
      "n_fields_linear_equal":sum(v["linear_minus_dark_seed_presence_fraction"]==0 for v in site),
      "seed_count_linear_minus_dark_positive_field_count":
          sum(v["linear_minus_dark_seed_count"]>0 for v in site),
      "not_independent_plant_randomized_light_clusters":True,
      "lighting_not_geometry_only":"26 bulbs linear vs 4 bulbs point vs none",
      "raw_diel_access_measured_in_same_treated_plants":False,
      "p3_outcome_any_seeds_per_plant_NOT_p1_same_endpoint":True,
      "causal_P1xP3_moderation_established":False
    }

def audit(raw_p1,raw_li,source_registry):
    require(hashlib.sha256(raw_p1).hexdigest()==SOURCE_SHA,"P1 author original sha mismatch")
    p1_raw=pd.read_csv(io.BytesIO(raw_p1))
    require(len(p1_raw)==1652,"P1 treatment rows changed")
    panels,_=match_panels(p1_raw)
    require(len(panels)==369 and panels.study_ID.nunique()==116,"P1 369 study panel backbone changed")
    repro=panels[panels.treatment_effectiveness_metric.isin(REPRO_METRICS)].copy()
    require(len(repro)==340 and repro.study_ID.nunique()==113,"P1 reproductive source backbone changed")

    cat=json.loads(Path(source_registry).read_text())
    src=cat["publications"]
    require(len(src)==9 and len({x["doi"] for x in src})==9,"P3 paper registry changed")
    taxa=sorted(set(t for x in src for t in x["taxa"]))
    exact=sorted(set(repro.plant_species)&set(taxa))
    require(exact==["Silene latifolia","Silene nutans"],"Direct species source intersection changed")
    p3_linked=[p for p in src if set(p["taxa"])&set(exact)]
    require({p["id"] for p in p3_linked}=={"macgregor2019","boom2020","giavi2020","li2026"},
        "P3 linked paper assignments changed")

    detail=[]
    for name in exact:
        q=repro[repro.plant_species.eq(name)].sort_values(["study_ID","year_start","treatment_effectiveness_metric"])
        by=[]
        for r in q.itertuples(index=False):
            v={"study_ID":str(r.study_ID),"study_year":int(r.year_start),
               "endpoint":str(r.treatment_effectiveness_metric),
               "source_units":str(r.effectiveness_units),
               "source_day_only":float(r.day_mean),
               "source_night_only":float(r.night_mean),
               "source_open":float(r.open_mean),
               "night_exclusion_relative_loss":float(r.night_loss),
               "night_day_ratio":float(r.night_mean/r.day_mean) if r.day_mean>0 else None}
            by.append(v)
        detail.append({"species":name,"p1_studies":q.study_ID.nunique(),
                       "p1_comparable_reproductive_panels":len(q),
                       "p1_study_effects":by,
                       "P3_publications":[{"id":x["id"],"doi":x["doi"],
                                            "design":x["design"],"p3_endpoint":x["endpoint"]}
                                           for x in p3_linked if name in x["taxa"]]})

    alias=repro[repro.plant_species.eq("Silene alba")]
    require(alias.study_ID.nunique()==2 and len(alias)==3,"possible synonym sensitivity cohort drift")
    aliases={"Silene alba":{"p1_papers":int(alias.study_ID.nunique()),
        "p1_reproduction_panels":int(len(alias)),
        "relationship":"taxonomic synonym candidate only; no botanical author in compiled P1 rows",
        "promoted_to_exact_species_level":False,
        "authority":"https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:921549-1"}}

    li=li_silene(li_original(raw_li))
    nut=next(x for x in detail if x["species"]=="Silene nutans")
    fruit=next(x for x in nut["p1_study_effects"] if x["endpoint"]=="fruit set")
    seed=next(x for x in nut["p1_study_effects"] if x["endpoint"]=="seed set")
    require(abs(fruit["night_exclusion_relative_loss"]-(1-.09/.7))<1e-9,"P1 SN fruit dependence")
    require(abs(seed["night_exclusion_relative_loss"]-(1-16.544/284.314))<1e-8,"P1 SN seed dependence")
    mac=[x for x in detail if x["species"]=="Silene latifolia"][0]
    require(any(x["study_ID"]=="Macgregor_2019" for x in mac["p1_study_effects"]),
           "same-paper P1 study no longer in source")
    return {
      "schema":"ALAN_P1xP3_source_taxon_match_and_directional_counterexample_20261115_v1",
      "source_P1_SHA256":SOURCE_SHA,"source_P3_Li_git_blob_SHA1":LI_SOURCE_GIT_BLOB,
      "P1_experimental_treatment_rows":len(p1_raw),"P1_reproductive_panels":len(repro),
      "P1_reproductive_papers":repro.study_ID.nunique(),
      "P3_catalog_publications":len(src),"P3_catalog_species_unique":len(taxa),
      "P3_catalog_species_assignments":sum(len(x["taxa"]) for x in src),
      "P1_P3_exact_shared_species":exact,"n_exact_species":len(exact),
      "n_P3_publications_with_any_exact_P1_species":len(p3_linked),
      "P3_papers_with_exact_P1_species":sorted(x["id"] for x in p3_linked),
      "P1_exact_overlap_studies":sorted(set(repro[repro.plant_species.isin(exact)].study_ID)),
      "P1_exact_overlap_panels":len(repro[repro.plant_species.isin(exact)]),
      "P1_by_exact_species":detail,"P1_taxonomic_alias_candidate_heldout":aliases,
      "same_publication_ALAN_and_diel_summary":["Macgregor_2019"],
      "same_publication_ALAN_is_NOT_stratified_by_light_in_P1_table":True,
      "same_plant_population_ALAN_response_and_P1_diel_effect_independently_estimated":0,
      "separate_newer_within_study_factorial_P3_paper":
          {"id":"lu2026jpe","species":["Pedicularis kansuensis","Euphrasia regelii"],
           "same_study_diel_x_light_factorial_reported":True,
           "source_2025_P1_table_includes_these_species":False,
           "original_numeric_plant_level_effect_and_variance_unavailable":True},
      "true_source_species_crosspaper_contrast_example":{
        "P1_Belgium_2018_Silene_nutans":
          {"fruit_set_night_exclusion_loss":fruit["night_exclusion_relative_loss"],
           "seeds_per_plant_night_exclusion_loss":seed["night_exclusion_relative_loss"],
           "p1_publication":"Vanderplank_2020","P1_ALAN_measured":False},
        "P3_UK_2022_Silene_nutans":li,
        "direction_contradicts_universal_high_night_dependence_implies_ALAN_seed_decrease":True,
        "not_formal_falsification_of_P3_same_population_effect_modification":True},
      "main_P1xP3_interaction_estimated":False,
      "reason":"0 matched same-population light-stratified diel reproduction effects in the original author P1 file and 9-paper P3 registry; only two exact plant species and incompatible units among treatments.",
      "uncertainty":"2 P1 species joined to 4 P3 papers is a coverage count, not k=4 independent moderation estimates; same-paper Macgregor access means aggregate light treatments.",
      "global_taxonomic_background_may_be_incomplete":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--p1-source",type=Path,required=True)
    ap.add_argument("--li-source",type=Path,required=True)
    ap.add_argument("--catalog",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    require(not args.out.exists(),"refuse to overwrite source receipt")
    result=audit(args.p1_source.read_bytes(),args.li_source.read_bytes(),args.catalog)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_P1_P3_REAL_EXACT_TAXON_FUNCTIONAL_JOIN",json.dumps({
        "P1_papers":result["P1_reproductive_papers"],"P1_panels":result["P1_reproductive_panels"],
        "P3_papers":result["P3_catalog_publications"],"exact_species":result["P1_P3_exact_shared_species"],
        "P3_linked":result["P3_papers_with_exact_P1_species"],
        "Li_SN":result["true_source_species_crosspaper_contrast_example"]["P3_UK_2022_Silene_nutans"]["linear_vs_dark_seed_presence_risk_difference"],
        "joint_interaction_estimable":result["main_P1xP3_interaction_estimated"]},ensure_ascii=False))
if __name__=="__main__":main()
