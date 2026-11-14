#!/usr/bin/env python3
"""Original-source ALAN lighting GEOMETRY reproduction analysis, independent of Lu et al.

Li et al. (2026) bioRxiv 10.64898/2026.07.10.737511; source public
author GitHub seedset.csv, original Git blob SHA pinned.
This reanalysis counts each true Fenswood experimental hayfield as an independent
block (three hayfields), never individual plant observations as randomized
independent lighting installations. No moth visitation or pollen delivery observed.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import math
import urllib.request
from pathlib import Path

SOURCE_REPO="Dongboli/experimental-data"
SOURCE_PATH="Artificial lights at night/seedset.csv"
SOURCE_URL="https://raw.githubusercontent.com/Dongboli/experimental-data/main/Artificial%20lights%20at%20night/seedset.csv"
GIT_BLOB_SHA1="b616aca5beeba1457549841477d61ca4f796a147"
GIT_AUTHOR_DATA_URL="https://github.com/Dongboli/experimental-data/blob/main/Artificial%20lights%20at%20night/seedset.csv"
TREATMENTS={"C":"no-light control","P":"isolated low-intensity point lights",
            "L":"25m low-intensity linear light array"}
SPECIES={"ML":"Matthiola longipetala","SN":"Silene nutans","ZC":"Zaluzianskya capensis"}

def require(v,msg):
    if not v: raise ValueError("ALAN_LI2026_RAW_GEOMETRY: "+msg)

def source_bytes():
    req=urllib.request.Request(SOURCE_URL,headers={"User-Agent":"ALAN-real-light-geometry-author-public-source/2026-11-14"})
    with urllib.request.urlopen(req,timeout=40) as resp:
        raw=resp.read(100_000)
    blob_header=f"blob {len(raw)}\0".encode()
    actual=hashlib.sha1(blob_header+raw).hexdigest()
    require(actual==GIT_BLOB_SHA1,"Original Li author Git blob changed: "+actual)
    return raw

def parse(raw):
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    required={"round","Treatment","replicate","species","pollinated","Seeds","Fruit","Label","Starting","No.uf"}
    require(len(rows)==79 and required<=set(rows[0]),"79 original records/columns changed")
    require({r["Treatment"] for r in rows}==set(TREATMENTS),"Treatment codes changed")
    require({r["species"] for r in rows}==set(SPECIES),"Plant codes changed")
    require({int(r["replicate"]) for r in rows}=={1,2,3},"Three true hayfield blocks missing")
    require({int(r["round"]) for r in rows}=={1,2,3},"Three experimental rounds missing")
    for r in rows:
        for col in ("Seeds","Fruit","pollinated"):
            if r[col]!="NA":
                x=float(r[col])
                require(x.is_integer() and x>=0, "Invalid original nonnegative outcome "+col)
        if r["pollinated"]!="NA":
            require(r["pollinated"] in {"0","1"},"Binary original seed-set changed")
            require(r["Seeds"]!="NA","Original binary reproductive outcome missing seed count")
            require(int(r["pollinated"])==int(float(r["Seeds"])>0),
                    "Binary seed set does not equal observed seed presence")
        else:
            require(r["Seeds"]=="NA","Missing binary but present seeds; inspect case")
    require(sum(x["pollinated"]=="NA" for x in rows)==1,"Original unknown outcome count changed")
    return rows

def summary(rows):
    valid=[r for r in rows if r["pollinated"]!="NA"]
    n=len(valid)
    yes=sum(int(r["pollinated"]) for r in valid)
    seeds=sum(int(r["Seeds"]) for r in valid)
    fruit=sum(int(r["Fruit"]) for r in valid)
    return {"all_original_source_rows":len(rows),
            "nonmissing_seedset_source_plants":n,
            "has_at_least_one_seed_plants":yes,
            "probability_observed":yes/n if n else None,
            "sum_real_seed_counts":seeds,
            "seeds_per_outcome_known_plant":seeds/n if n else None,
            "sum_real_fruit_counts":fruit,
            "fruits_per_outcome_known_plant":fruit/n if n else None}

def arm_groups(rows,groupfields):
    keys=sorted(set(tuple(r[k] for k in groupfields) for r in rows))
    out=[]
    for key in keys:
        rs=[r for r in rows if tuple(r[k] for k in groupfields)==key]
        out.append({"group":dict(zip(groupfields,key)),**summary(rs)})
    return out

def contrast(rows,treat,control="C"):
    a=summary([r for r in rows if r["Treatment"]==treat])
    b=summary([r for r in rows if r["Treatment"]==control])
    require(a["nonmissing_seedset_source_plants"]>0 and b["nonmissing_seedset_source_plants"]>0,
            "Treatment/control not observed")
    pa=a["probability_observed"]; pb=b["probability_observed"]
    yes1=a["has_at_least_one_seed_plants"]; no1=a["nonmissing_seedset_source_plants"]-yes1
    yes0=b["has_at_least_one_seed_plants"]; no0=b["nonmissing_seedset_source_plants"]-yes0
    return {"treatment":treat,"control":control,
            "raw_seed_presence_risk_difference":pa-pb,
            "raw_unadjusted_plants_as_rows_odds_ratio":yes1*no0/(no1*yes0)
                if no1*yes0>0 else None,
            "raw_seed_count_per_plant_difference":
                a["seeds_per_outcome_known_plant"]-b["seeds_per_outcome_known_plant"],
            "raw_fruit_count_per_plant_difference":
                a["fruits_per_outcome_known_plant"]-b["fruits_per_outcome_known_plant"],
            "treatment_source_summary":a,"control_source_summary":b,
            "not_independent_plants_as_randomized_field_replicates":True}

def field_comparisons(rows,treat):
    fields=sorted(set(int(r["replicate"]) for r in rows))
    out=[]
    for field in fields:
        subset=[r for r in rows if int(r["replicate"])==field]
        item=contrast(subset,treat)
        out.append({"independent_original_hayfield":field,
                    "seed_presence_difference":item["raw_seed_presence_risk_difference"],
                    "seed_mean_difference":item["raw_seed_count_per_plant_difference"],
                    "lit_positive":item["treatment_source_summary"]["has_at_least_one_seed_plants"],
                    "lit_n":item["treatment_source_summary"]["nonmissing_seedset_source_plants"],
                    "control_positive":item["control_source_summary"]["has_at_least_one_seed_plants"],
                    "control_n":item["control_source_summary"]["nonmissing_seedset_source_plants"]})
    signed=[v["seed_presence_difference"] for v in out]
    k=sum(x>0 for x in signed)
    m=sum(x<0 for x in signed)
    nposneg=k+m
    # Diagnostic for signed directions only, not a valid exact assignment test:
    # midseason lighting positions changed and n_sites=3.
    sign_fraction=2*sum(math.comb(nposneg,j) for j in range(0,min(k,m)+1))/(2**nposneg) if nposneg else None
    return {"field_rows":out,"independent_hayfields":len(fields),
            "hayfields_positive_difference":k,
            "hayfields_negative_difference":m,
            "hayfields_equal_difference":len(signed)-nposneg,
            "unweighted_hayfield_mean_risk_difference":sum(signed)/len(signed),
            "two_sided_sign_fraction_descriptive_not_original_assignment_p":sign_fraction,
            "leave_one_whole_hayfield_out_pooled_plant_risk_differences":[
                {"omitted_field":f,"difference":contrast(
                    [r for r in rows if int(r["replicate"])!=f],treat)["raw_seed_presence_risk_difference"]}
                for f in fields],
            "field_not_plant_is_primary_independent_replication":True,
            "midseason_field_treatment_location_reshuffle_cannot_reconstruct_exact_assignment":True}

def complete_plant_field_checks(rows):
    out=[]
    for species in SPECIES:
        per_species=[r for r in rows if r["species"]==species]
        for field in (1,2,3):
            subset=[r for r in per_species if int(r["replicate"])==field]
            arms={r["Treatment"] for r in subset}
            if {"L","C"}<=arms:
                eff=contrast(subset,"L")
                out.append({"plant_species":species,"field":field,
                            "unadjusted_array_minus_dark_seed_presence":
                                eff["raw_seed_presence_risk_difference"],
                            "unadjusted_array_minus_dark_seeds_per_plant":
                                eff["raw_seed_count_per_plant_difference"],
                            "complete_both_light_treatments":True})
            else:
                out.append({"plant_species":species,"field":field,
                            "complete_both_light_treatments":False,
                            "missing_arms":sorted({"L","C"}-arms)})
    positive=[x for x in out if x["complete_both_light_treatments"]
              and x["unadjusted_array_minus_dark_seed_presence"]>0]
    negative=[x for x in out if x["complete_both_light_treatments"]
              and x["unadjusted_array_minus_dark_seed_presence"]<0]
    ties=[x for x in out if x["complete_both_light_treatments"]
          and x["unadjusted_array_minus_dark_seed_presence"]==0]
    return {"strata":out,"supported_strata":len(positive)+len(negative)+len(ties),
            "positive_supported_strata":len(positive),
            "negative_supported_strata":len(negative),
            "tied_supported_strata":len(ties),
            "strata_share_hayfield_and_not_independent_plant_replicates":True}

def data_quality(rows):
    known=[r for r in rows if r["pollinated"]!="NA"]
    labeled=[r for r in rows if r["Label"]!="NA"]
    d=collections.Counter((r["species"],r["round"],r["Label"]) for r in labeled)
    return {"missing_binary_seedset_rows":len(rows)-len(known),
            "source_rows_with_missing_label":len(rows)-len(labeled),
            "source_rows_with_missing_id":sum(r["id"]=="NA" for r in rows),
            "source_species_round_label_keys_duplicated":sum(v>1 for v in d.values()),
            "source_rows_in_duplicated_species_round_label_keys":sum(v for v in d.values() if v>1),
            "label_keys_are_not_verified_independent_biological_plant_IDs":True,
            "no_imputation_for_one_unknown_point_light_outcome":True}

def main(out):
    out=Path(out)
    require(not out.exists(),"Do not overwrite an audit")
    raw=source_bytes()
    rows=parse(raw)
    overall=arm_groups(rows,["Treatment"])
    array=contrast(rows,"L")
    point=contrast(rows,"P")
    fields_array=field_comparisons(rows,"L")
    fields_point=field_comparisons(rows,"P")
    only_ML_SN=[r for r in rows if r["species"] in {"ML","SN"}]
    result={
        "schema":"ALAN_Li2026_real_public_79plants_three_field_blocks_geometry_seedset_20261114_v1",
        "source":{"published_preprint_DOI":"10.64898/2026.07.10.737511",
                  "preprint_not_peer_reviewed":True,
                  "original_source_repository":SOURCE_REPO,
                  "original_source_file":GIT_AUTHOR_DATA_URL,
                  "author_public_original_git_blob_sha1":GIT_BLOB_SHA1,
                  "downloaded_original_csv_byte_sha256":hashlib.sha256(raw).hexdigest(),
                  "original_data_folder_says_peer_review_only":True,
                  "actual_plant_source_rows":len(rows),
                  "experimental_field_blocks":3,
                  "original_lighting_classes":TREATMENTS,
                  "source_species_abbreviations":SPECIES,
                  "study_environment":"Bristol UK phytometer experimental hayfields 2022, not island or continental repeated study"},
        "data_quality":data_quality(rows),
        "treatment_summary":overall,
        "raw_effects":{"linear_array_versus_unlit":array,
                       "isolated_point_versus_unlit":point},
        "site_block_robustness":{"linear_array_versus_unlit":fields_array,
                                 "isolated_point_versus_unlit":fields_point},
        "plant_species_summary":arm_groups(rows,["species","Treatment"]),
        "array_versus_control_field_by_species":complete_plant_field_checks(rows),
        "both_all_three_fields_species_only_ML_SN":{
            "species_included":["Matthiola longipetala","Silene nutans"],
            "linear_array_versus_unlit":contrast(only_ML_SN,"L"),
            "site_block":field_comparisons(only_ML_SN,"L")},
        "functional_inference_guard":{
            "actual_seed_presence_and_seed_fruit_count_outcomes":True,
            "source_polllinated_column_is_seed_greater_than_zero_NOT_stigma_pollen":True,
            "true_nocturnal_moth_visits_observed":False,
            "nocturnal_and_diurnal_pollination_separately_excluded":False,
            "treatment_geometry_vs_lamp_count_identifiable":False,
            "three_fields_not_79_randomized_independent_lighting_replications":True,
            "point_lights_vs_linear_array_differ_total_bulb_count":True,
            "source_original_randomization_permutation_unavailable_within_season_relocation":True,
            "no_pollen_rescue_or_plant_past_island_shift_evidence":True,
            "plant_origin_is_not_island_mainland ALAN effect":True,
            "not_generalizable_global_moth_syndrome_coefficient":True,
        }
    }
    def check_close(a,b):
        require(abs(a-b)<1e-12,"Source unexpectedly moved from preaudited real counts")
    arms={x["group"]["Treatment"]:x for x in overall}
    require([(t,arms[t]["nonmissing_seedset_source_plants"],
              arms[t]["has_at_least_one_seed_plants"])
             for t in ("L","C","P")]==[("L",27,20),("C",26,8),("P",25,12)],
             "Original empirical 79-plant treatment/outcome counts drifted")
    check_close(array["raw_seed_presence_risk_difference"],20/27-8/26)
    require(fields_array["hayfields_positive_difference"]==3,
            "Expected original 3/3 independently positive field directions changed")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_LI2026_ORIGINAL_PUBLIC_79_PLANTS",json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    main(a.out)
