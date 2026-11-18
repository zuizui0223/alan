#!/usr/bin/env python3
"""Reanalyse actual published Földesi et al 2020 Dryad SVD meta data (111 records).

Original article Földesi et al. 2021 J Appl Ecol DOI 10.1111/1365-2664.13798
Original Dryad https://doi.org/10.5061/dryad.05qfttf1n
Official API source version 91142, authoritative file Metadata.csv.
The published article references 127 source observations/28 studies/30
plants incl 33 flies; CURRENT published source version only 111 rows,
27 studies, 29 plants, 30 Diptera, 81 Hymenoptera. Report, do not repair
or silently inflate this difference.

Raw outcome lnRdelta = corrected log single-visit stigmatic pollen
deposition for focal fly or wild bee versus SAME-STUDY honey bee Apis
mellifera. It is NOT fly/moth, nor absolute grains per stigmatic visit.
For overlapping records in a study, correlations from shared honeybee
controls are unobserved. Therefore aggregate original study means first,
then resample STUDIES (not 111 taxon records) to get descriptive CIs.
No plant-specific light effect, wetland, flowering-time or moth efficiency
will be manufactured from cross-species honeybee comparisons.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import io
import json
import math
import random
import statistics
import urllib.request
import zipfile
from pathlib import Path

DOI="10.5061/dryad.05qfttf1n"
URL="https://datadryad.org/api/v2/versions/91142/download"
EXPECTED_SHA256="a46ca27b969962da46e579323879d721576c9b3abd10068052e0366d0b42a051"
SCHEMA="alan_foldesi2020_original_111_svd_study_unit_calibration_gate_v1"
ROWS=111
SEED=20261118
NBOOT=7999
ORIGINAL_PUBLICATION={"reported_observations":127,"reported_studies":28,
    "reported_plant_species":30,"reported_Diptera":33,"reported_wild_bees":94}
LU_SAME_FLOWER_SPECIES=[
    "Pedicularis kansuensis","Carum carvi",
    "Leontopodium leontopodioides","Anaphalis hancockii",
    "Dasiphora fruticosa","Cirsium setosum","Geum aleppicum"
]

def require(cond,why):
    if not bool(cond):raise ValueError("ALAN_FOLDESI_ORIGINAL_SVD_STOP: "+str(why))

def source_rows():
    request=urllib.request.Request(URL,headers={
       "User-Agent":"Open-science stigma deposition quantitative audit",
       "Accept":"application/zip,application/octet-stream,*/*"})
    with urllib.request.urlopen(request,timeout=65) as r:
        payload=r.read(15<<20)
    require(zipfile.is_zipfile(io.BytesIO(payload)),
            "Dryad published original must be true ZIP, not anti-bot HTML")
    with zipfile.ZipFile(io.BytesIO(payload)) as z:
        require(z.testzip() is None,"Original Dryad compressed CRC failed")
        src=[x for x in z.namelist() if x.endswith("Metadata.csv")]
        require(len(src)==1,"Original Dryad contains ambiguous file source")
        raw=z.read(src[0])
    require(hashlib.sha256(raw).hexdigest()==EXPECTED_SHA256,
            "Original exact source dataset changed version")
    import csv
    txt=raw.decode("utf-8-sig")
    frame=list(csv.DictReader(io.StringIO(txt),delimiter=";"))
    require(len(frame)==ROWS,"Original dataset must be 111 records, never pretend 127")
    expected={"studyID","Plant_sp","Plant_family","type","Pollinator",
        "Poll_order","Approach","lnR","V_lnR","lnRdelta","Var_lnRdelta"}
    require(expected.issubset(frame[0]),"No author study/plant/order/effect metadata")
    studies=set()
    plants=set()
    guilds=collections.Counter()
    rows=[]
    for r in frame:
        studies.add(r["studyID"])
        plants.add(r["Plant_sp"])
        guilds[r["Poll_order"]]+=1
        beta=float(r["lnRdelta"])
        variance=float(r["Var_lnRdelta"])
        require(math.isfinite(beta) and variance>0 and math.isfinite(variance),
                "Missing original corrected pollen effect/variance")
        rows.append({
           "study":r["studyID"],"plant":r["Plant_sp"],
           "guild":r["Poll_order"],
           "pollinator":r["Pollinator"],
           "plant_type":r["type"],
           "approach":r["Approach"],
           "effect_ln_focal_vs_honeybee":beta,
           "reported_ln_variance":variance
        })
    require(len(studies)==27 and len(plants)==29,"Published source study/plant set changed")
    require(dict(guilds)=={"Diptera":30,"Hymenoptera":81},
            "Published original source fly/bee counts differ from expected")
    return rows,{
        "original_DOI":DOI,"source_API_version":91142,
        "source_filename":"Metadata.csv",
        "original_metadata_sha256":EXPECTED_SHA256,
        "source_observations":len(rows),"source_studies":len(studies),
        "source_plant_taxa":len(plants),
        "source_effect_counts_by_guild":dict(guilds),
        "original_publication_description":ORIGINAL_PUBLICATION,
        "source_minus_published":{
           "n_observations":len(rows)-ORIGINAL_PUBLICATION["reported_observations"],
           "n_studies":len(studies)-ORIGINAL_PUBLICATION["reported_studies"],
           "n_plants":len(plants)-ORIGINAL_PUBLICATION["reported_plant_species"],
           "n_Diptera":guilds["Diptera"]-ORIGINAL_PUBLICATION["reported_Diptera"],
           "n_wild_bees":guilds["Hymenoptera"]-ORIGINAL_PUBLICATION["reported_wild_bees"]
        }
    }

def percentiles(values):
    vs=sorted(values)
    return [vs[int(.025*(len(vs)-1))],vs[int(.975*(len(vs)-1))]]

def group_summary(rows,group,rng):
    x=[v for v in rows if v["guild"]==group]
    by=collections.defaultdict(list)
    for v in x:by[v["study"]].append(v)
    effects=[statistics.fmean(r["effect_ln_focal_vs_honeybee"] for r in y)
             for y in by.values()]
    require(len(effects)>=3,"Too few original studies for source-study bootstrap")
    samples=[
       statistics.fmean(rng.choices(effects,k=len(effects)))
       for _ in range(NBOOT)
    ]
    c=percentiles(samples)
    vals=[r["effect_ln_focal_vs_honeybee"] for r in x]
    return {
        "guild":group,"source_original_obs":len(x),
        "independent_original_studies":len(by),
        "original_plant_species":len({r["plant"] for r in x}),
        "source_taxa_positive_effect_vs_honeybee":sum(v>0 for v in vals),
        "source_taxa_negative_effect_vs_honeybee":sum(v<0 for v in vals),
        "source_taxa_effect_min_max":[min(vals),max(vals)],
        "source_original_study_means_corrected_lnRatio":
             statistics.fmean(effects),
        "source_original_study_median_corrected_lnRatio":
             statistics.median(effects),
        "study_equal_geomean_pollen_deposition_vs_honeybee":
             math.exp(statistics.fmean(effects)),
        "study_cluster_bootstrap_95CI_ratio_vs_honeybee":[math.exp(z) for z in c],
        "study_equal_method_not_an_independent_new_meta_analysis":True,
        "study_IDs_are_resampling_units_not_taxon_entries":True,
        "source_effect_is_corrected_lnRdelta_focal_vs_Apis_mellifera":True
    }

def within_study_original_plant_pairs(rows):
    groups=collections.defaultdict(lambda:{"Diptera":[],"Hymenoptera":[]})
    for r in rows:groups[(r["study"],r["plant"])][r["guild"]].append(r)
    matched=[]
    for (study,plant),v in groups.items():
        if not v["Diptera"] or not v["Hymenoptera"]:continue
        # Study x plant shared comparator covariance cancels *in contrast*
        # of two log ratios, but the covariance needed for SE is not in
        # the original CSV. Hence do not invent SE/p for this match.
        a=statistics.fmean(r["effect_ln_focal_vs_honeybee"] for r in v["Diptera"])
        b=statistics.fmean(r["effect_ln_focal_vs_honeybee"] for r in v["Hymenoptera"])
        matched.append({
          "source_study":study,
          "source_plant_taxon":plant,
          "fly_original_taxon_obs":len(v["Diptera"]),
          "wild_bee_original_taxon_obs":len(v["Hymenoptera"]),
          "descriptive_ratio_fly_vs_wild_bee":
              math.exp(a-b),
          "not_fly_vs_moth_effectiveness":True,
          "no_valid_SE_due_shared_honeybee_controls":True
        })
    return {"source_study_plant_guild_pair_count":len(matched),
      "independent_studies_with_both_fly_and_wild_bee":
        len({r["source_study"] for r in matched}),
      "plant_taxa_with_both_fly_and_wild_bee":len({r["source_plant_taxon"] for r in matched}),
      "study_plant_pairs":matched,
      "no_inference_on_moth_pollen_deposition":True
    }

def match_lu(rows):
    plants={r["plant"] for r in rows}
    canonical={x.replace("_"," ").strip().lower() for x in plants}
    lu=[s.lower() for s in LU_SAME_FLOWER_SPECIES]
    return {
      "lu_seven_same_flower_moth_loss_fly_gain_species":LU_SAME_FLOWER_SPECIES,
      "exact_Lu_plant_name_matches_in_Dryad":[s for s in LU_SAME_FLOWER_SPECIES
                         if s.lower() in canonical],
      "number_exact_species_overlap":sum(s in canonical for s in lu),
      "genus_overlap":[s for s in LU_SAME_FLOWER_SPECIES
             if any(p.startswith(s.split()[0].lower()+" ") for p in canonical)],
      "fly_moth_same_plant_SVD_effects_available":False
    }

def execute():
    rows,source=source_rows()
    random_state=random.Random(SEED)
    guild={
        g:group_summary(rows,g,random_state)
        for g in ("Diptera","Hymenoptera")
    }
    return {
       "schema":SCHEMA,
       "source":source,
       "source_group_study_unit_descriptive":guild,
       "same_original_study_same_plant_fly_bee_strata":
             within_study_original_plant_pairs(rows),
       "match_original_Lu_7_plant_species":
             match_lu(rows),
       "original_pollinator_orders_present":sorted({r["guild"] for r in rows}),
       "fly_moth_direct_efficiency_missing":True,
       "original_raw_grains_per_visit_available":False,
       "meta_pollen_delivery_to_seed_set_reproduction_linked":False,
       "methodological_caveats":[
          "Published 127 vs publicly archived 111 rows; missing 16 not inferred from article",
          "Original Fly/Apis and wild bee/Apis log effect, NOT per visit absolute deposition or moth comparison",
          "Study-level equal weight descriptive study-resampling, no precision pooled random effects",
          "Underlying studies are not an ALAN light-exposure intervention",
          "Original fly observations share honeybee controls and some plants with other records; do not use simple 111-independent-effect CI",
          "No per-guild species effectiveness can be transferred to specific Lu plants without matched flowers and independent measurements",
          "Moth per-visit stigma pollen deposition is not in this publicly archived meta dataset"
       ],
       "original_publication_not_reproduced_in_full_due_source_archive_difference":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    result=execute()
    f=Path(args.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("FOLDESI_TRUE_ORIGINAL_STUDY_UNIT_FLY_POLLEN_DEPOSITION",json.dumps({
       "source":result["source"],
       "guilds":result["source_group_study_unit_descriptive"],
       "same_plant_study_support":{
          key:val for key,val in
          result["same_original_study_same_plant_fly_bee_strata"].items()
          if key!="study_plant_pairs"},
       "Lu_species_support":result["match_original_Lu_7_plant_species"]
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
