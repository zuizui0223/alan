#!/usr/bin/env python3
"""Genuine original Macgregor2019 per-flower success->seed-damage stage ledger.

Only aggregate numbers and probabilities reported in peer-reviewed open PDF.
NO light-by-cage treatment arm cell counts, no inferred larval ALAN effects,
no fake joint correlation or pseudo-independent 417 randomized light plots.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,re,urllib.request,unicodedata
from pathlib import Path

URL="https://eprints.whiterose.ac.uk/id/eprint/149121/1/Macgregor_et_al_2019_Ecosphere.pdf"
SHA256="e0f83cc2086f2f494f175b0d1eee8495cf9153dd62f9ca7f35338ff9e5f1ad1c"
DOI="10.1002/ecs2.2550"
COUNTS={"flower_observed":417,"female_plants":169,"experimental_runs":4,
        "spatial_plot_positions":5,"seed_present_flowers":194,
        "larval_infested_seed_present_flowers":46,"uninfested_seed_capsules_quantified":148,
        "infested_day_access_only":7,"infested_fully_caged":2}
DAY=.376;NIGHT=.528;OPEN=.720

def ensure(x,reason):
    if not x:raise ValueError("ALAN_MACGREGOR2019_SOURCE_STAGE_GATE: "+reason)

def compute():
    z=COUNTS
    ensure(z["seed_present_flowers"]-z["larval_infested_seed_present_flowers"]==
           z["uninfested_seed_capsules_quantified"],"194-46=148 source mass balance")
    pred=1-(1-DAY)*(1-NIGHT)
    return {
      "schema":"ALAN_Macgregor2019_true_staged_reproduction_original_pdf_v1",
      "source_DOI":DOI,
      "verified_original_publication_pdf":URL,
      "required_original_pdf_sha256":SHA256,
      "published_primary_counts":z,
      "genuine_staged_descriptive_proportions":{
        "any_seeds_per_total_mature_flower":z["seed_present_flowers"]/z["flower_observed"],
        "not_seed_present_mature_flowers":z["flower_observed"]-z["seed_present_flowers"],
        "larval_infestation_per_seed_positive_flower":z["larval_infested_seed_present_flowers"]/z["seed_present_flowers"],
        "larval_infestation_per_total_mature_flower":z["larval_infested_seed_present_flowers"]/z["flower_observed"],
        "seed_count_measured_per_total_mature_flower":z["uninfested_seed_capsules_quantified"]/z["flower_observed"],
        "seed_count_measured_per_seed_positive_flower":z["uninfested_seed_capsules_quantified"]/z["seed_present_flowers"],
        "infested_day_only_of_all_infested":z["infested_day_access_only"]/z["larval_infested_seed_present_flowers"],
        "infested_fully_caged_of_all_infested":z["infested_fully_caged"]/z["larval_infested_seed_present_flowers"]
      },
      "published_original_diel_pollination_success":{
        "day_only":DAY,"night_only":NIGHT,"open":OPEN,
        "nocturnal_capability_over_open":NIGHT/OPEN,
        "marginal_night_access_loss_1_minus_day_over_open":1-DAY/OPEN,
        "model_independent_union_from_day_and_night":pred,
        "observed_open_minus_independence_union":OPEN-pred,
        "independence_union_NOT_a_statistical_test":True,
        "pooled_over_light_conditions_not_light_specific":True},
      "published_primary_results":{
        "light_x_cage_pollination_success":"NON_SIGNIFICANT_REPORTED_NO_EXACT_BETA_SE",
        "caging_vs_larval_infestation_chi_square":3.77,
        "caging_vs_larval_infestation_df":3,
        "caging_vs_larval_infestation_p":.288,
        "larvae_confirmed_to_be_adult_Hadena_species":False,
        "one_raised_larva_identified_as":"Lacanobia oleracea",
        "adult_fly_moth_visitor_counts_by_light_available":False
      },
      "strong_larval_source_warning":{
        "larvae_may_move_between_capsules_after_glasshouse_return":True,
        "no_separate_original_light_cage_infestation_analysis":True,
        "not_a_measured_field_light_effect_on_seed_predation":True,
        "conditioning_seed_quality_on_uninfested_capsules":True,
        "seed_count_for_infested_capsules_MISSING_not_zero":True,
        "full_bloom_per_flower_total_seed_output_NOT_ESTIMATED":True
      },
      "replication_caveat":{
        "flowers_not_independent_light_randomizations":True,
        "female_plants_reused_between_experimental_runs":True,
        "multiple_lamps_and_shared_plot_run_allocations":True,
        "true_level_of_light_randomization_spatial_plot_by_run":True},
      "P1_x_P3_new_original_interaction_estimates":0,
      "what_is_tested":"published null interaction light x caging in original paper; no refit",
      "note":"All reported percentages pooled over lighting, unknown treatment-specific effective pollen fractions and supplementary SE."}

def acquire_and_verify():
    req=urllib.request.Request(URL,headers={"User-Agent":"ALAN-primary-bioassay-staged-source-2026-11-15"})
    with urllib.request.urlopen(req,timeout=55) as r: raw=r.read(1_800_000)
    ensure(raw.startswith(b"%PDF"),"not genuine original PDF")
    ensure(hashlib.sha256(raw).hexdigest()==SHA256,"actual source PDF drift")
    import fitz
    pdf=fitz.open(stream=raw,filetype="pdf")
    pages=len(pdf);t="\n".join(p.get_text() for p in pdf)
    pdf.close()
    ensure(pages==17,"original PDF pagination changed")
    # The CC-BY original PDF uses Unicode ligatures (e.g. \uFB02owers),
    # and column-wrapped / hyphenated lines. Normalize typography before
    # checking exact, SHA-gated source anchor text.
    norm=" ".join(unicodedata.normalize("NFKC",t).lower().split())
    norm=re.sub(r"(?<=\\w)-\\s+(?=\\w)","",norm)
    checks={
      "417_flower_sample":bool(re.search(r"\b417\b.{0,100}\bflowers\b",norm)),
      "169_female_plants":bool(re.search(r"\b169\b.{0,55}\bplants\b",norm)),
      "194_pollinated_flowers":bool(re.search(r"\b194\b.{0,120}\bflowers\b",norm)),
      "46_infestations":bool(re.search(r"\b46\b.{0,60}(?:flowers|parasitized)\b|\b46 out of 194\b",norm)),
      "148_seeds_counted":bool(re.search(r"\b148\b.{0,100}\bflowers\b",norm)),
      "caging_light_interaction_non_significant":bool(re.search(r"no significant interaction.{0,80}between the two",norm)),
      "post_greenhouse_larval_cross_flower_movement":bool(re.search(r"move between.{0,50}flowers",norm)),
      "author_reported_pooled_diel_percentages":("52.8%" in norm and "37.6%" in norm and "72.0%" in norm)
    }
    ensure(all(checks.values()),"original PDF semantic anchors missing "+json.dumps(checks))
    return {"author_pdf_SHA256":hashlib.sha256(raw).hexdigest(),
            "pdf_pages":pages,"semantic_anchors":checks,
            "true_published_paper_retrieved":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    ensure(not args.out.exists(),"refuse replace")
    d=compute();d["original_pdf_integrity"]=acquire_and_verify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_MAC2019_REAL_SOURCE_REPRO_STAGE",json.dumps({
        "sha":d["original_pdf_integrity"]["author_pdf_SHA256"],
        "f":d["genuine_staged_descriptive_proportions"],
        "day_night":d["published_original_diel_pollination_success"],
        "interactions_new":0},ensure_ascii=False))
if __name__=="__main__":main()
