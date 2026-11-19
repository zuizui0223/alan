#!/usr/bin/env python3
"""Stage-C original plant-specific moth-to-fly replacement REQUIRED q threshold.

Lu et al. (2026) DOI 10.1002/ece3.74368,
the already audited source in PR #65 (DO NOT COUNT ANOTHER study).
Europe PMC original nested supplementary Data workbook, original XLSX SHA
ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f.

This pipeline does NOT measure single-visit conspecific compatible pollen,
pollen limitation, fruit success or selfing. A q threshold is NECESSARY
conditional on visit counts and constant flower availability, NOT proof of
successful substitution, and 11 plants are dependent members of one 5-block
experiment. None of the 6 plants with both guild shifts is independent
lighting experiments.
"""
import argparse
from collections import defaultdict
import hashlib
import io
import json
import math
import urllib.request
import zipfile
from pathlib import Path

SOURCE="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles"
NESTED="ECE3-16-e74368-s001.zip"
WORKBOOK="Date.xlsx"
AUTHOR_SHA="ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f"
SCHEMA="alan_20261119_lu2026_original_plant_specific_stageC_moth_fly_required_q_not_observed_pollination_v1"

def source_bytes():
    req=urllib.request.Request(SOURCE,headers={
        "User-Agent":"ALAN-source-real-singlevisit-pollen-compensation-gates/2026-11-19",
        "Accept":"application/zip,*/*"})
    with urllib.request.urlopen(req,timeout=50) as f:
        data=f.read(20_000_000)
    with zipfile.ZipFile(io.BytesIO(data)) as outer:
        if NESTED not in outer.namelist():
            raise ValueError("Original Lu author supplement nested archive absent")
        inner=outer.read(NESTED)
    with zipfile.ZipFile(io.BytesIO(inner)) as actual:
        if WORKBOOK not in actual.namelist():
            raise ValueError("Original Lu 2026 source Date.xlsx missing")
        raw=actual.read(WORKBOOK)
    if hashlib.sha256(raw).hexdigest()!=AUTHOR_SHA:
        raise ValueError("Original Lu 2026 workbook identity SHA drift")
    return raw


def read_sheet(wb,name,max_rows=2500):
    if name not in wb.sheetnames:
        raise ValueError("Actual paper original workbook sheet missing")
    sh=wb[name]
    iterator=sh.iter_rows(min_row=1,max_row=min(max_rows,sh.max_row),values_only=True)
    header=next(iterator)
    header=[str(x).strip() if x is not None else "" for x in header]
    if len(set(k for k in header if k))!=sum(bool(k) for k in header):
        raise ValueError("Original workbook has repeated named source columns")
    rows=[]
    for raw in iterator:
        if raw[0] is None:
            continue
        row={k:v for k,v in zip(header,raw) if k}
        rows.append(row)
    if len(rows)>=max_rows-1:
        raise ValueError("Source original worksheet may be truncated")
    return rows


def num(v):
    if v is None or v=="":
        raise ValueError("True source visit count missing")
    value=float(v)
    if not math.isfinite(value) or value<0 or value!=int(value):
        raise ValueError("Invalid original visitor event count")
    return int(value)


def source_plant_q_thresholds(plant_rows,visit_rows,strict=True):
    necessary={"plant","sampling_time","treatment","block",
               "Lepidoptera","Diptera","all"}
    if not plant_rows or not necessary.issubset(plant_rows[0]):
        raise ValueError("Original Lu plant×visitor source schema missing")
    if not visit_rows or not {"sampling_time","treatment","L_abundance","D_abundance","block"}.issubset(visit_rows[0]):
        raise ValueError("Original Lu all-visit source schema missing")
    if strict and (len(plant_rows)!=1038 or len(visit_rows)!=400):
        raise ValueError("Original 1038 plant rows and 400 survey rows changed")
    plants=defaultdict(lambda:defaultdict(lambda:{"moth":0,"fly":0,"all":0,"blocks":set()}))
    raw_plant_names=set()
    for r in plant_rows:
        name=str(r.get("plant") or "").strip()
        if not name:raise ValueError("Original focal plant identity missing")
        raw_plant_names.add(name)
        period=r.get("sampling_time")
        if period not in ("Nocturnal","Diurnal"):raise ValueError("Source diurnal period changed")
        t=r.get("treatment")
        if t not in ("Lit","Unlit"):raise ValueError("Source light treatment changed")
        counts=(num(r["all"]),num(r["Lepidoptera"]),num(r["Diptera"]))
        if r["sampling_time"]!="Nocturnal":continue
        block=int(r["block"])
        if block not in (1,2,3,4,5):raise ValueError("Five original paired field blocks changed")
        x=plants[name][t]
        for key,value in zip(("all","moth","fly"),counts):
            x[key]+=value
        x["blocks"].add(block)
    if strict and len(raw_plant_names)!=11:
        raise ValueError("Original paper focal 11 documented plant taxa changed")
    original=defaultdict(lambda:{"moth":0,"fly":0,"n":0})
    for r in visit_rows:
        if r["sampling_time"]!="Nocturnal":continue
        t=r["treatment"]
        if t not in ("Lit","Unlit"):raise ValueError("Published original source light labels changed")
        original[t]["moth"]+=num(r["L_abundance"])
        original[t]["fly"]+=num(r["D_abundance"])
        original[t]["n"]+=1
    if strict and (original["Unlit"]["moth"],original["Lit"]["moth"],
                  original["Unlit"]["fly"],original["Lit"]["fly"])!=(838,357,54,151):
        raise ValueError("Original Lu 5-block guild count source drift")
    if strict and (original["Unlit"]["n"]!=100 or original["Lit"]["n"]!=100):
        raise ValueError("Original 100 night surveys per light arm changed")
    threshold=[]
    observed_both=0
    for name in sorted(raw_plant_names):
        groups=plants[name]
        if "Lit" not in groups or "Unlit" not in groups:
            # Don't manufacture zero for undocumented plant×light visits.
            continue
        lit=groups["Lit"];dark=groups["Unlit"]
        loss=dark["moth"]-lit["moth"]
        gain=lit["fly"]-dark["fly"]
        eligible=loss>0 and gain>0
        if eligible:observed_both+=1
        if eligible:
            needed=loss/gain
            threshold.append({
                # Taxon identity public but not needed for aggregate inference;
                # do not duplicate source original rare taxon x plot rows.
                "original_source_plant_anon_rank":None,
                "loss_original_moth_visitor_events":loss,
                "gain_original_fly_visitor_events":gain,
                "required_fly_over_moth_single_visit_effectiveness_ratio_if_other_factors_constant":needed,
                "equal_q_is_theoretically_sufficient_for_lost_moth_event_count":needed<=1,
                "plant_q_measurement_observed":False,
                "light_exposed_flower_availability_verified":False,
                "no_actual_pollen_fruit_seed_compensation_verified":True,
            })
    if strict and observed_both!=6:
        raise ValueError("Original 6 of 11 focal plants both guild signs changed")
    threshold.sort(key=lambda x:x["required_fly_over_moth_single_visit_effectiveness_ratio_if_other_factors_constant"])
    for i,x in enumerate(threshold,1):x["original_source_plant_anon_rank"]=i
    total_loss=original["Unlit"]["moth"]-original["Lit"]["moth"]
    total_gain=original["Lit"]["fly"]-original["Unlit"]["fly"]
    if total_loss!=481 or total_gain!=97:
        raise ValueError("Original community moth and fly changes mismatch")
    needs=total_loss/total_gain
    ratios=[x["required_fly_over_moth_single_visit_effectiveness_ratio_if_other_factors_constant"] for x in threshold]
    def quantile(p):
        if not ratios:return None
        u=(len(ratios)-1)*p
        lo=int(math.floor(u)); hi=int(math.ceil(u))
        return (1-(u-lo))*ratios[lo]+(u-lo)*ratios[hi]
    return {
        "schema":SCHEMA,
        "original_article_DOI":"10.1002/ece3.74368",
        "source_original_workbook_sha256":AUTHOR_SHA,
        "source_total_survey_rows":len(visit_rows),
        "source_plant_positive_visitation_rows":len(plant_rows),
        "source_field_blocks":5,
        "source_original_focal_plant_species":len(raw_plant_names),
        "original_documented_nocturnal_moth_down_fly_up_plant_n":observed_both,
        "community_moth_visits_dark_and_lit":[original["Unlit"]["moth"],original["Lit"]["moth"]],
        "community_fly_visits_dark_and_lit":[original["Unlit"]["fly"],original["Lit"]["fly"]],
        "community_visit_event_losses_moth":total_loss,
        "community_visit_event_gains_fly":total_gain,
        "community_conditional_qfly_qmoth_threshold_if_all_visits_act_on_same_equivalent_plants":needs,
        "focal_six_plant_required_efficiency_ratio_descriptive_summary":{
            "n_plant_taxa_with_both_source_changes":len(ratios),
            "minimum":min(ratios) if ratios else None,
            "median":quantile(.5),
            "maximum":max(ratios) if ratios else None,
            "count_where_equal_q_visits_would_suffice":sum(v<=1 for v in ratios),
            "count_where_fly_q_must_exceed_moth_q":sum(v>1 for v in ratios),
            "theoretical_only_NOT_compensated_plant_proportion":True,
        },
        "focal_plant_thresholds_aggregated_no_original_species_or_site_ID":threshold,
        "guard":{
            "source_11_plant_taxa_not_all_29_nocturnal_bloom_taxa":True,
            "six_plant_thresholds_one_five_block_experiment_NOT_six_independent_light_studies":True,
            "not_equal_per_flower_visitation_effort_or_flower_availability":True,
            "cross_plant_community_q_threshold_not_biological_same_plant_pollination":True,
            "guild_Lepidoptera_not_identical_to_single_moth_species":True,
            "source_q_SVD_and_stigma_pollen_not_measured":True,
            "no_pollen_limitation_selfing_or_mature_seed_data":True,
            "no_causal_pollination_rescue_fraction_reported":True,
            "no_new_metaregression_or_statistical_test":True,
        }
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    from openpyxl import load_workbook
    b=source_bytes()
    wb=load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    plant=read_sheet(wb,"species visitation")
    visitation=read_sheet(wb,"visitation")
    wb.close()
    z=source_plant_q_thresholds(plant,visitation,strict=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(z,ensure_ascii=False,indent=2)+"\n")
    print("LU2026_ORIGINAL_SINGLE_FIELD_STUDY_COMMUNITY_CONDITIONAL_Q_RATIO",z["community_conditional_qfly_qmoth_threshold_if_all_visits_act_on_same_equivalent_plants"])
    print("LU2026_ORIGINAL_6_PLANT_REQUIRED_Q_DISTRIBUTION",z["focal_six_plant_required_efficiency_ratio_descriptive_summary"])
    print("LU2026_NO_ACTUAL_POLLINATION_COMPENSATION",z["guard"])

if __name__=="__main__":
    main()
