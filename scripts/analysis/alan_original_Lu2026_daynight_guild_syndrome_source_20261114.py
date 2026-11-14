#!/usr/bin/env python3
"""Source-locked ALAN field experiment: real day/night visitor-guild reshuffling.

Lu et al. Ecology and Evolution 2026, DOI 10.1002/ece3.74368. Raw original
Date.xlsx is inside Europe PMC supplementaryFiles outer ZIP and the paper's
ECE3-16-e74368-s001.zip inner ZIP. Original workbook's sha256 preverified.
Primary analysis:
  (i) genuine paired experiment at block×round×diel scale (not pseudo-matched repeat),
  (ii) exact 2^5 block-sign randomization on mean visit differences,
  (iii) weighted raw visitation redistribution, and
  (iv) plant-specific positive documentation overlap (no false plant absence).
Actual plant visits do NOT identify effective pollen deposition, fruit set,
moth-evolved floral syndrome, predator behavior, island shift or phenology.
"""
from __future__ import annotations
import argparse, collections, hashlib, io, itertools, json, math, urllib.request, zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

OUTER="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles"
OUTER_SHA="7a6a4d9cea40cedd87a2122df6f3a20b1303926c5841445e376021b3e377c635"
INNER="ECE3-16-e74368-s001.zip"
XLSX="Date.xlsx"
XLSX_SHA="ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f"
CORE=["time","sampling_time","treatment","block","repeat"]
GUILDS={"Diptera":"D_abundance","Hymenoptera":"H_abundance",
        "Lepidoptera":"L_abundance","Coleoptera":"C_abundance","Others":"O_abundance"}
SOURCE_GUILDS=list(GUILDS)

def require(cond,msg):
    if not bool(cond): raise ValueError("ALAN_LU_ORIGINAL: "+msg)

def sha(data):return hashlib.sha256(data).hexdigest()

def download():
    req=urllib.request.Request(OUTER,headers={"User-Agent":"ALAN-original-Lu-ECE-2026-licensed-open-data"})
    with urllib.request.urlopen(req,timeout=75) as r:
        b=r.read(20_000_000)
    # EuropePMC dynamically assembles the outer supplementary archive; its ZIP
    # receipt may differ byte-wise on repeated requests even with identical
    # published files. Pin the original nested XLSX SHA and ZIP member name.
    require(zipfile.is_zipfile(io.BytesIO(b)), "EuropePMC supplied non-ZIP bytes")
    with zipfile.ZipFile(io.BytesIO(b)) as f:
        require(INNER in f.namelist(),"Original source ECE nested zip absent")
        nested=f.read(INNER)
    with zipfile.ZipFile(io.BytesIO(nested)) as z:
        require(XLSX in z.namelist(),"Original Date.xlsx workbook absent")
        raw=z.read(XLSX)
    require(sha(raw)==XLSX_SHA,"Original published Date.xlsx changed")
    return raw,sha(b)

def sheet(w,name,maxrows=2500):
    require(name in w.sheetnames,"Original author worksheet unavailable: "+name)
    s=w[name]
    # Max_row incorrectly reports 1,048,576 due to blank formatted spreadsheet
    # cells in a few sheets. Strictly bounded source extraction.
    rows=list(s.iter_rows(min_row=1,max_row=min(s.max_row,maxrows),
                          values_only=True))
    header=[str(x).strip() if x is not None else "" for x in rows[0]]
    require(len(set(h for h in header if h))==len([h for h in header if h]),
            "Duplicate original worksheet fields")
    out=[]
    for a in rows[1:]:
        if a[0] is None:continue
        out.append({k:v for k,v in zip(header,a) if k})
    require(len(out)<maxrows,"Potential source worksheet truncation")
    return pd.DataFrame(out)

def exact_paired_block_sign(values):
    """Five independent experimental blocks: sign-flip 32 possible assignments.

    This is a transparent *descriptive randomization sensitivity*; original
    experimental randomization was at plot/block design, not fully documented
    symmetric sign-flip, thus not a formal field-study exact causal test.
    """
    v=np.asarray(values,dtype=float)
    require(len(v)==5 and np.isfinite(v).all(),"Need all 5 actual original block contrasts")
    obs=abs(v.mean())
    null=np.array([abs(np.mean(v*np.array(s)))
                   for s in itertools.product((-1,1),repeat=5)])
    return {"block_mean_difference":float(v.mean()),
            "five_block_values_all_same_sign":bool(np.all(v>0) or np.all(v<0)),
            "two_sided_complete_32_signflip_p_DESCRIPTIVE":float(
                np.mean(null>=obs-1e-12)),
            "not_claimed_exact_original_assignment_test":True}

def original_visitation_design(a):
    require(len(a)==400,"400 original visit survey design rows expected")
    require(set(CORE+["abundance"]+list(GUILDS.values()))<=set(a),
            "Original measured visit columns missing")
    require(set(a.sampling_time)=={"Diurnal","Nocturnal"} and
            set(a.treatment)=={"Lit","Unlit"},
            "Original source field night/light labels changed")
    require(set(a.time.astype(int))=={1,2,3,4} and
            set(a.block.astype(int))=={1,2,3,4,5} and
            set(a.repeat.astype(int))=={1,2,3,4,5},
            "Original 4 rounds, 5 blocks, 5 replicates changed")
    # Source repeat labels are not asserted unique: check block×round means
    # rather than invent a session-level lit/dark match from reused labels.
    for c in ["abundance"]+list(GUILDS.values()):
        a[c]=pd.to_numeric(a[c],errors="raise")
        require(a[c].ge(0).all() and np.isfinite(a[c]).all(),
                "Visit counts cannot be missing or negative")
    sums=a[list(GUILDS.values())].sum(axis=1)
    # Source visit columns do not always conserve to total 'abundance'.
    # Preserve a transparent discrepancy; no false exact guild mass balance.
    a["source_sum_minus_total"]=sums-a.abundance
    return a

def original_design_model(a):
    a=original_visitation_design(a.copy())
    cols={"Total":"abundance",**GUILDS}
    b=a.set_index(CORE)
    outcome={}
    for period in ("Nocturnal","Diurnal"):
        source=a[a.sampling_time.eq(period)]
        require(len(source)==200 and
                (source.groupby("treatment").size()==100).all(),
                "Original within-period trial count drift")
        each={}
        for label,col in cols.items():
            # Original source duplicates some repeat labels. Equal-weight
            # within 20 block×round units, then 5 source blocks independent.
            unit=source.groupby(["time","block","treatment"],observed=True)[col].mean().unstack("treatment")
            require(unit.shape==(20,2) and unit.notna().all().all(),
                    "Original 20 block-round illuminated/dark comparisons not matched")
            blockdelta=(unit["Lit"]-unit["Unlit"]).groupby(level="block").mean()
            require(len(blockdelta)==5,"Original block factor missing")
            by=source.groupby("treatment")[col].agg(["sum","mean"])
            dark=float(by.loc["Unlit","sum"]); lit=float(by.loc["Lit","sum"])
            difference=lit-dark
            each[label]={
                "original_dark_total_visit_events":dark,
                "original_lit_total_visit_events":lit,
                "net_difference_original_100_surveys_each":difference,
                "lit_dark_visit_ratio":float(lit/dark) if dark>0 else None,
                "relative_change_percent":float(100*difference/dark) if dark>0 else None,
                "block_level_difference":exact_paired_block_sign(blockdelta.to_numpy()),
                "n_paired_block_rounds":int(len(unit)),
                "source_repeat_labels_not_forced_to_be_unique":True}
        total_d=each["Total"]["net_difference_original_100_surveys_each"]
        group_d=sum(each[x]["net_difference_original_100_surveys_each"]
                    for x in SOURCE_GUILDS)
        source_errors=source["source_sum_minus_total"].to_numpy(float)
        conservation={
            "source_row_sums_match_total_n":int(np.sum(np.isclose(source_errors,0,atol=1e-8))),
            "source_original_200_rows":len(source),
            "original_group_sum_vs_total_delta":float(group_d-total_d),
            "sum_discrepancies_min_max":[float(np.min(source_errors)),float(np.max(source_errors))],
            "sum_discrepancies_total":float(np.sum(source_errors)),
            "group_composition_complete":bool(np.allclose(source_errors,0,atol=1e-8))
        }
        d=each["Diptera"]["net_difference_original_100_surveys_each"]
        m=each["Lepidoptera"]["net_difference_original_100_surveys_each"]
        # Signed visit ratio describes observed replacement OPPORTUNITY only.
        replacement=(d/-m) if d>0 and m<0 else None
        outcome[period]={"n_visitation_source_records":len(source),
            "independent_field_blocks":5,"survey_rounds":4,
            "nominal_replicate_levels":5,
            "actual_block_round_light_record_count_range":[int(source.groupby(["time","block","treatment"],observed=True).size().min()),int(source.groupby(["time","block","treatment"],observed=True).size().max())],
            "source_repeated_full_design_keys":int(source.duplicated(CORE).sum()),
            "experimental_survey_pair_group_effects":each,
            "fly_visit_gain_divided_by_moth_visit_loss_if_opposite_signs":replacement,
            "original_source_visit_count_conservation_audit":conservation,
            "actual_pollen_replacement_measured":False,
            "per_flower_visitation_denom_NOT_available":True}
    return outcome

def species_descriptive(s):
    expected=set(CORE+["plant","all","Lepidoptera","Coleoptera",
                       "Diptera","Hymenoptera","Others"])
    require(expected<=set(s),"Original plant identity × insect-order counts missing")
    require(len(s)==1038,"Original species-visitation source row count changed")
    for c in ("all","Lepidoptera","Coleoptera","Diptera","Hymenoptera","Others"):
        s[c]=pd.to_numeric(s[c],errors="raise")
        require(s[c].ge(0).all(),"Original source species visit count invalid")
    source_species_error=(s[["Lepidoptera","Coleoptera","Diptera","Hymenoptera","Others"]].sum(axis=1)-s["all"])
    species_conserved=int(np.isclose(source_species_error,0,atol=1e-8).sum())
    # Not all plant species necessarily flower at every survey; unlisted rows
    # remain UNKNOWN and will never be assigned real flowering availability 0.
    total=s.plant.nunique()
    support={"all_original_plant_positive_observation_rows":len(s),
      "plant_names_with_any_recorded_visits":int(total),
      "species_rows_with_five_order_sum_matching_all":species_conserved,
      "species_guild_sum_to_total_discrepancy_range":[float(source_species_error.min()),float(source_species_error.max())],
      "species_counts_are_not_assumed_complete_partition":species_conserved!=len(s),
      "no_missing_plant_survey_row_coded_zero":True,
      "flower_availability_by_plant_survey_unmeasured":True}
    groups={}
    for per in ("Nocturnal","Diurnal"):
        q=s[s.sampling_time.eq(per)].copy()
        by=q.groupby(["plant","treatment"],observed=True).agg(
            raw_visits=("all","sum"),
            moth=("Lepidoptera","sum"),
            fly=("Diptera","sum")).reset_index()
        if not {"Lit","Unlit"} <= set(by.treatment):
            groups[per]={"status":"NOT_ESTIMABLE_NO_POSITIVE_SOURCE_COMMON_SUPPORT"}
            continue
        w=by.pivot(index="plant",columns="treatment",
                   values=["raw_visits","moth","fly"])
        supported=w["raw_visits"].dropna()
        common=supported.index
        q2=w.loc[common].copy()
        if len(q2):
            dark=q2[("raw_visits","Unlit")].astype(float)
            moth=q2[("moth","Unlit")].astype(float)
            fly=q2[("fly","Unlit")].astype(float)
            lit=q2[("raw_visits","Lit")].astype(float)
            majority=moth.gt(dark/2)
            moth_down=q2[("moth","Lit")].astype(float).lt(moth)
            fly_up=q2[("fly","Lit")].astype(float).gt(fly)
            # Source-documented visitor proportions, not a plant evolutionary
            # pollination syndrome or a causal service-loss coefficient.
            def subgroup(label,mask):
                mask=np.asarray(mask,dtype=bool)
                z=q2.loc[mask]
                vals={}
                for name in ("raw_visits","moth","fly"):
                    dark_count=float(z[(name,"Unlit")].sum())
                    lit_count=float(z[(name,"Lit")].sum())
                    vals[name]={
                        "unlit_source_positive_visit_counts":dark_count,
                        "lit_source_positive_visit_counts":lit_count,
                        "lit_minus_unlit_visit_counts":lit_count-dark_count,
                        "lit_over_unlit_visit_ratio":float(lit_count/dark_count) if dark_count>0 else None,
                    }
                return {"source_control_composition_group":label,
                    "plant_taxa_with_both_lit_unlit_positive_visit_records":int(len(z)),
                    "treatment_visit_counts_UNADJUSTED_FOR_TRUE_FLOWER_AVAILABILITY":vals,
                    "not_causal_effective_pollen_or_seed":True}
            groups[per]={
              "realized_control_moth_visitation_majority_groups":{
                "source_unlit_moth_fraction_above_half":subgroup("moth_over_half",majority),
                "source_unlit_moth_fraction_not_over_half":subgroup("moth_not_over_half",~majority)},
              "unique_plant_names_with_source_visits":int(q.plant.nunique()),
              "plants_positive_in_both_original_lit_and_dark_observations":int(len(q2)),
              "dark_moth_dominant_documented_plants":int(majority.sum()),
              "dark_moth_dominant_with_lit_moth_visits_lower":int(
                  (majority&moth_down).sum()),
              "both_lit_moth_lower_and_lit_fly_higher_documented_plants":int(
                  (moth_down&fly_up).sum()),
              "plants_moth_visits_lower":int(moth_down.sum()),
              "plants_fly_visits_higher":int(fly_up.sum()),
              "all_plant_level_population_effort_and_flower_availability_not_equal":True,
              "all_plant_guild_labels_descriptive_NOT_effective_pollination":True
            }
        else:
            groups[per]={"status":"NOT_ESTIMABLE_NO_POSITIVE_SOURCE_COMMON_SUPPORT"}
    return {"source":support,"by_period":groups}

def main(out):
    out=Path(out)
    require(not out.exists(),"Original source aggregate exists; refuse overwrite")
    raw,outer_actual_sha=download()
    w=load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
    vis=sheet(w,"visitation")
    sp=sheet(w,"species visitation")
    total=original_design_model(vis)
    species=species_descriptive(sp)
    w.close()
    result={
       "schema":"ALAN_original_Lu2026_manipulated_field_syndrome_turnover_daynight_v1",
       "status":"AUTHENTIC_OPEN_PUBLISHED_FIELD_VISIT_COUNTS_BLOCK_PAIRED_NO_SEED_CLAIM",
       "source":{"article_DOI":"10.1002/ece3.74368",
         "PMCID":"PMC13645961","original_workbook_title":"Date.xlsx",
         "first_EuropePMC_assembled_outer_zip_receipt_SHA256":OUTER_SHA,
         "actual_outer_zip_dynamic_receipt_SHA256":outer_actual_sha,
         "inner_original_workbook_SHA_pinned_instead_of_dynamic_outer":True,
         "actual_workbook_sha256":XLSX_SHA,
         "verified_original_400_experiment_survey_rows":len(vis),
         "verified_original_1038_plant_visitation_source_rows":len(sp),
         "source_independent_field_blocks":5,
         "original_light_treatments":["Lit","Unlit"],
         "original_day_night":["Diurnal","Nocturnal"],
         "original_four_survey_rounds":4,
         "original_replicates_per_block_round":5},
       "guild_visitation":total,
       "documented_plant_visit_opportunities":species,
       "interpretation_limits":{
          "real_ALAN_manipulation":True,
          "one_original_field_experiment_not_100_independent_site_studies":True,
          "guild_replacement_only_VISITS_not_effective_pollination":True,
          "source_label_LEPIDOPTERA_includes_moths_and_butterflies":True,
          "no_local_island_mainland_contrast":True,
          "not_a_predator_attack_experiment":True,
          "source_missing_plant_visit_not_plant_absence_or_first_flowering":True,
          "no_direct_stigma_pollen_seed_fruit_or_selfing_measured":True,
          "no_plant_trait_verified_evolutionary_syndrome":True,
          "no_raw_plant_species_x_plot_id_site_rows_exported":True,
          "five_block_signflip_not_strict_permutation_randomized_trial":True,
       }}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ALAN_REAL_LU2026_DAYNIGHT_GUILD_VISITS",
          json.dumps({"status":result["status"],"source":result["source"],
                      "guild_visitation":result["guild_visitation"],
                      "plant":result["documented_plant_visit_opportunities"]},
                      ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    main(a.out)
