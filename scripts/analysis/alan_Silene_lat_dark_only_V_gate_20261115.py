#!/usr/bin/env python3
"""Source-faithful S. latifolia diel-night dependence study-pooling audit.

Existing alan results/night_exclusion/experiments.csv contains both natural
day/open data and one Macgregor2019 pooled-across-light-treatment contrast.
For a 'dark-only' nocturnal dependence moderator, the pooled experimental
light data must be omitted. This code NEVER edits the original results.
It emits an explicit source-assignment sensitivity panel.

Source counts and loss values are exact original catalog source extracts,
not source-raw individual plant datasets; source SEMs unavailable.
"""

import collections
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
ORIG=ROOT/"results/night_exclusion/experiments.csv"
SPEC=ROOT/"results/night_exclusion/species.csv"
SOURCE_MIXED={"Macgregor_2019"}
STRICT_OUT=ROOT/"results/night_exclusion/S_latifolia_dark_only_source_gate_20261115.json"
STRICT_CSV=ROOT/"results/night_exclusion/strict_night_V_by_species_20261115.csv"

def load_rows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def study_means(rows):
    groups=collections.defaultdict(list)
    for r in rows:
        if r["used"]!="True":continue
        loss=float(r["loss"])
        if not math.isfinite(loss):
            raise ValueError("Nonfinite experiment-level night dependence")
        denom=float(r["open_value"])
        day=float(r["day_value"])
        if denom<=0 or not math.isclose(loss,1-day/denom,abs_tol=1e-9):
            raise AssertionError("Original extraction D/O not source-consistent")
        groups[(r["plant"],r["study_ID"])].append(loss)
    out=collections.defaultdict(list)
    details=[]
    for (species,study),values in sorted(groups.items()):
        mean=sum(values)/len(values)
        out[species].append((study,mean))
        details.append({"taxon":species,"study_ID":study,"n_contrasts":len(values),
                        "mean_within_study_V":mean,
                        "pooled_across_light_treatments":study in SOURCE_MIXED,
                        "is_ALAN_light_free_V_source":study not in SOURCE_MIXED})
    return out,details

def equally_weight_studies(items):
    return None if not items else sum(x[1] for x in items)/len(items)

def compute(rows,species_table):
    if len(rows)!=42:
        raise AssertionError("Unexpected source version; re-audit full 42 rows")
    groups,details=study_means(rows)
    if len(groups)!=12:
        raise AssertionError("Study/plant mapping changed")
    spp={r["plant"]:r for r in species_table}
    if set(groups)!=set(spp):
        raise AssertionError("Original reported species table mismatch")
    results=[]
    for taxon,items in sorted(groups.items()):
        old=equally_weight_studies(items)
        strict=[x for x in items if x[0] not in SOURCE_MIXED]
        fixed=equally_weight_studies(strict)
        if not math.isclose(old,float(spp[taxon]["loss"]),abs_tol=1e-10):
            raise AssertionError("Cannot replicate original registered species mean: "+taxon)
        values=[r for r in rows if r["used"]=="True" and r["plant"]==taxon
                and r["study_ID"] not in SOURCE_MIXED]
        source_values=[float(x["loss"]) for x in values]
        results.append({
            "plant":taxon,
            "original_V_with_all_used_studies":old,
            "strict_no_assigned_lights_V":fixed,
            "difference_strict_minus_original":fixed-old,
            "n_studies_original":len(items),
            "n_studies_strict":len(strict),
            "n_source_contrasts_strict":len(values),
            "strict_source_contrast_min":min(source_values) if source_values else None,
            "strict_source_contrast_max":max(source_values) if source_values else None,
            "source_mean_is_not_same_population_V_for_ALAN_experiment":True,
            "no_unreported_variance_imputed":True
        })
    s=next(x for x in results if x["plant"]=="Silene latifolia")
    detail=[d for d in details if d["taxon"]=="Silene latifolia"]
    if not (math.isclose(s["original_V_with_all_used_studies"],.6718269901631905,abs_tol=1e-10)
            and math.isclose(s["strict_no_assigned_lights_V"],.7236903716106219,abs_tol=1e-10)
            and s["n_studies_strict"]==4):
        raise AssertionError("Critical source correction not reproduced")
    # Descriptive species-level OLS only; this is NOT a refit of the prior
    # 78-region cluster-robust observational spatial model.
    gaps={r["plant"]:float(r["gap"]) for r in species_table}
    def descriptive_slope(field):
        x=[r[field] for r in results]
        y=[gaps[r["plant"]] for r in results]
        mx=sum(x)/len(x); my=sum(y)/len(y)
        return sum((a-mx)*(b-my) for a,b in zip(x,y))/sum((a-mx)**2 for a in x)
    original_gap_slope=descriptive_slope("original_V_with_all_used_studies")
    strict_gap_slope=descriptive_slope("strict_no_assigned_lights_V")
    return {
        "schema":"ALAN_B_Silene_latifolia_dark_only_gate_v1",
        "source":"results/night_exclusion/experiments.csv (42 rows); results/night_exclusion/species.csv",
        "stage":"source-level moderator definition correction, not a new original study",
        "original_rows":len(rows),
        "species_n":len(results),
        "source_confounded_with_light":["Macgregor_2019"],
        "source_pooled_lighting_note":"Macgregor2019 day/open from Figure2 is pooled across dark,HPS full/part-night,LED full/part-night, so using it as DARK-only V is forbidden.",
        "s_lat_original_V":s["original_V_with_all_used_studies"],
        "s_lat_strict_V":s["strict_no_assigned_lights_V"],
        "s_lat_absolute_change":s["difference_strict_minus_original"],
        "s_lat_individual_contrasts_without_mac_range":[s["strict_source_contrast_min"],s["strict_source_contrast_max"]],
        "strict_is_not_independently_measured_local_ALAN_V":True,
        "original_species_geographical_model_not_refitted":True,
        "external_species_mean_may_not_transport_to_other_regions":True,
        "new_reproductive_ALAN_effects":0,
        "global_beta_V":None,
        "exploratory_12species_gap_OLS_per_10pp_V_before":original_gap_slope*0.1,
        "exploratory_12species_gap_OLS_per_10pp_V_strict":strict_gap_slope*0.1,
        "full_78region_model_refitted":False,
        "exploratory_OLS_does_NOT_have_cluster_SE_or_p":True,
        "studies":details,
        "species":results
    }

def main():
    r=compute(load_rows(ORIG),load_rows(SPEC))
    STRICT_OUT.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    with STRICT_CSV.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(r["species"][0]))
        writer.writeheader();writer.writerows(r["species"])
    print(json.dumps({k:r[k] for k in
       ("species_n","s_lat_original_V","s_lat_strict_V",
        "s_lat_absolute_change","s_lat_individual_contrasts_without_mac_range",
        "global_beta_V","original_species_geographical_model_not_refitted",
        "exploratory_12species_gap_OLS_per_10pp_V_before",
        "exploratory_12species_gap_OLS_per_10pp_V_strict")},ensure_ascii=False))

if __name__=="__main__":
    main()
