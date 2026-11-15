#!/usr/bin/env python3
"""ALAN A: hard source-extraction readiness census, NEVER a fabricated meta slope.

The 14-entry study registry has 12 reproductive publications (9 with
direct manipulation, 3 observational gradients) and two mechanistic-only
literature records. A manuscript must NOT confuse these with 9 equivalent
fruit-set log-RRs or 9 independent measures of nocturnal dependence.

Cross-study moderator V needs dark/day-only and dark/open outcomes measured
within the SAME source population and a clustered covariance with the
lit/open reproductive outcome. A night-only seed loss in Lu is not an
unrestricted/total seed loss. No missing outcome/SE is made equal to zero.
"""
from __future__ import annotations
import csv
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
REG=ROOT/"config/alan_A_reproductive_screen_20261115.json"
OUT=ROOT/"results/alan_A_endpoint_independence_readiness_20261115.json"
CSV=ROOT/"results/alan_A_endpoint_independence_readiness_20261115.csv"

# Manual source-adjudicated rules, frozen in the same repository and
# checked against the independently sourced 14-paper registry.
SOURCE_GATE={
 "knop2017":("direct","achenes/flowerhead, all-visitor bagging","no","no","no","5 lit+5 dark meadows, 100 plants, 200 capitula; light x ALL-visitor bagging ≠ nighttime access"),
 "macgregor2019":("direct","at least one seed/original flower","yes","no","no","5 plots x 4 rotations, dark-specific access×lighting cells unavailable; pooled diel exclusion not V_dark"),
 "boom2020":("direct","fruit/fertilized flower and model-estimated seeds/flower","no","no","no","one spectrum experiment, shared dark controls; published seeds/flower model ANOVA p=.63"),
 "giavi2020":("direct","logistic odds of seed/capsule and nursery infestation","no","no","no","4 lit+4 dark meadows; adjacent dark treated within lit meadows, not independent controls"),
 "bariles2021":("observational","fruit/seed and pollen","no","no","no","six-site skyglow/radiance+tree cover observation, no assigned experimental light"),
 "wilson2021":("direct","mature fruits/plant","no","no","no","5 field sites × 15 randomized transects 8 lit/7 dark; skyglow gradient separate; daytime yucca moth"),
 "dzul2022":("observational","fruit, seed, stigma pollen","no","no","no","light/noise gradient, no assigned lamps"),
 "li2026":("direct","seed-positive plant, fruits and seeds/starting bud","no","no","no","3 fields and 3 moth-syndrome plants; original seed data available, but no matched night/day exclusion"),
 "lu2026jpe":("direct","seed/ovule among fruit (NIGHT ONLY: -56%), pollen open/night matched","yes","no","no","5 paired lit/dark plots, 2 plants; open seed ratio NOT recovered, original 16 cells+SE missing"),
 "xiao2025":("direct","fruiting inflorescence count, reproduction phenology","no","no","no","not original flower successful fruit rate; no measured nighttime insect access dependence"),
 "grognuz2026pre":("direct","target pot/inflorescence seed count after DAYTIME 90-min visits","no","no","no","28 daylight choice blocks, repeated pots under shared lamps; pre-exposure night flowers bagged"),
 "soteras2025":("observational","fruit set under urbanization gradient","no","no","no","observational urban light/impervious cover, no assigned ALAN treatment"),
 "lu2026evo":("not_reproductive","day/night insect visits ONLY","no","no","no","do not double-count Lu JPE reproduction; 2026 Ecol Evol is visit-network source"),
 "grognuz2026pheno":("not_reproductive","flowering time/phenology ONLY","no","no","no","no measured fruit or mature seed output; same author group does not imply same experiment")
}

def derive(reg:dict)->dict:
    entries=reg["entries"]
    if len(entries)!=14:
        raise AssertionError("Source registry identity must be explicitly reviewed before count change")
    ids=[r["id"] for r in entries]
    if len(ids)!=len(set(ids)) or set(ids)!=set(SOURCE_GATE):
        raise AssertionError("Mismatched or duplicated source IDs; manually update proof gate")
    out=[]
    for e in entries:
        kind,endpoint,factorial,V,meta,notes=SOURCE_GATE[e["id"]]
        if kind=="direct" and not e["light_design"].startswith("manipulated"):
            raise AssertionError("Direct light label changed without source review "+e["id"])
        if kind=="observational" and e["light_design"]!="observational":
            raise AssertionError("Observational light gradient was promoted "+e["id"])
        if kind=="not_reproductive" and e["a_primary_stage"]!="exclude_A_no_reproduction":
            raise AssertionError("A mechanism-only paper may not add reproductive effects")
        if e["verified_joint_effect_moderator_se"]:
            raise AssertionError("Source registry claims a completed V×y covariance unsupported "+e["id"])
        if e["id"]=="lu2026jpe":
            src=e["source_audit"]
            if src["Pedicularis_seed_lit_dark_night_only"]!=.44 or src["Pedicularis_seed_lit_dark_open"] is not None:
                raise AssertionError("Lu -56% was wrongly promoted to full open reproduction")
            if src["seed_endpoint"]!="developed seeds / total ovules per formed fruit":
                raise AssertionError("Lu seed denominator changed")
        if e["id"]=="macgregor2019" and e["external_night_dependence"] not in ("same_study_exclosure_not_external","dark_only_access_cells_missing"):
            raise AssertionError("Macgregor dark-only V is not recovered")
        out.append({
          "source_id":e["id"],"doi":e["doi"],"year":e["year"],
          "design_class":kind,"outcome_grain":endpoint,
          "same_study_light_by_access_factorial":factorial,
          "dark_specific_V_and_SE_recovered":V,
          "qualified_primary_repro_lnRR_and_cluster_SE":"no",
          "joint_covariance_V_and_lnRR_recovered":"no",
          "confirmed_reproductive_meta_beta_V_ready":meta,
          "author_or_original_data_status":e["screen_status"],
          "reason":notes
        })
    kinds=Counter(x["design_class"] for x in out)
    if kinds!={"direct":9,"observational":3,"not_reproductive":2}:
        raise AssertionError("Source candidate breakdown changed; explicit review required")
    factorial=sum(x["same_study_light_by_access_factorial"]=="yes" for x in out)
    if factorial!=2:
        raise AssertionError("ALAN x access original factorial source paper count changed")
    blockers=Counter()
    for r in out:
        if r["design_class"]!="direct":continue
        if r["qualified_primary_repro_lnRR_and_cluster_SE"]=="no":blockers["unusable_primary_outcome_and_precision"]+=1
        if r["dark_specific_V_and_SE_recovered"]=="no":blockers["dark_only_moderator_missing"]+=1
        if r["joint_covariance_V_and_lnRR_recovered"]=="no":blockers["full_joint_covariance_missing"]+=1
    return {
      "schema":"ALAN_A_20261115_strict_effect_meta_readiness_v1",
      "source_base":"zuizui0223/alan config/alan_A_reproductive_screen_20261115.json, source-checked by user IWE original accepted manuscript audit",
      "candidate_publications_total":len(out),
      "candidate_reproduction_papers":kinds["direct"]+kinds["observational"],
      "confirmed_direct_lighting_reproduction_publications":kinds["direct"],
      "observational_gradient_reproduction_publications":kinds["observational"],
      "known_mechanism_only_papers":kinds["not_reproductive"],
      "known_light_by_daynight_access_factorial_papers":factorial,
      "open_reproduction_lnRR_with_clustered_SE_and_independent_V":0,
      "global_meta_beta_V":None,
      "no_nuisance_seed_OR_to_lnRR_conversion":True,
      "key_previously_invalid_inference":"Lu 56% seed/ovule decline is NIGHT-ONLY; open unrestricted seed change numerically unknown. Do not treat it as pooled all-channel ALAN reproductive loss.",
      "additional_direct_papers_found_by_this_step":0,
      "actual_new_author_raw_data_files_acquired_by_this_step":0,
      "blockers_by_direct_study":dict(blockers),
      "source_statuses":out
    }

def main():
    registry=json.loads(REG.read_text(encoding="utf-8"))
    summary=derive(registry)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    with CSV.open("w",newline="",encoding="utf-8") as out:
        columns=list(summary["source_statuses"][0])
        writer=csv.DictWriter(out,fieldnames=columns)
        writer.writeheader();writer.writerows(summary["source_statuses"])
    print(json.dumps({k:summary[k] for k in (
        "candidate_publications_total","candidate_reproduction_papers",
        "confirmed_direct_lighting_reproduction_publications",
        "known_light_by_daynight_access_factorial_papers",
        "open_reproduction_lnRR_with_clustered_SE_and_independent_V",
        "global_meta_beta_V","blockers_by_direct_study"
    )},ensure_ascii=False))

if __name__=="__main__":
    main()
