#!/usr/bin/env python3
"""Audit and import Lu 2026 ORIGINAL-PAPER reported access-specific ratios.

IMPORTANT CROSS-REPO CORRECTION from user's IWE source audit:
Lu JPE 2026's -56% SEED/OVULE response applies to NIGHT-ONLY pollination,
NOT open day-and-night pollination.  The observed pollen responses -23.09%
open and -53.80% night-only identify a descriptive ratio of response ratios
but not a fitted intervention interaction with a standard error.

Reads the provenance-pinned user-owned IWE source-transcription registry
directly. Does NOT claim the table is the authors' raw data or that the
publication PDF/appendix was downloaded. Original accepted manuscript was
inspected in IWE by cached full-text (33 pp, Results p14 and Discussion p16).
"""

from __future__ import annotations
import csv
import hashlib
import io
import json
import math
import pathlib
import urllib.request

DOI="10.1093/jpe/rtag160"
SOURCE_REPOSITORY="zuizui0223/iwe"
SOURCE_COMMIT="cd0130b4730735f0ec90b1162780ad67b19a222a"
SOURCE_PATH="data/registry/lu2026_access_scoped_contrasts_v01.csv"
SOURCE_SHA256="b8f18e68f3efcbdc108b59fcd79df9551fd6651cea6a34a0f6ae9fdf5d7b6b43"
SOURCE_URL=f"https://raw.githubusercontent.com/{SOURCE_REPOSITORY}/{SOURCE_COMMIT}/{SOURCE_PATH}"
AUTHORS_PDF="https://academic.oup.com/jpe/advance-article-pdf/doi/10.1093/jpe/rtag160/68754751/rtag160.pdf"

PLANTS=("Pedicularis kansuensis","Euphrasia regelii")
ENDPOINTS=("conspecific_pollen","seed_set")
ACCESSES=("open","day_only","night_only")
DEFINITIONS={
 "conspecific_pollen":"conspecific_pollen_grains_per_stigma",
 "seed_set":"developed_seeds_divided_by_total_ovules_per_fruit",
}
NUMERIC={
 ("Pedicularis kansuensis","conspecific_pollen","open"):(-23.09,.7691),
 ("Pedicularis kansuensis","conspecific_pollen","night_only"):(-53.8,.4620),
 ("Pedicularis kansuensis","seed_set","night_only"):(-56.0,.4400),
}
FIELDS={"doi","plant","endpoint","pollination_access","reported_effect",
        "percent_change","lit_dark_ratio","endpoint_definition","source_pages",
        "sem","sampling_variance","p_value","raw_data_obtained"}


def get_original() -> list[dict]:
    req=urllib.request.Request(SOURCE_URL,headers={
        "User-Agent":"ALAN-IWE-crossrepo-proof/2026-11-15 (https://github.com/zuizui0223/alan)"
    })
    with urllib.request.urlopen(req,timeout=45) as response:
        raw=response.read(25_000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA256:
        raise ValueError("IWE original source-transcription source SHA mismatch")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    return rows


def validate(rows:list[dict]) -> dict:
    if len(rows)!=12:
        raise ValueError("Original source crosswalk must have exact 2×2×3 combinations")
    if any(not FIELDS.issubset(set(r)) for r in rows):
        raise ValueError("Missing provenance or uncertainty fields")
    keys=[(r["plant"],r["endpoint"],r["pollination_access"]) for r in rows]
    expected={(p,e,a) for p in PLANTS for e in ENDPOINTS for a in ACCESSES}
    if len(set(keys))!=12 or set(keys)!=expected:
        raise ValueError("Original source scope has missing/duplicated access arm")
    numeric={}
    for row in rows:
        k=(row["plant"],row["endpoint"],row["pollination_access"])
        if row["doi"]!=DOI or row["raw_data_obtained"]!="no":
            raise ValueError("DOI mismatch or falsely promoted author raw data")
        if row["endpoint_definition"]!=DEFINITIONS[row["endpoint"]]:
            raise ValueError("Outcome denominator mismatch")
        if not row["source_pages"].strip():
            raise ValueError("Missing primary manuscript page provenance")
        if any(row[field].strip() for field in ("sem","sampling_variance","p_value")):
            raise ValueError("Original text did not publish effect-level precision")
        if k in NUMERIC:
            pct,ratio=NUMERIC[k]
            if row["reported_effect"]!="decrease_reported":
                raise ValueError("Published direction drift")
            if not math.isclose(float(row["percent_change"]),pct,abs_tol=1e-10):
                raise ValueError("Published percentage drift")
            if not math.isclose(float(row["lit_dark_ratio"]),ratio,abs_tol=1e-10):
                raise ValueError("Published light/dark ratio drift")
            numeric[k]=ratio
        else:
            if (row["reported_effect"]!="no_significant_effect_reported"
                or row["percent_change"].strip() or row["lit_dark_ratio"].strip()):
                raise ValueError("Do not set nonsignificant outcomes to RR=1")
            numeric[k]=None

    rr_open=numeric[("Pedicularis kansuensis","conspecific_pollen","open")]
    rr_night=numeric[("Pedicularis kansuensis","conspecific_pollen","night_only")]
    B=math.log(rr_open/rr_night)
    matched=[]
    for plant in PLANTS:
        for endpoint in ENDPOINTS:
            a=numeric[(plant,endpoint,"open")]
            b=numeric[(plant,endpoint,"night_only")]
            matched.append({
                "plant":plant,
                "endpoint":endpoint,
                "open_lit_dark_ratio":a,
                "night_only_lit_dark_ratio":b,
                "ratio_of_response_ratios":None if a is None or b is None else a/b,
                "log_ratio_of_response_ratios":None if a is None or b is None else math.log(a/b),
                "variance_and_CI":None,
                "interpretation":"Pair of rounded author-reported percentages; no fitted interaction or SE."
            })

    if not math.isclose(B,.509656102,abs_tol=1e-6):
        raise ValueError("Published pollen access contrast changed")
    if numeric[("Pedicularis kansuensis","seed_set","open")] is not None:
        raise ValueError("The -56% night seed outcome may not be mislabelled OPEN")

    return {
      "schema":"ALAN_Lu2026_author_text_access_scope_crossrepo_20261115_v1",
      "source_paper_doi":DOI,
      "accepted_manuscript_url":AUTHORS_PDF,
      "author_fulltext_not_acquired_in_this_execution":True,
      "source_from_user_repository":SOURCE_REPOSITORY,
      "source_frozen_commit":SOURCE_COMMIT,
      "source_git_path":SOURCE_PATH,
      "source_transcription_sha256":SOURCE_SHA256,
      "source_is_authored_IWE_transcription_not_raw_trial_data":True,
      "original_studies":1,
      "site_level_matched_lit_dark_pairs_reported":5,
      "plants":list(PLANTS),
      "endpoint_definitions":DEFINITIONS,
      "source_combinations_count":12,
      "source_nonmissing_numeric_ratios":3,
      "nonsignificant_but_unquantified_ratios":9,
      "numeric_open_vs_night_pollen_pairs":1,
      "numeric_open_vs_night_seed_pairs":0,
      "source_rounded_point_pairs":matched,
      "Pedicularis_pollen_open_vs_night": {
         "RR_lit_dark_open":rr_open,
         "RR_lit_dark_night_only":rr_night,
         "ratio_of_response_ratios":rr_open/rr_night,
         "log_ratio_of_response_ratios":B,
         "CI95":None,
         "sampling_covariance":None,
         "causal_daytime_compensation_inferred":False
       },
      "Pedicularis_seed_night_only_lit_dark_ratio":.44,
      "Pedicularis_seed_open_lit_dark_ratio":None,
      "Pedicularis_seed_NIGHT_ONLY_minus_56_percent":True,
      "Pedicularis_seed_undeclared_full_open_decline":True,
      "seed_set_source_denominator":"developed seeds / ovules among formed fruits (not seed count per initial flower)",
      "independent_light_by_access_reproduction_meta_effects":0,
      "global_ALAN_night_dependence_beta_V":None,
      "invalid_legacy_inference":"Do NOT compare the 56% night-only seed decrease with 35.71% night share as if both were open-pollination whole-plant seed outcomes. Different experimental access and denominator.",
      "important_limitation":"Open-arm seed-set reported non-significant is not an observed zero nor proof of perfect compensation. Original Table S4b and arm-wise SE and covariance are absent."
    }


def main():
    out=validate(get_original())
    dest=pathlib.Path("results/alan_Lu2026_IWE_access_scoped_verified_20261115.json")
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
       "rows":out["source_combinations_count"],
       "numeric_ratios":out["source_nonmissing_numeric_ratios"],
       "night_only_seed_rr":out["Pedicularis_seed_night_only_lit_dark_ratio"],
       "open_seed_rr":out["Pedicularis_seed_open_lit_dark_ratio"],
       "pollen_B":out["Pedicularis_pollen_open_vs_night"]["log_ratio_of_response_ratios"],
       "paired_pollen":out["numeric_open_vs_night_pollen_pairs"],
       "paired_seed":out["numeric_open_vs_night_seed_pairs"],
       "beta_V":out["global_ALAN_night_dependence_beta_V"],
    },ensure_ascii=False))


if __name__=="__main__":
    main()
