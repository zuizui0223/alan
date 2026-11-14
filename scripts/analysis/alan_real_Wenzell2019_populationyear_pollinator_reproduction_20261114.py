#!/usr/bin/env python3
"""Reproduce a real study population-year visitor x reproduction join.

Bridge the source-authenticated 301-case Kendall/Nicholson global diel
fruit/seed reproductive meta-dataset to independent, public plant-level
Castilleja exclusion and same-population/year actual visitor counts from
Wenzell et al. 2024 (Ecology and Evolution, DOI 10.1002/ece3.11443).

Wenzell plant exclusion data are a source of seven 2019 population-year units:
three nocturnal-dominant and four diurnal-dominant. Exactly SIX of the seven
can be linked to the previously SHA-authenticated global 301-case study,
five within 1km and one 1-2km (provisional). The seventh is a genuinely
independently available original population that is NOT inserted as an extra
global meta-analysis case.

This is one independent SOURCE PAPER, NOT six/seven independent studies;
visitation is during the same population-year, not necessarily the same
individual plant, specific exclusion day, or evidence of stigma pollen
deposition. NO original geographic or species-specific source rows exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ORIGINAL_SOURCE_SHA256 = "21b1fef1beb9e6238cfa47101b991815d4c0f33787b9737f834426dc3bab15af"
ORIGINAL_AUTHOR_BLOB = "39a82d966afcdde670603add1c881ace8791ffe9"
SOURCE_STUDY = "Wenzell_2021"
DIRECT_SOURCE_DOI = "10.1002/ece3.11443"
FILES = {
    "direct_raw": "85491786cd9381120bea598b9a9d73ae08acd5ad",
    "direct_population": "6a620db6f3ffb4425ff48b114056ecf9907dbb5e",
    "direct_visitors": "050f2db829a7ce06c7007d85fbb8c4c1010122d1",
}
SPECIES_CODES = {
    "S": "Castilleja sessiliflora",
    "C": "Castilleja citrina",
    "L": "Castilleja lindheimeri",
}
VISITOR_COUNTS = [
    "ct_otherD", "ct_hawkmoth", "ct_small_bee",
    "ct_otherN", "ct_bumblebee", "ct_hummingbird",
]


def require(ok, msg):
    if not bool(ok):
        raise ValueError("ALAN_DIRECT_POPULATION_YEAR_GUILD: " + str(msg))


def gitblob(path):
    b = Path(path).read_bytes()
    return hashlib.sha1(b"blob " + str(len(b)).encode() + b"\0" + b).hexdigest()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def haversine_km(lat1, lon1, lat2, lon2):
    r = np.pi / 180.
    dlat = (float(lat2) - float(lat1)) * r
    dlon = (float(lon2) - float(lon1)) * r
    q = (np.sin(dlat / 2)**2 +
         np.cos(float(lat1) * r) * np.cos(float(lat2) * r) *
         np.sin(dlon / 2)**2)
    return float(2*6371.0088*np.arctan2(np.sqrt(q), np.sqrt(max(0, 1-q))))


def original_direct_population(raw, pop, visitors):
    require(set(["plantID", "pop", "year", "trt_Open", "nFlrs", "nFrts"])
            <= set(raw.columns), "Missing direct individual fruit-set fields")
    require(set(["pop", "year", "sp", "majorPollPeriod", "latN", "lonW",
                 "visTotal"]) <= set(pop.columns), "Missing public population metadata")
    require(set(["pop", "year"]+VISITOR_COUNTS) <= set(visitors.columns),
            "Missing public real visitor group event count fields")
    require(len(raw) == 525 and len(pop) == 10 and len(visitors) == 17,
            "Independent original direct author release rows altered")
    sel = pop.loc[pop.year.eq(2019)].copy()
    require(sel["pop"].nunique() == 7 and len(sel) == 7,
            "Expected exactly seven original 2019 direct field populations")
    require(sel.sp.isin(SPECIES_CODES).all() and
            sel.majorPollPeriod.isin(["diurnal", "nocturnal"]).all(),
            "Missing original visitor group or plant mapping")
    result = []
    used_plant_ids = set()
    for row in sel.itertuples(index=False):
        group = raw.loc[raw.year.eq(2019) & raw["pop"].eq(row.pop)].copy()
        require(len(group) > 0 and group.plantID.is_unique,
                "Source population-year individual plant data missing/duplicated")
        require(not set(group.plantID) & used_plant_ids,
                "Original direct individual plants repeated across populations")
        used_plant_ids.update(group.plantID)
        require(group.trt_Open.isin(["day", "night", "open"]).all(),
                "Unexpected 2019 reproductive original intervention type")
        require(set(group.trt_Open) == {"day", "night", "open"},
                "Actual day/night/open outcomes missing for a population-year")
        require(group.majorPollPeriod.eq(row.majorPollPeriod).all(),
                "Original plant-level source period label disagrees with field metadata")
        outcomes = {}
        for arm in ["day", "night", "open"]:
            x = group.loc[group.trt_Open.eq(arm)].copy()
            nflowers = pd.to_numeric(x.nFlrs, errors="raise").to_numpy(float)
            nfruit = pd.to_numeric(x.nFrts, errors="raise").to_numpy(float)
            require(len(x) >= 7 and np.all(np.isfinite(nflowers)) and
                    np.all(nflowers > 0) and np.all(np.isfinite(nfruit)) and
                    np.all((nfruit >= 0) & (nfruit <= nflowers)),
                    "Original independent fruit-set treatment invalid")
            outcomes[arm] = {
                "n_original_plants": int(len(x)),
                "fruit_set": float(nfruit.sum()/nflowers.sum())
            }
        frac = outcomes["night"]["fruit_set"] / (
            outcomes["night"]["fruit_set"] + outcomes["day"]["fruit_set"])
        require(np.isfinite(frac) and 0 <= frac <= 1,
                "Invalid original fruit-set night fraction")
        vr = visitors.loc[visitors.year.eq(2019) &
                          visitors["pop"].eq(row.pop)]
        require(len(vr) == 2 and vr.dataset.nunique() == 2,
                "Original narrow+wide visitor datasets not both present")
        vcounts = {k: int(pd.to_numeric(vr[k], errors="raise").sum())
                   for k in VISITOR_COUNTS}
        total = sum(vcounts.values())
        require(total == int(row.visTotal) and total > 0,
                "Actual source visitor counts do not match original population metadata")
        require((vcounts["ct_hawkmoth"]+vcounts["ct_otherN"] >
                vcounts["ct_otherD"]+vcounts["ct_small_bee"]+
                vcounts["ct_bumblebee"]+vcounts["ct_hummingbird"]) ==
                (row.majorPollPeriod == "nocturnal"),
                "Published most-common visitor period disagrees with group counts")
        result.append({
            "pop": str(row.pop), "year": 2019,
            "species": SPECIES_CODES[row.sp],
            "lat": float(row.latN), "lon": float(row.lonW),
            "night_fraction": float(frac),
            "night_vs_day_stronger": bool(
                outcomes["night"]["fruit_set"] > outcomes["day"]["fruit_set"]),
            "dominant_period": row.majorPollPeriod,
            "night_visitor_count": vcounts["ct_hawkmoth"]+vcounts["ct_otherN"],
            "day_visitor_count": total-vcounts["ct_hawkmoth"]-vcounts["ct_otherN"],
            "visitor_records": total,
            "original_individuals": sum(
                v["n_original_plants"] for v in outcomes.values())
        })
    result = pd.DataFrame(result)
    require((len(result), result["pop"].nunique(),
             result.species.nunique(),
             int(result.original_individuals.sum())) == (7, 7, 3, 204),
            "Original independent 2019 population-year / individual support altered")
    require((result.dominant_period.eq("nocturnal").sum(),
             result.dominant_period.eq("diurnal").sum()) == (3,4),
            "Direct field visitor night/day dominance support altered")
    return result


def source_global_match(preflight, author, direct):
    require(sha256(preflight) == ORIGINAL_SOURCE_SHA256,
            "Original 301 source SHA256 not reproduced")
    require(gitblob(author) == ORIGINAL_AUTHOR_BLOB,
            "Original published 1652 treatment data Git blob changed")
    original = pd.read_csv(preflight, dtype={"study_key": str, "species": str})
    require((len(original), original.study_key.nunique(),
             original.species.nunique()) == (301, 118, 123),
            "Original independent study site/year full cohort changed")
    a = pd.read_csv(author, dtype=str, keep_default_na=False)
    require(len(a) == 1652 and a.study_ID.nunique() == 140,
            "Original published author study file changed")
    w = original.loc[original.study_key.eq(SOURCE_STUDY)].copy()
    require((len(w), w.species.nunique(), w.year_exp.nunique()) == (6, 3, 1),
            "Original Wenzell 6-site/3-species/1-year subcohort changed")
    require(w.year_exp.eq(2019).all() and w.metric.eq("fruit set").all(),
            "Original matched Wenzell true 2019 fruit-set source changed")
    require(w.species.isin(direct.species).all(),
            "Original Wenzell botanical binomial not matched")
    out=[]
    for row in w.itertuples(index=False):
        candidates = direct.loc[direct.species.eq(row.species) &
                                direct.year.eq(row.year_exp)].copy()
        candidates["distance_km"] = candidates.apply(
            lambda z: haversine_km(row.lat, row.lon, z.lat, z.lon), axis=1)
        near = candidates.sort_values("distance_km")
        require(len(near) >= 1, "No independent botanical visitor candidate")
        p = near.iloc[0]
        require(p.distance_km <= 2.,
                "Source spatial join exceeds a provisional 2 km buffer")
        if len(near) > 1:
            require(float(near.iloc[1].distance_km) >
                    float(p.distance_km) + 1.,
                    "Ambiguous original site-year population choice")
        d0=float(row.night_reproduction)
        d1=float(row.day_reproduction)
        require(d0 >= 0 and d1 >= 0 and d0+d1 > 0,
                "Original meta-analysis night/day reproductive values invalid")
        global_frac = d0/(d0+d1)
        out.append({
            "population": p["pop"],
            "distance_km": float(p.distance_km),
            "original_source_night_share": global_frac,
            "direct_data_night_share": float(p.night_fraction),
            "period": p.dominant_period,
            "global_day_night_sign_matches_direct": bool(
                (d0 > d1) == bool(p.night_vs_day_stronger))
        })
    matched = pd.DataFrame(out)
    require((len(matched), matched.population.nunique()) == (6,6),
            "Original six source rows not matched one-to-one to field population-year")
    within1 = int(matched.distance_km.le(1.).sum())
    within2 = int(matched.distance_km.le(2.).sum())
    require((within1,within2) == (5,6),
            "Original spatial join no longer 5 strict / 6 provisional")
    require(int(matched.global_day_night_sign_matches_direct.sum()) == 6,
            "Original macro source treatment-sign does not match source direct fruit-set")
    return matched


def summarize_period(frame):
    grp = {}
    for period in ("nocturnal", "diurnal"):
        z = frame.loc[frame.dominant_period.eq(period)]
        grp[period] = {
            "population_years": int(len(z)),
            "mean_population_night_reproductive_share": float(z.night_fraction.mean()),
            "original_individual_plants_all_three_treatments": int(z.original_individuals.sum()),
            "original_visitor_events_with_distinct_narrow_and_wide_protocols":
                int(z.visitor_records.sum()),
            "nighttime_visitor_events": int(z.night_visitor_count.sum()),
            "daytime_visitor_events": int(z.day_visitor_count.sum()),
        }
    gap = (grp["nocturnal"]["mean_population_night_reproductive_share"] -
           grp["diurnal"]["mean_population_night_reproductive_share"])
    return {
        "by_actual_original_dominant_visitor_period": grp,
        "nocturnal_minus_diurnal_population_equal_night_share_difference": float(gap),
        "source_field_population_years": int(len(frame)),
        "population_year_night_day_fruit_direction_agrees_with_visitor_majority": int(
            ((frame.dominant_period.eq("nocturnal")) ==
             frame.night_vs_day_stronger).sum()),
        "unmatched_direction_population_years": int(
            ((frame.dominant_period.eq("nocturnal")) !=
             frame.night_vs_day_stronger).sum()),
    }


def run(preflight, author, wraw, wpop, wvis, output):
    output = Path(output)
    require(not output.exists(), "Do not overwrite a verified original source result")
    require(gitblob(wraw) == FILES["direct_raw"], "Direct Wenzell individual original changed")
    require(gitblob(wpop) == FILES["direct_population"], "Wenzell population original changed")
    require(gitblob(wvis) == FILES["direct_visitors"], "Wenzell visitor count original changed")
    raw = pd.read_csv(wraw, dtype={"year": int})
    pop = pd.read_csv(wpop, dtype={"year": int})
    vis = pd.read_csv(wvis, dtype={"year": int})
    direct = original_direct_population(raw, pop, vis)
    matched = source_global_match(preflight, author, direct)
    in_meta = direct.loc[direct["pop"].isin(set(matched.population))]
    strict_pops = matched.loc[matched.distance_km.le(1.),"population"]
    strict = direct.loc[direct["pop"].isin(set(strict_pops))]
    require((len(in_meta), len(strict)) == (6,5), "Meta source match eligibility changed")
    author_source = {
        "original_301_paired_source_SHA256": ORIGINAL_SOURCE_SHA256,
        "original_1652_author_treatment_git_blob": ORIGINAL_AUTHOR_BLOB,
        "independent_wenzell_data_git_blobs": FILES,
        "original_Castilleja_global_plant_fruitset_pairs": 6,
        "original_meta_independent_source_papers": 1,
        "source_original_publication_key": SOURCE_STUDY,
        "direct_published_article_DOI": DIRECT_SOURCE_DOI,
    }
    direct_all = summarize_period(direct)
    direct_inmeta = summarize_period(in_meta)
    direct_strict = summarize_period(strict)
    m = matched.merge(direct[["pop","dominant_period","night_fraction"]],
                      left_on="population",right_on="pop",validate="one_to_one")
    discrepancies = abs(m.original_source_night_share-m.direct_data_night_share)
    require(discrepancies.max() <= .08,
            "Original source field reproductive means mismatch direct published source")
    result = {
        "schema": "alan_real_original_Wenzell_populationyear_actual_visitors_reproduction_v1",
        "status": "AUTHENTIC_ONE_PAPER_2019_POPULATION_YEAR_VISITOR_REPRODUCTIVE_LINK",
        "source": author_source,
        "ecological_units": {
            "original_direct_2019_population_years": 7,
            "original_direct_2019_individual_treatment_plants": 204,
            "original_direct_2019_botanical_species": 3,
            "original_301_original_global_meta_rows_siteyear_matched": 6,
            "distance_strict_at_most_1km_source_matches": int(
                matched.distance_km.le(1.).sum()),
            "distance_extra_1_to_2km_provisional_source_matches": int(
                matched.distance_km.gt(1.).sum()),
            "distance_max_km": float(matched.distance_km.max()),
            "meta_source_and_direct_fruitset_sign_agreement": int(
                matched.global_day_night_sign_matches_direct.sum()),
            "meta_vs_independent_direct_reproductive_night_share_max_absolute_gap":
                float(discrepancies.max()),
            "one_direct_2019_population_NOT_in_frozen_original_301_source": 1,
        },
        "direct_seven_actual_field_population_years": direct_all,
        "direct_six_original_global_matched_population_years": direct_inmeta,
        "direct_five_within_1km_strict_population_years": direct_strict,
        "source_matched_six_original_meta_night_share_by_visitor_dominance": {
            "nocturnal_population_mean": float(
                m.loc[m.period.eq("nocturnal"),
                      "original_source_night_share"].mean()),
            "diurnal_population_mean": float(
                m.loc[m.period.eq("diurnal"),
                      "original_source_night_share"].mean()),
        },
        "epistemic_limits": {
            "one_independent_original_published_study_not_seven_studies": True,
            "observed_field_visitors_same_population_year_but_not_identical_trial_plant_or_day": True,
            "visitor_event_presence_not_equivalent_to_stigma_pollen_transfer": True,
            "visitor_count_protocols_distinct_narrow_wide_not_causal_visitation_rate": True,
            "one_original_source_georeference_1_to_2km_provisional": True,
            "one_direct_population_not_added_to_frozen_global_301_denominator": True,
            "original_source_night_share_vs_new_direct_ratio_different_aggregation": True,
            "plant_function_study_already_published_not_independent_novel_finding": True,
            "not_independent_test_of_fully_global_pollinator_shift": True,
            "no_Artificial_Night_Light_effect_fitted_or_claimed": True,
            "original_source_site_taxon_event_rows_NOT_exported": True,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print("ALAN_AUTHENTIC_DIRECT_2019_VISITOR_REPRODUCTION_BRIDGE",
          json.dumps(result, ensure_ascii=False), flush=True)
    return result


if __name__ == "__main__":
    a=argparse.ArgumentParser(description=__doc__)
    for k in ("preflight","author","wraw","wpop","wvis","output"):
        a.add_argument("--"+k,type=Path,required=True)
    p=a.parse_args()
    run(p.preflight,p.author,p.wraw,p.wpop,p.wvis,p.output)
