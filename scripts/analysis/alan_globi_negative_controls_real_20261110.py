#!/usr/bin/env python3
"""Source-matched negative controls for ALAN × published moth floral interactions.

Input: same 85,318 0.05° island site×year WCVP species denominators in the
SHA-verified 93,209,351 original-photo source and real 1992–2024 PCNL.
Outcomes: documented moth, butterfly, their UNION, overlap, record-only
categories and contemporaneous paired outcome contrasts.

GloBI cache contains positive MOTH/BUTTERFLY floral relationships, not all
insect guilds; missing evidence is UNKNOWN rather than nonpollination.
No photo colour, no proof of nocturnal pollination or plant decline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import alan_full93m_moth_PCNLPrior_two_way_FE_20261109 as base

SCHEMA = "ALAN_real93m_GloBI_moth_butterfly_negative_controls_20261110_v1"
CANONICAL_PCNL_SHA = "440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c"
SAMPLE_SITES = 85318
SOURCE_NAMES = {
    "moth_positive": "moth",
    "butterfly_positive": "butterfly",
    "both_documented": "both",
    "moth_without_butterfly_record": "moth_record_no_butterfly",
    "butterfly_without_moth_record": "butterfly_record_no_moth",
    "neither_documented": "neither_documented",
    "moth_database_pollinates": "moth_pollinates",
    "nocturnal_moth_pollen_carriage": "night_pollen",
}


def require(condition, reason):
    if not condition:
        raise ValueError("ALAN_SOURCE_NEGATIVE_CONTROL: " + str(reason))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as src:
        for block in iter(lambda: src.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verified_source(pcnl, positive_guilds):
    require(digest(pcnl) == CANONICAL_PCNL_SHA,
            "Changed original 85,318-row PCNL source; must use frozen original ZIP")
    p = pd.read_csv(pcnl)
    g = pd.read_csv(positive_guilds)
    for q in (p, g):
        q["island"] = q.island.astype(str)
        require(not q[["island", "latbin", "lonbin", "year"]].duplicated().any(),
                "Duplicate plant site/year key")
    require(len(p) == len(g) == SAMPLE_SITES, "Original GSHHG island site/year support changed")
    for col in ("all_accepted_recorded_species", "moth_visitor_positive_recorded_species"):
        require(col in p.columns, "Original WCVP denominator unavailable: " + col)
    needed = {"island", "latbin", "lonbin", "year", "all_species", *SOURCE_NAMES.values()}
    require(needed.issubset(g.columns), "Not the exact original GloBI moth/butterfly source")
    cols = ["island", "latbin", "lonbin", "year", "all_species", *SOURCE_NAMES.values()]
    d = p.merge(g[cols], on=["island","latbin","lonbin","year"],
                how="left", validate="one_to_one", indicator=True)
    require(d._merge.eq("both").all(), "Incomplete same-site positive-source evidence match")
    require(d.all_species.eq(d.all_accepted_recorded_species).all(),
            "Accepted WCVP denominator not conserved")
    require(d.moth.eq(d.moth_visitor_positive_recorded_species).all(),
            "Positive GloBI moth plant numerator not conserved")
    require((d.both+d.moth_record_no_butterfly+
             d.butterfly_record_no_moth+d.neither_documented).eq(d.all_species).all(),
            "Mutually exclusive GloBI documentation cells not closed")
    require((d.both+d.moth_record_no_butterfly).eq(d.moth).all(),
            "The original moth coefficient would not be reconstructed")
    require((d.both+d.butterfly_record_no_moth).eq(d.butterfly).all(),
            "Positive butterfly documentation not conserved")
    require(d.moth_pollinates.le(d.moth).all() and
            d.night_pollen.le(d.moth).all(), "Pollen evidence cannot exceed positive moth records")
    d["any_moth_or_butterfly"] = d.moth+d.butterfly_record_no_moth
    d["butterfly_minus_moth"] = d.butterfly-d.moth
    d["any_minus_moth"] = d.any_moth_or_butterfly-d.moth
    return d.drop(columns="_merge")


def run(pcnl, guild, out):
    data = verified_source(pcnl, guild)
    supported = base.prepare(data)
    require(len(supported) == 77647 and supported.island.nunique() == 562,
            "PR24 source-support gate changed: must retain 562 original islands / 77,647 sites")
    names = dict(SOURCE_NAMES)
    names.update({"any_moth_or_butterfly": "any_moth_or_butterfly",
                  "butterfly_minus_moth": "butterfly_minus_moth",
                  "any_minus_moth": "any_minus_moth"})
    records = {}
    for label, column in names.items():
        d = supported.copy()
        d["frac_pp"] = 100*d[column]/d.all_accepted_recorded_species
        fitted = {
            "island_equal": base.fit(d, "island_equal"),
            "site_equal": base.fit(d, "site_equal"),
            "island_equal_plus_species_richness": base.fit(d, "island_equal", True),
        }
        require(all(v["status"] == "estimated_exploratory_observational"
                    for v in fitted.values()), "Source support or model rank insufficient")
        records[label] = {
            "n_accepted_species_site_records": int(d[column].sum()),
            "models": {k: v["terms"][base.LIGHT] for k, v in fitted.items()}
        }
    beta = lambda name: records[name]["models"]["island_equal"]["beta"]
    require(abs(beta("moth_positive") - 0.44933020729608203) < 1e-6,
            "Existing PR24 empirical moth association not exactly reproduced")
    require(np.isclose(beta("any_moth_or_butterfly")-beta("moth_positive"),
                           beta("any_minus_moth"), atol=1e-8),
            "A paired contrast is not the same as the butterfly-only source category")
    require(np.isclose(beta("moth_positive"), beta("both_documented")+
                           beta("moth_without_butterfly_record"), atol=1e-8),
            "Moth-positive decomposition does not exactly reconstruct parent model")
    require(np.isclose(beta("butterfly_positive")-beta("moth_positive"),
                           beta("butterfly_minus_moth"), atol=1e-8),
            "Butterfly–moth paired contrast not reproducible")
    result = {
        "schema": SCHEMA,
        "input_PCNL_source_SHA256": digest(pcnl),
        "input_GloBI_source_SHA256": digest(guild),
        "original_photo_ID_source": 93209351,
        "original_genuine_PCNL_site_years": len(data),
        "original_within_island_sites": len(supported),
        "original_within_island_islands": int(supported.island.nunique()),
        "outcome_definition": "recorded accepted plant species with published GloBI positive moth and/or butterfly flower interactions per source cell/year",
        "positive_GloBI_cache_scope": "moth-family and butterfly-family source relations only, NOT all insects",
        "paired_contrast_CI_uses_same_observations": True,
        "missing_record_is_unknown_NOT_nonpollination": True,
        "cannot_assert_moth_dependence_or_causal_light_effect": True,
        "not_a_native_only_or_uncultivated_only_sample": True,
        "results": records,
    }
    dest = Path(out)
    require(not dest.exists(), "Refuse to overwrite an earlier signed sensitivity")
    dest.mkdir(parents=True)
    (dest/"ALAN_real_full93m_GloBI_negative_controls_20261110.json").write_text(
        json.dumps(result, indent=2, sort_keys=True)+"\n")
    print("ALAN_GLOBI_REAL_SAME_SOURCE_NEGATIVE_CONTROLS",
          json.dumps({k:v["models"]["island_equal"] for k,v in records.items()},
                     sort_keys=True), flush=True)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--pcnl",type=Path,required=True)
    parser.add_argument("--guild",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    run(args.pcnl,args.guild,args.out)


if __name__=="__main__": main()
