#!/usr/bin/env python3
"""GBIF 2013-2025 genuine animal photo counts × prior PCNL: source-effort safeguards.

New analyses on EXACT preceding ALAN source cohort. Species non-record is NOT
an absence. PHOTO ratios are neither abundance nor local predation/pollination.
No raw GBIF cell/coordinate/species×site records are exported.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, hstack
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsmr
from scipy.stats import t as student_t

from alan_13year_moth_gecko_repeat_site_FE_20261112 import read_pair, repeated_sites
from alan_13year_moth_gecko_same_island_source_ALAN_20261112 import cohort


def require(test, why):
    if not bool(test):
        raise ValueError(why)


def fit_photo_share(data, response, weight="island_equal", repeat=False):
    """FWL weighted regression with geographic island-cluster CR1.

    Count weights are a sensitivity design rather than a new biological unit.
    """
    d = data.copy()
    n = len(d)
    islands = d.island.nunique()
    if n < 70 or islands < 12:
        return {"estimable": False, "reason": "low_geographic_source_support",
                "n": n, "islands": islands}

    if response == "gecko":
        y = 100 * d.gecko_share
    elif response == "gecko_jeffreys":
        y = 100 * (d.gecko_source_photo_count + .5) / (
            d.original_total_reptile_photo_source_events + 1)
    elif response == "moth":
        y = d.moth_share
    else:
        raise ValueError(response)

    cols = ["loglight", "log_insect", "log_reptile"]
    if not repeat:
        cols += ["latitude"]
    if weight == "island_equal":
        w = 1 / d.groupby("island").island.transform("size").to_numpy(float)
    elif weight == "reptile_capped_count":
        w = (1 / d.groupby("island").island.transform("size").to_numpy(float)
             * np.minimum(d.original_total_reptile_photo_source_events, 10))
    elif weight == "insect_capped_count":
        w = (1 / d.groupby("island").island.transform("size").to_numpy(float)
             * np.minimum(d.all_four_family_records, 40))
    elif weight == "cell_equal":
        w = np.ones(n)
    else:
        raise ValueError(weight)

    sw = np.sqrt(w)
    Zs = []
    dim, codes = [], []
    for term in (("site", "islandyear") if repeat else ("islandyear",)):
        c, vals = np.unique(d[term], return_inverse=True)
        codes.append(vals)
        dim.append(len(c))
        Zs.append(csc_matrix((sw, (np.arange(n), vals)), shape=(n, len(c))))
    Z = hstack(Zs, format="csc")
    V = np.column_stack((d[cols].to_numpy(float), y.to_numpy(float))) * sw[:, None]
    require(np.isfinite(V).all(), "Nonfinite source share or predictor")
    for j in range(V.shape[1]):
        coeff, status, *_ = lsmr(Z, V[:, j], atol=1e-11, btol=1e-11,
                                 maxiter=8000)
        require(status != 7, "Source FE projection failed to converge")
        V[:, j] -= Z @ coeff

    X, Y = V[:, :-1], V[:, -1]
    beta, _, rank, _ = np.linalg.lstsq(X, Y, rcond=1e-11)
    if rank != len(cols):
        return {"estimable": False, "reason": "fixed_effect_rank_deficient",
                "n": n, "islands": islands}

    if repeat:
        ns, ny = dim
        site, island_year = codes
        graph = coo_matrix((np.ones(2*n), (
            np.r_[site, ns+island_year],
            np.r_[ns+island_year, site]
        )), shape=(ns+ny, ns+ny)).tocsr()
        components = connected_components(graph, directed=False)[0]
        fe_rank = ns + ny - components
    else:
        components = 1
        fe_rank = dim[0]
    df = n - fe_rank - len(cols)
    require(df > 0, "Source model has no residual degrees of freedom")
    residual = Y - X @ beta
    cluster, keys = np.unique(d.island, return_inverse=True)
    G = len(cluster)
    scores = np.column_stack([
        np.bincount(keys, weights=X[:,j]*residual, minlength=G)
        for j in range(len(cols))
    ])
    bread = np.linalg.inv(X.T @ X)
    C = bread @ (scores.T @ scores) @ bread * (
        G/(G-1) * (n-1)/df)
    se = math.sqrt(max(0.0, float(C[0,0])))
    critical = float(student_t.ppf(.975, G-1))
    return {
        "estimable": True, "n": n, "islands": G,
        "sites": int(d.site.nunique()), "islandyears": int(d.islandyear.nunique()),
        "fixed_effects_rank": int(fe_rank),
        "residual_df": int(df), "fe_components": int(components),
        "response": response, "weight": weight,
        "effect": {"beta": float(beta[0]),
                   "ci95": [float(beta[0]-critical*se),
                            float(beta[0]+critical*se)],
                   "p_nominal": float(2*student_t.sf(abs(beta[0]/se),G-1))
                   if se > 0 else None}
    }


def gate(d, minimum_reptile, repeat):
    a = d[d.original_total_reptile_photo_source_events.ge(minimum_reptile)].copy()
    if repeat:
        return repeated_sites(a, .02)
    g = a.groupby("islandyear").loglight.agg(["size", "min", "max"])
    valid = g.index[g["size"].ge(2) & (g["max"]-g["min"]).ge(.02)]
    return a[a.islandyear.isin(valid)].copy()


def run(insects, reptiles):
    a = read_pair(insects, reptiles)
    spatial = cohort(a)
    repeated = repeated_sites(spatial, .02)
    require((len(spatial), spatial.island.nunique()) == (3424, 50),
            "Prior 50-island source geography changed")
    require((len(repeated), repeated.site.nunique(), repeated.island.nunique())
            == (1979, 454, 30),
            "Prior 454-source-site temporal cohort changed")
    models = {}
    for label, d, is_repeated in (
        ("spatial50", spatial, False),
        ("repeat30", repeated, True)
    ):
        for response in ("gecko", "gecko_jeffreys", "moth"):
            for weight in ("island_equal",
                           "insect_capped_count" if response == "moth"
                           else "reptile_capped_count"):
                models[f"{label}__{response}__{weight}"] = fit_photo_share(
                    d, response, weight, repeat=is_repeated)
        for min_n in (2,3,5,10):
            selected = gate(d, min_n, is_repeated)
            models[f"{label}__gecko__reptile_photo_min_{min_n}"] = fit_photo_share(
                selected, "gecko", "island_equal", repeat=is_repeated)

    checks = {}
    for key, d in (("spatial50", spatial), ("repeat30", repeated)):
        share = d.gecko_share
        reptile = d.original_total_reptile_photo_source_events
        checks[key] = {
            "site_years": len(d),
            "islands": int(d.island.nunique()),
            "reptile_photo_opportunities_median": float(reptile.median()),
            "source_exact_gecko_0_or_1_fraction": float(share.isin((0,1)).mean()),
            "source_exact_gecko_0_or_1_when_reptile_photos_ge5_fraction":
                float(share[reptile.ge(5)].isin((0,1)).mean()),
            "original_source_reptile_photographs_total": int(reptile.sum()),
            "original_source_insect_photographs_total": int(d.all_four_family_records.sum()),
            "source_siteyears_before_2018_fraction": float(d.year.lt(2018).mean()),
            "year_counts": {str(k):int(v) for k,v in d.groupby("year").size().items()}
        }
    return {
        "schema": "alan_13yr_original_GBIF_gecko_moth_photo_effort_reliability_v1",
        "source_years": list(range(2013,2026)),
        "frozen_insect_SHA256":
            "64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852",
        "frozen_reptile_SHA256":
            "03aa855c9e7691855bfee4baa099f0d3a31db8d53066ec36a95afda48c3ba339",
        "models":models,
        "source_sampling_diagnostics":checks,
        "warnings":{
            "photographed_guild_ratio_is_NOT_population_size":True,
            "record_zero_is_NOT_species_absence":True,
            "source_photo_opportunity_gate_is_selection_not_experiment":True,
            "capped_count_weights_not_true_independent_sample_design":True,
            "Jeffreys_0p5_smoothed_record_fraction_is_sensitivity_only":True,
            "no_gecko_moth_predation_or_local_moth_pollination_observed":True,
            "no_causal_nightlight_effect":True,
            "all_exploratory_p_values_unadjusted":True,
            "prior_README_PL64_positive_unmodified":True,
            "no_original_source_precise_sites_uploaded":True
        }
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--insects",required=True,type=Path)
    parser.add_argument("--reptiles",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("Refuse overwrite")
    result = run(args.insects,args.reptiles)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("ALAN_13YR_SOURCE_EFFORT_RELIABILITY",json.dumps({
        k:v["effect"] for k,v in result["models"].items()
        if k.startswith("spatial50__gecko") and v["estimable"]},sort_keys=True))
