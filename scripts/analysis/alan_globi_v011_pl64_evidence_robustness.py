#!/usr/bin/env python3
"""GloBI v0.11 evidence update of the frozen 64-plant pollen-limitation × ALAN contrast.

All same 10,698 original plant-photograph cell-years and all covariates are
restored by the frozen original source code. This tests SOURCE DOCUMENTATION
stability, not local moth pollination or light-caused plant decline.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from alan_gecko_primary_fiveyear import primary

ROOT = Path(__file__).resolve().parents[2]
BASELINE_BETA = -0.019380045756016737
REVIEWED_POSITIVE = frozenset((
    "Capparis spinosa", "Lantana camara", "Lonicera japonica",
    "Malus domestica", "Mirabilis jalapa", "Vincetoxicum hirundinaria"
))


def weighted_species_response(matrix: pd.DataFrame,
                              denominator: np.ndarray,
                              PL: pd.Series) -> np.ndarray:
    """Each photographed botanical species is present once per site/year."""
    assert 3 <= len(PL) <= 2000 and PL.index.is_unique
    assert PL.notna().all() and np.isfinite(PL.to_numpy(dtype=float)).all()
    assert len(matrix) == len(denominator)
    centre = PL - PL.mean()
    assert float(centre.pow(2).sum()) > 0
    weights = centre / centre.pow(2).sum()
    return 100 * (matrix.reindex(columns=weights.index, fill_value=0)
                   .to_numpy(dtype=float) @ weights.to_numpy(dtype=float)
                  ) / denominator


def make_fit(s: pd.DataFrame):
    """Same region×year land-class balanced WLS + 5° clustered SE as README."""
    assert len(s) == 10698 and int(s.island.sum()) == 933
    assert s.region5.nunique() == 78 and s.region_year.nunique() == 195
    cols = ["island", "loglight", "light_island", "logpop", "logphotos",
            "logspecies", "bio5", "lograin", "elevation_km"]
    X = s[cols].to_numpy(float)
    assert np.isfinite(X).all()
    n = s.groupby(["region_year", "land_class"], observed=True).size()
    ix = pd.MultiIndex.from_frame(s[["region_year", "land_class"]])
    w = (n.groupby(level=0).min().reindex(s.region_year).to_numpy(float)
         / n.reindex(ix).to_numpy(float))
    assert np.isfinite(w).all() and (w > 0).all()
    group, levels = pd.factorize(s.region_year, sort=True)
    ws = np.bincount(group, weights=w, minlength=len(levels))
    sw = np.sqrt(w)
    for j in range(X.shape[1]):
        mean = np.bincount(group, weights=w * X[:, j], minlength=len(levels)) / ws
        X[:, j] = (X[:, j] - mean[group]) * sw
    bread = np.linalg.inv(X.T @ X)
    proj = bread @ X.T
    cl, ck = pd.factorize(s.region5, sort=True)
    G = len(ck)
    tc = float(student_t.ppf(0.975, G - 1))
    degree = (G / (G - 1)) * ((len(s) - 1) / (len(s) - len(cols)))

    def fit(y: np.ndarray) -> dict:
        y = np.asarray(y, float)
        assert len(y) == len(s) and np.isfinite(y).all()
        mu = np.bincount(group, weights=w * y, minlength=len(levels)) / ws
        response = (y - mu[group]) * sw
        beta = proj @ response
        score = X * (response - X @ beta)[:, None]
        sums = np.column_stack([
            np.bincount(cl, weights=score[:, j], minlength=G)
            for j in range(len(cols))
        ])
        cov = bread @ sums.T @ sums @ bread * degree
        k = cols.index("light_island")
        se = math.sqrt(max(0.0, cov[k, k]))
        assert se > 0
        return {
            "beta": float(beta[k]),
            "CI95": [float(beta[k] - tc * se), float(beta[k] + tc * se)],
            "p_nominal": float(2 * student_t.sf(abs(beta[k] / se), G - 1)),
            "se": float(se),
            "island_beta": float(beta[k] + beta[cols.index("loglight")]),
            "mainland_beta": float(beta[cols.index("loglight")]),
        }
    return fit


def run(inputs: Path, signed: Path) -> dict:
    sites, source_key_sha, plants, _ = primary(inputs, signed)
    assert len(sites) == 10698
    traits = pd.read_csv(ROOT / "results/moth_pollen_evidence/shared_544_evidence.csv").set_index("plant")
    assert len(traits) == 544 and traits.index.is_unique
    traits = traits.loc[traits.glopl_median_log_response_ratio.notna()].copy()
    assert len(traits) == 64

    flags = pd.read_csv(ROOT / "results/globi_v011_fiveyear/source_moth_species.csv")
    assert not flags.accepted_name.duplicated().any()
    flags = flags.set_index("accepted_name")
    assert set(traits.index).issubset(flags.index)
    for tag in ("legacy_wcvp", "v011_wcvp"):
        traits[tag] = flags.loc[traits.index, tag].astype(str).str.lower().eq("true").to_numpy()
    assert int(traits.v011_wcvp.sum()) == 30
    assert int(traits.legacy_wcvp.sum()) == 54

    positive = traits.index[traits.v011_wcvp].tolist()
    unrecorded = traits.index[~traits.v011_wcvp].tolist()
    initial_transport = set(traits.index[traits.moth_transport_documented.astype(str).str.lower().eq("true")])
    assert len(initial_transport) == 11
    reviewed_transport = sorted(initial_transport | REVIEWED_POSITIVE)
    assert len(reviewed_transport) == 17
    occ = plants.loc[plants.wcvp_accepted_name.isin(traits.index),
                     ["key", "wcvp_accepted_name"]].drop_duplicates()
    matrix = pd.crosstab(occ.key, occ.wcvp_accepted_name)
    matrix = matrix.reindex(index=sites.key, columns=traits.index, fill_value=0)
    assert ((matrix.to_numpy() == 0) | (matrix.to_numpy() == 1)).all()
    fit = make_fit(sites)
    denom = sites.n_species.to_numpy(dtype=float)

    def y(names):
        return weighted_species_response(matrix, denom,
                                         traits.loc[list(names), "glopl_median_log_response_ratio"])
    def estimate(names):
        names = list(names)
        return {
            "species": len(names),
            "positive_siteyears": int(matrix[names].any(axis=1).sum()),
            **fit(y(names))
        }
    models = {
        "original_PL64": estimate(traits.index),
        "updated_v011_positive_PL30": estimate(positive),
        "no_v011_positive_record_PL34": estimate(unrecorded),
        "old_taxonomically_resolved_PL54": estimate(traits.index[traits.legacy_wcvp]),
        "reviewed_moth_pollen_transport_PL17": estimate(reviewed_transport),
        "updated_and_reviewed_transport_PL8": estimate(
            [n for n in reviewed_transport if bool(traits.loc[n, "v011_wcvp"])]
        ),
        "original_PL_omit_abs_PL_gt2": estimate(
            traits.index[traits.glopl_median_log_response_ratio.abs().le(2)]),
        "v011_PL_omit_abs_PL_gt2": estimate(
            traits.index[traits.v011_wcvp &
                         traits.glopl_median_log_response_ratio.abs().le(2)]),
        "direct_v011_documented_minus_undocumented_PL_slope": {
            "species_positive": len(positive), "species_no_record": len(unrecorded),
            **fit(y(positive) - y(unrecorded)),
        }
    }
    assert abs(models["original_PL64"]["beta"] - BASELINE_BETA) < 1e-9
    assert abs(models["original_PL64"]["p_nominal"] - .043993409104237075) < 1e-7

    support = {}
    for label, names in (
        ("PL64", list(traits.index)), ("v011_PL30", positive),
        ("unrecorded_PL34", unrecorded), ("reviewed_transport_PL17", reviewed_transport)
    ):
        idx = matrix[names].any(axis=1).to_numpy()
        p = sites.iloc[np.flatnonzero(idx)]
        geo = p.groupby("region5").land_class.nunique()
        support[label] = {
            "island_positive_photo_siteyears": int(p.island.sum()),
            "mainland_positive_photo_siteyears": int((1-p.island).sum()),
            "geographic_blocks_with_both_classes": int(geo.eq(2).sum())
        }

    sensitivity = {}
    for label, names in (
        ("PL64", list(traits.index)), ("v011_PL30", positive),
        ("reviewed_transport_PL17", reviewed_transport)
    ):
        effects = [(n, estimate([t for t in names if t != n])) for n in names]
        effects.sort(key=lambda z: z[1]["beta"])
        sensitivity[label] = {
            "source_species_omissions": len(effects),
            "all_slope_negative": all(x["beta"] < 0 for _, x in effects),
            "beta_range": [effects[0][1]["beta"], effects[-1][1]["beta"]],
            "nominal_CIs_excluding_zero": sum(
                x["CI95"][1] < 0 or x["CI95"][0] > 0 for _, x in effects),
            "least_negative_if_removed": effects[-1][0],
            "most_negative_if_removed": effects[0][0],
        }

    return {
        "schema": "ALAN_GloBI_v011_readme_PL64_source_update_v1",
        "status": "real_signed_source_recomputed_exploratory",
        "years": [2019, 2021, 2023, 2024, 2025],
        "target_key_set_sha256": source_key_sha,
        "design": {
            "site_years": len(sites), "island_site_years": int(sites.island.sum()),
            "mainland_site_years": int((1 - sites.island).sum()),
            "geographic_blocks": 78, "year_region_groups": 195,
        },
        "source_label_counts": {
            "PL64": 64, "v011_moth_positive": 30, "v011_no_record": 34,
            "old_WCVP_moth_positive": 54, "reviewed_transport": 17,
        },
        "models": models,
        "positive_documentation_support": support,
        "leave_one_species_out_overlapping_sensitivities": sensitivity,
        "limitations": {
            "v011_no_record_is_NOT_biological_no_moth": True,
            "moth_pollen_carriage_is_not_reproductive_pollination": True,
            "species_level_GloPL_median_not_measured_at_photo_site": True,
            "original_5year_photos_are_not_flora_censuses": True,
            "nominal_uncorrected_tests_after_exploration": True,
            "no_causal_ALAN_effect_or_local_extinction": True,
        },
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs", type=Path, help="Previously SHA-signed five-year source folders")
    ap.add_argument("signed", type=Path, help="Previously SHA-signed canonical photos/WCVP files")
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    if args.output.exists():
        raise ValueError("Refuse overwriting a previous research result")
    result = run(args.inputs, args.signed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"original": result["models"]["original_PL64"],
                      "updated": result["models"]["updated_v011_positive_PL30"]}))
