#!/usr/bin/env python3
"""Predeclared ALAN-versus-built-change common-support audit (NOT causal proof).

Use the frozen 1999-2023 paired wet-plant site panel in historical commit 6a9a484.
No outcome enters matching. Records are selected wet-decline candidates, so
the contrast is exploratory and cannot validate the postselected main claim.
All reported effect estimates are conditional on common support.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

REPO = Path(__file__).resolve().parents[2]
SOURCE_COMMIT = "6a9a484"
SOURCE_PATH = "results/all_species_forest/wet_water/coordinate_trends.csv.gz"
SEED = 20261119
BASE = ["built", "built0", "light0", "start", "span", "temp_trend",
        "rain_trend", "effort", "days", "spring", "summer", "freq0", "area0"]
WATER_EXTRA = ["freq_trend", "area_trend"]
REQUIRED = set(BASE + WATER_EXTRA + ["unit_id", "block", "light",
                                     "rate_difference", "candidate_rate",
                                     "reference_rate"])

def read_source():
    raw = subprocess.check_output(
        ["git", "show", f"{SOURCE_COMMIT}:{SOURCE_PATH}"], cwd=REPO)
    df = pd.read_csv(io.BytesIO(raw), compression="gzip")
    missing = REQUIRED - set(df)
    if missing:
        raise RuntimeError("FROZEN_SOURCE_SCHEMA_STOP: " + repr(sorted(missing)))
    if df.unit_id.duplicated().any():
        raise RuntimeError("FROZEN_SOURCE_DUPLICATE_UNIT_STOP")
    return df, hashlib.sha256(raw).hexdigest()

def robust_scale(s):
    v = s.dropna().astype(float)
    if len(v) < 3:
        return 1.
    q = float(v.quantile(.75) - v.quantile(.25))
    sd = float(v.std())
    return max(q / 1.349, sd * .10, 1e-5)

def maxgap(s, iqrs=1.0):
    a = s.dropna().astype(float)
    q = float(a.quantile(.75) - a.quantile(.25))
    sd = float(a.std())
    return max(iqrs * q, sd * .10, 1e-5)

def pair_sites(df, water_matched=False):
    """Exact 100-km block, no replacement; all calipers frozen before outcome."""
    cols = BASE + (WATER_EXTRA if water_matched else [])
    f = df.dropna(subset=cols + ["unit_id", "block", "light"]).copy()
    f = f[np.isfinite(f[cols + ["light"]]).all(axis=1)]
    f = f.sort_values("unit_id").reset_index(drop=True)
    scales = {k: robust_scale(f[k]) for k in cols}
    limits = {k: maxgap(f[k], .5 if k == "built" else 1.) for k in cols}
    limits["start"] = 5.
    limits["span"] = 10.
    # Distinct source-era high/low exposure: increasing vs non-increasing light.
    f["treated"] = f.light.gt(0.)
    pairs = []
    for block, z in f.groupby("block", sort=True):
        high = z[z.treated].sort_values("unit_id")
        low = z[~z.treated].sort_values("unit_id")
        if not len(high) or not len(low):
            continue
        d2 = np.zeros((len(high), len(low)), dtype=float)
        acceptable = np.ones_like(d2, dtype=bool)
        for k in cols:
            diff = high[k].to_numpy(float)[:, None] - low[k].to_numpy(float)[None, :]
            d2 += (diff / scales[k]) ** 2
            if k in ("built", "built0", "light0", "start", "span"):
                acceptable &= np.abs(diff) <= limits[k]
            if water_matched and k in WATER_EXTRA:
                acceptable &= np.abs(diff) <= limits[k]
        # Avoid trivially indistinguishable light slopes.
        acceptable &= (high.light.to_numpy()[:, None] -
                       low.light.to_numpy()[None, :]) >= 1.
        # Dummy columns allow unmatched treated sites; maximize eligible pairs first.
        penalty = 1.e6
        cost = np.where(acceptable, d2, penalty * 2.)
        cost = np.column_stack([cost, np.full((len(high), len(high)), penalty)])
        rr, cc = linear_sum_assignment(cost)
        for ri, ci in zip(rr, cc):
            if ci >= len(low) or not acceptable[ri, ci]:
                continue
            a, b = high.iloc[ri], low.iloc[ci]
            item = {"block": block, "light_difference": float(a.light - b.light)}
            for k in cols:
                item["d_" + k] = float(a[k] - b[k])
            # Audit water changes even if they are not matching covariates.
            for k in WATER_EXTRA:
                item["d_" + k] = float(a[k] - b[k])
            # Outcome is accessed only after pair assignments are frozen.
            item["outcome_difference"] = float(
                a.rate_difference - b.rate_difference)
            item["candidate_difference"] = float(
                a.candidate_rate - b.candidate_rate)
            item["reference_difference"] = float(
                a.reference_rate - b.reference_rate)
            pairs.append(item)
    return pd.DataFrame(pairs), f, scales

def summarize_pairs(pairs, df, scales):
    if pairs.empty:
        return {"status": "NO_COMMON_SUPPORT", "pairs": 0, "blocks": 0}
    p = pairs.dropna(subset=["outcome_difference", "light_difference"]).copy()
    nblocks = int(p.block.nunique())
    x, y = p.light_difference.to_numpy(), p.outcome_difference.to_numpy()
    slope = float(np.dot(x, y) / np.dot(x, x))
    mean_delta = float(np.mean(y))
    # Balance of paired pre-exposure covariates; diagnostics, not post-hoc filter.
    balance = {k[2:]: float(p[k].mean() / scales[k[2:]])
               for k in p.columns if k.startswith("d_")}
    # Cluster bootstrap resamples spatial blocks, not individual plants/sites.
    rng = np.random.default_rng(SEED)
    groups = {k: g for k, g in p.groupby("block")}
    keys = list(groups)
    draws_mean, draws_slope = [], []
    for _ in range(2500):
        picked = rng.choice(len(keys), size=len(keys), replace=True)
        sampled = pd.concat([groups[keys[i]] for i in picked], ignore_index=True)
        sx = sampled.light_difference.to_numpy()
        sy = sampled.outcome_difference.to_numpy()
        draws_mean.append(float(sy.mean()))
        draws_slope.append(float(np.dot(sx, sy) / np.dot(sx, sx)))
    low, high = np.quantile(draws_mean, [.025, .975])
    slo, shi = np.quantile(draws_slope, [.025, .975])
    # Structural stop is independent of the sign or p-value of the outcome.
    # Baseline open water and weather must ALSO balance. A built-up-only
    # diagnostic is not source-supported if wet-habitat baseline is different.
    safe = (len(p) >= 30 and nblocks >= 10 and
            all(abs(v) <= .25 for k, v in balance.items()
                if k in ("built", "built0", "light0", "area0", "freq0",
                         "temp_trend", "rain_trend")))
    return {
        "status": "EXPLORATORY_COMMON_SUPPORT" if safe
                  else "SOURCE_SUPPORT_STOP",
        "matched_pairs": int(len(p)), "independent_spatial_blocks": nblocks,
        "source_eligible_sites": int(len(df)),
        "light_slope_contrast_median_DN_per_decade":
           float(p.light_difference.median()),
        "built_change_SMD": balance.get("built"),
        "baseline_built_SMD": balance.get("built0"),
        "baseline_light_SMD": balance.get("light0"),
        "all_matched_covariate_SMD": balance,
        "candidate_minus_reference_trend_pair_difference_pp_per_decade":
            mean_delta,
        "mean_pair_difference_cluster_bootstrap_95pct": [float(low), float(high)],
        "slope_per_10_DN_brightening_pp_per_decade": 10 * slope,
        "slope_per_10_DN_cluster_bootstrap_95pct": [10 * float(slo), 10 * float(shi)],
        "interpretation": "postselected wet-decline groups and historical site trends; "
                          "not ALAN causal effect or new independent validation",
    }

def same_cohort_water_adjustment(pairs):
    """Post-hoc sensitivity: identical matched pairs, water-trend adjustment.

    Baseline open-water and climate are held as covariates in BOTH models.
    The only contrast is adding water-area and detection-frequency trends.
    This is NOT prespecified, and must not be described as causal mediation.
    """
    terms = ["light_difference", "d_built", "d_built0", "d_light0",
             "d_temp_trend", "d_rain_trend", "d_freq0", "d_area0"]
    added = ["d_freq_trend", "d_area_trend"]
    f = pairs.dropna(subset=terms + added + ["outcome_difference"]).copy()
    if len(f) < 30 or f.block.nunique() < 10:
        return {"status": "SOURCE_SUPPORT_STOP", "pairs": len(f),
                "blocks": int(f.block.nunique())}
    def slope(df, cols):
        x = df[cols].to_numpy(float)
        y = df.outcome_difference.to_numpy(float)
        if np.linalg.matrix_rank(x) < x.shape[1]:
            return np.nan
        return 10. * float(np.linalg.lstsq(x, y, rcond=None)[0][0])
    before, after = slope(f, terms), slope(f, terms + added)
    rng = np.random.default_rng(SEED + 1)
    groups = {k: z for k, z in f.groupby("block")}
    keys = list(groups)
    draws = []
    for _ in range(1800):
        selected = rng.choice(len(keys), len(keys), replace=True)
        z = pd.concat([groups[keys[i]] for i in selected], ignore_index=True)
        draws.append((slope(z, terms), slope(z, terms + added)))
    arr = np.asarray(draws)
    valid = np.isfinite(arr).all(axis=1)
    ci = (np.quantile(arr[valid], [.025, .975], axis=0)
          if valid.sum() > 500 else np.full((2, 2), np.nan))
    return {
        "status": "POSTHOC_SENSITIVITY_ONLY",
        "same_matched_pairs": int(len(f)),
        "independent_spatial_blocks": int(f.block.nunique()),
        "model_without_water_change_slope_per_10_DN": before,
        "model_without_water_change_cluster_95pct": ci[:, 0].tolist(),
        "model_with_water_change_slope_per_10_DN": after,
        "model_with_water_change_cluster_95pct": ci[:, 1].tolist(),
        "sensitivity_difference_not_causal_water_mediation": after - before,
        "both_models_include_baseline_water_and_climate": True,
        "bootstrap_resamples_with_complete_rank": int(valid.sum()),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    data, sha = read_source()
    runs = {}
    for label, water in [("matched_development", False),
                         ("matched_development_and_water_change", True)]:
        pairs, eligible, scales = pair_sites(data, water)
        runs[label] = summarize_pairs(pairs, eligible, scales)
        if not water:
            runs[label]["same_cohort_regression_sensitivity"] = (
                same_cohort_water_adjustment(pairs))
    report = {
        "source": f"{SOURCE_COMMIT}:{SOURCE_PATH}",
        "source_sha256": sha,
        "frozen_matching_algorithm": "within-block Hungarian no-replacement; "
                                     "built change within 0.5 IQR; light difference >=1 DN/decade; "
                                     "baseline built and light within IQR; start <=5y; span <=10y",
        "sample_description": "same existing 58 postselected moist-decliners and 305 "
                              "wet comparison species; no new plant records",
        "run": runs
    }
    s = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True)
    print("ALAN_WET_LIGHT_DEVELOPMENT_OVERLAP_AUDIT\n" + s, flush=True)
    if args.out:
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(s + "\n")

if __name__ == "__main__":
    main()
