#!/usr/bin/env python3
"""Source-verified additive attribution of old vs GloBI v0.11 plant labels.

Target: EXACT frozen 10,698 source photo-cell years, old -> updated moth evidence
- same photo and accepted species denominator
- same 4-month, same-year VIIRS/population/environment and geography
- same year x geographic-block FE, island/mainland balance and cluster design

Every change is a DOCUMENTATION LABEL CHANGE, not a biological loss/colonization.
Only aggregate species-name, genus, hemisphere and year results may be exported;
original photo locations or species x geographic-cell identities remain private.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import primary
from alan_globi_crossland_light_robust_20261111 import fit

ROOT = Path(__file__).resolve().parents[2]
ID = "wcvp_accepted_plant_name_id"
BETA = "beta_island_minus_mainland_slope_pp_per_log1p"
CLASSES = ("common", "old_only", "new_only", "neither")
DESIGN = ("island", "loglight", "light_island", "logpop", "logphotos",
          "logspecies", "bio5", "lograin", "elevation_km")


def light_island_influence(s: pd.DataFrame) -> np.ndarray:
    """FWL weighted regression influence of *raw percentage response* on slope.

    beta[k] = sum_i influence[i] * Y[i]; does not require centering Y because
    the residualized X is weighted-orthogonal to every year x region intercept.
    """
    if not s.target_id.is_unique or s.region_year.nunique() < 2:
        raise ValueError("Expected unique matched source sites in multiple FE groups")
    g, levels = pd.factorize(s.region_year, sort=True)
    sizes = s.groupby(["region_year", "land_class"]).size()
    ix = pd.MultiIndex.from_frame(s[["region_year", "land_class"]])
    sizes_min = sizes.groupby(level=0).min()
    w = sizes_min.reindex(s.region_year).to_numpy(float) / sizes.reindex(ix).to_numpy(float)
    if not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError("Invalid island/mainland class balance")
    sw = np.sqrt(w)
    X = s[list(DESIGN)].to_numpy(float).copy()
    ws = np.bincount(g, weights=w, minlength=len(levels))
    for j in range(X.shape[1]):
        mu = np.bincount(g, weights=w * X[:, j], minlength=len(levels)) / ws
        X[:, j] = (X[:, j] - mu[g]) * sw
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError("Rank deficient frozen panel")
    return (np.linalg.solve(X.T @ X, X.T)[DESIGN.index("light_island"), :] * sw)


def classify_species(sp: pd.DataFrame, old: pd.DataFrame, flag_csv: Path) -> pd.DataFrame:
    """Old published name equality versus new resolved WCVP accepted-ID label."""
    flag = pd.read_csv(flag_csv, dtype={"accepted_id": str})
    assert flag.accepted_id.is_unique and flag.accepted_name.is_unique
    assert set(flag.v011_wcvp.unique()) <= {False, True}
    new_ids = set(flag.loc[flag.v011_wcvp, "accepted_id"])
    old_names = set(old.loc[old.guild.eq("moth_family_visit"), "plant_name"])
    z = sp[["key", ID, "wcvp_accepted_name"]].copy()
    assert not z.duplicated(["key", ID]).any()
    z["old_positive"] = z.wcvp_accepted_name.isin(old_names)
    z["new_positive"] = z[ID].astype(str).isin(new_ids)
    z["class"] = np.select(
        [z.old_positive & z.new_positive,
         z.old_positive & ~z.new_positive,
         ~z.old_positive & z.new_positive],
        ["common", "old_only", "new_only"],
        default="neither")
    return z


def supported_fits(s: pd.DataFrame, observations: pd.DataFrame):
    counts = (observations.groupby(["key", "class"]).size()
              .unstack(fill_value=0).reindex(index=s.key, columns=CLASSES, fill_value=0))
    assert len(counts) == len(s)
    n = s.n_species.to_numpy(dtype=int)
    assert np.array_equal(counts.sum(axis=1).to_numpy(), n)
    old = (counts.common + counts.old_only).to_numpy()
    new = (counts.common + counts.new_only).to_numpy()
    assert np.array_equal(old, s.all_moth.to_numpy())
    responses = {
        "common": counts.common.to_numpy(),
        "old_only": counts.old_only.to_numpy(),
        "new_only": counts.new_only.to_numpy(),
        "published_old": old,
        "updated_v011": new,
        "updated_minus_published": new - old,
    }
    models = {}
    for name, y in responses.items():
        # The original regression can fit a signed pseudo-count for a direct
        # paired response contrast; this is not a biological signed frequency.
        models[name] = fit(s.assign(all_moth=y), "all", "moth", balance=True)
        if not models[name].get("estimable"):
            raise ValueError(f"Original frozen source fit failed: {name}")
    b = lambda name: float(models[name][BETA])
    np.testing.assert_allclose(b("published_old"),
                               b("common") + b("old_only"), atol=1e-9, rtol=0)
    np.testing.assert_allclose(b("updated_v011"),
                               b("common") + b("new_only"), atol=1e-9, rtol=0)
    np.testing.assert_allclose(b("updated_minus_published"),
                               b("new_only") - b("old_only"), atol=1e-9, rtol=0)
    return counts, responses, models


def aggregate_influence(s: pd.DataFrame, sp: pd.DataFrame, response: dict,
                        models: dict) -> dict:
    influence = light_island_influence(s)
    n = s.n_species.to_numpy(dtype=float)
    for name, count in response.items():
        pred = float(np.sum(100 * influence * count / n))
        np.testing.assert_allclose(pred, models[name][BETA], atol=1e-8, rtol=0)
    support = pd.DataFrame({
        "key": s.key.to_numpy(),
        "site_influence": 100 * influence / n,
        "year": s.year.to_numpy(),
        "land_class": s.land_class.to_numpy(),
        "region5": s.region5.to_numpy(),
    })
    support["lat_midpoint"] = support.region5.str.split("_").str[-1].astype(int) * 5 - 87.5
    lat_bins = [-90, -60, -30, 0, 30, 60, 90]
    lat_labels = ["S_60plus", "S_30to60", "S_0to30", "N_0to30", "N_30to60", "N_60plus"]
    support["latitude_band"] = pd.cut(
        support.lat_midpoint, bins=lat_bins, labels=lat_labels, include_lowest=True)
    joined = sp.merge(support, on="key", how="inner", validate="many_to_one")
    assert len(joined) == len(sp)
    joined["delta_sign"] = np.select(
        [joined["class"].eq("old_only"), joined["class"].eq("new_only")],
        [-1, 1], default=0)
    joined["delta_contribution"] = joined.site_influence * joined.delta_sign
    joined["old_only_contribution"] = np.where(
        joined["class"].eq("old_only"), joined.site_influence, 0.0)
    joined["new_only_contribution"] = np.where(
        joined["class"].eq("new_only"), joined.site_influence, 0.0)
    diff_beta = float(models["updated_minus_published"][BETA])
    np.testing.assert_allclose(joined.delta_contribution.sum(), diff_beta, atol=1e-8, rtol=0)
    for col, model in [("old_only_contribution", "old_only"),
                       ("new_only_contribution", "new_only")]:
        np.testing.assert_allclose(joined[col].sum(), models[model][BETA], atol=1e-8, rtol=0)

    unique_species = joined.groupby([ID, "wcvp_accepted_name", "class"], observed=True).agg(
        observed_site_years=("key", "size"),
        years_recorded=("year", "nunique"),
        island_site_years=("land_class", lambda v: int(v.eq("island").sum())),
        mainland_site_years=("land_class", lambda v: int(v.eq("mainland").sum())),
        contribution_to_new_minus_old=("delta_contribution", "sum")
    ).reset_index()
    unique_species["genus"] = unique_species.wcvp_accepted_name.str.split().str[0]
    delta_taxa = unique_species.loc[unique_species["class"].isin(["old_only", "new_only"])].copy()
    sort = delta_taxa.sort_values("contribution_to_new_minus_old")
    columns = [ID, "wcvp_accepted_name", "class", "genus", "observed_site_years",
               "island_site_years", "mainland_site_years", "contribution_to_new_minus_old"]

    genus = delta_taxa.groupby("genus").agg(
        n_changed_species=(ID, "nunique"),
        absolute_site_years=("observed_site_years", "sum"),
        contribution_to_new_minus_old=("contribution_to_new_minus_old", "sum")
    ).reset_index()
    genus.sort_values("contribution_to_new_minus_old", inplace=True)
    # Each attribution column sums to global model beta; does not estimate
    # within-year/within-latitude biological ALAN response.
    rows = []
    for levels in [["year"], ["year", "land_class"], ["latitude_band"], ["latitude_band", "land_class"]]:
        agg = (joined.groupby(levels, observed=True)
               .agg(changed_species_site_years=("delta_sign", lambda x: int(x.ne(0).sum())),
                    contribution_to_new_minus_old=("delta_contribution", "sum"))
               .reset_index())
        rows.append({"stratification": "_by_".join(levels),
                     "nonindependent_additive_contributions": agg.to_dict("records")})
        np.testing.assert_allclose(
            float(agg.contribution_to_new_minus_old.sum()), diff_beta, atol=1e-8, rtol=0)

    by_class = (unique_species.groupby("class", observed=True)
                .agg(species=(ID, "nunique"), site_species_years=("observed_site_years", "sum"),
                     island_site_species_years=("island_site_years", "sum"),
                     mainland_site_species_years=("mainland_site_years", "sum"))
                .reindex(CLASSES, fill_value=0))
    assert by_class.loc["common", "species"] == 612
    assert by_class.loc["old_only", "species"] == 713
    assert by_class.loc["new_only", "species"] == 132
    assert int(by_class.species.sum()) == 8581

    # Compare concentration: are a few ubiquitous species dominating the delta?
    ranked = delta_taxa.assign(abs_contribution=delta_taxa.contribution_to_new_minus_old.abs())
    ranked = ranked.sort_values("abs_contribution", ascending=False)
    top10_abs = float(ranked.abs_contribution.head(10).sum())
    total_abs = float(ranked.abs_contribution.sum())
    return {
        "species_support_by_record_class": by_class.reset_index(names="class").to_dict("records"),
        "source_total_accepted_species": int(by_class.species.sum()),
        "source_total_unique_species_site_years": int(by_class.site_species_years.sum()),
        "top20_species_increasing_delta": sort.tail(20).iloc[::-1][columns].to_dict("records"),
        "top20_species_decreasing_delta": sort.head(20)[columns].to_dict("records"),
        "top15_genera_increasing_delta": genus.tail(15).iloc[::-1].to_dict("records"),
        "top15_genera_decreasing_delta": genus.head(15).to_dict("records"),
        "sum_absolute_taxon_contributions": total_abs,
        "top10_absolute_taxa_share_of_absolute_contributions": top10_abs / total_abs,
        "additive_attribution_by_stratum": rows,
        "note": ("Signed contributions are algebraic changes in a fixed observational "
                 "photography-regression coefficient, NOT local colonization, loss, "
                 "pollinator replacement, causal mediation or species-level ALAN effects.")
    }


def old_only_evidence_sources(old: pd.DataFrame, sp: pd.DataFrame) -> dict:
    """Audit available old SOURCE titles, without assuming they cause GloBI removal."""
    only_names = set(sp.loc[sp["class"].eq("old_only"), "wcvp_accepted_name"])
    a = old.loc[old.guild.eq("moth_family_visit") & old.plant_name.isin(only_names)].copy()
    if a.empty:
        return {"matched_original_positive_rows": 0}
    title = a.study_title.fillna("").str.lower()
    a["inat"] = title.str.contains("inaturalist", regex=False)
    a["trophich"] = title.str.contains("trophich", regex=False)
    studies = a.groupby("plant_name").agg(
        n_source_rows=("plant_name", "size"),
        all_iNaturalist=("inat", "all"),
        any_non_iNaturalist=("inat", lambda v: bool((~v).any())),
        all_trophiCH=("trophich", "all")
    )
    return {
        "old_only_accepted_names_with_exact_old_source_record": len(studies),
        "old_positive_source_rows_for_old_only_names": len(a),
        "old_only_names_with_only_iNaturalist_titled_rows": int(studies.all_iNaturalist.sum()),
        "old_only_names_with_any_non_iNaturalist_titled_row": int(studies.any_non_iNaturalist.sum()),
        "old_only_names_with_only_trophiCH_titled_rows": int(studies.all_trophiCH.sum()),
        "caveat": ("Exact accepted-name records only, not comprehensive synonym provenance; "
                   "title substring is a study-origin proxy. Missing in new GloBI "
                   "does not imply a false old interaction or real biological loss.")
    }


def run(inputs: Path, signed: Path, outfile: Path) -> dict:
    if outfile.exists():
        raise FileExistsError(f"Refusing to overwrite source output: {outfile}")
    s, key_sha, accepted, legacy = primary(inputs, signed)
    assert len(s) == 10698 and s.region5.nunique() == 78 and int(s.island.sum()) == 933
    assert s.region_year.nunique() == 195
    z = classify_species(
        accepted, legacy, ROOT / "results/globi_v011_fiveyear/source_moth_species.csv")
    counts, responses, models = supported_fits(s, z)
    reference = json.loads(
        (ROOT / "results/globi_v011_fiveyear/comparison.json").read_text())
    for name, previous in [("published_old", "published"), ("updated_v011", "updated"),
                           ("updated_minus_published", "updated_minus_published")]:
        np.testing.assert_allclose(
            [models[name][BETA], *models[name]["95CI_cluster_5deg"]],
            [reference["fits"][previous][BETA], *reference["fits"][previous]["95CI_cluster_5deg"]],
            atol=1e-8, rtol=0)
    species = aggregate_influence(s, z, responses, models)
    result = {
        "schema": "ALAN_v011_old713_new132_common612_same10698_site_delta_attribution_v1",
        "status": "AUTHENTIC_SOURCE_RECOMPUTED_EXPLORATORY",
        "parent_frozen_target_hash": key_sha,
        "source_input_sha_manifest": "config/alan_v011_fiveyear_input_sha256.json",
        "source_label_table": "results/globi_v011_fiveyear/source_moth_species.csv",
        "design": {
            "original_same_photo_site_years": 10698,
            "island_site_years": 933, "mainland_site_years": 9765,
            "geographic_clusters": 78, "year_by_region_groups": 195,
            "model": "Same fixed original 2019/2021/2023/2024/2025 source, same 4 sampled VIIRS months, controls, class weights, cluster intervals",
            "no_origin_location_or_species_x_location_rows_public": True,
        },
        "additive_component_fits": models,
        "additive_contribution_diagnostics": species,
        "published_old_only_study_origin_proxy": old_only_evidence_sources(legacy, z),
        "interpretation": {
            "contribution_measure": ("Regression-geometry-weighted species-presence "
                                     "contribution to observed photographic-GloBI evidence "
                                     "island-minus-mainland nightlight coefficient"),
            "no_actual_local_visitor_decrease_proven": True,
            "no_true_moth_absence_inferred_from_missing_source": True,
            "exploratory_and_posthoc_no_confirmatory_causal_effect": True,
            "per_species_contribution_not_species_specific_regression": True,
            "GloBI_interaction_source_disjoint_from_photo_observer_not_proven": True,
        }
    }
    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                  allow_nan=False, default=str) + "\n")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("signed", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    output = run(args.inputs, args.signed, args.output)
    print(json.dumps({
        "status": output["status"],
        "target_hash": output["parent_frozen_target_hash"],
        "fits": {name: {k: obj[k] for k in (BETA, "95CI_cluster_5deg", "p_nominal")}
                 for name, obj in output["additive_component_fits"].items()},
        "species_support": output["additive_contribution_diagnostics"]["species_support_by_record_class"],
    }, ensure_ascii=False, allow_nan=False))
