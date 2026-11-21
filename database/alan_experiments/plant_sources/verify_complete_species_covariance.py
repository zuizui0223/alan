"""Independent source-table/Jacobian verification of the Penryn four-species analysis.

Usage:
    python verify_penryn_complete_species.py --analysis-dir PATH [--out JSON]

Requires numpy, pandas and scipy. Does not import or execute the analysis
implementation. Source tables are read from PATH's parent directory. This is
an AI calculation check, not human duplicate extraction of the original paper.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


SPECIES = [
    "Agrostis tenuis",
    "Anthoxanthum odoratum",
    "Hypochaeris radicata",
    "Lotus pedunculatus",
]
FOCAL = "Lotus pedunculatus"
TROPHIC = ["T1", "T2", "T3"]
LIGHTS = ["L1", "L2", "L3"]
ATOL = 1e-12
RTOL = 1e-12


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_numeric(actual, stored, label):
    actual = np.asarray(actual, dtype=float)
    stored = np.asarray(stored, dtype=float)
    require(actual.shape == stored.shape, f"{label}: shape mismatch")
    require(np.isfinite(actual).all() and np.isfinite(stored).all(),
            f"{label}: nonfinite value")
    require(np.allclose(actual, stored, atol=ATOL, rtol=RTOL),
            f"{label}: numerical mismatch")
    return float(np.max(np.abs(actual - stored)))


def read_covariance(path, ids):
    frame = pd.read_csv(path)
    require(not frame.duplicated(["effect_id_i", "effect_id_j"]).any(),
            f"{path.name}: duplicate covariance cells")
    require(len(frame) == len(ids) ** 2, f"{path.name}: covariance cell count")
    require(set(frame.effect_id_i) == set(ids) == set(frame.effect_id_j),
            f"{path.name}: unexpected effect IDs")
    return frame.pivot(index="effect_id_i", columns="effect_id_j",
                       values="covariance").loc[ids, ids].to_numpy()


def verify(analysis_dir):
    analysis_dir = analysis_dir.resolve()
    base = analysis_dir.parent
    source_names = ["species_biomass_source_rows.csv",
                    "species_survival_reconciliation.csv",
                    "anic2022_community_units.csv"]
    output_names = ["selected_source_rows.csv", "species_effects.csv",
                    "species_covariance.csv", "contrast_weights.csv",
                    "composition_and_trophic_contrasts.csv", "contrast_covariance.csv",
                    "mean_log_estimand_sensitivity.csv", "comparator_sensitivity.csv",
                    "analysis_summary.json"]
    d, survival, allocation = [pd.read_csv(base / name) for name in source_names]
    require(not d.duplicated(["species", "mesocosm"]).any(),
            "Duplicate source species-by-mesocosm rows")
    require(len(allocation) == 54 and allocation.mesocosm.is_unique,
            "Expected 54 unique allocation records")
    cell_sizes = allocation.groupby(["trophic", "light"]).size()
    require(set(cell_sizes.index) == {(tr, li) for tr in TROPHIC for li in LIGHTS}
            and cell_sizes.eq(6).all(), "Expected nine cells of six mesocosms")
    require(not survival.duplicated(["species", "light", "trophic"]).any(),
            "Duplicate source survival-table cells")

    selected = []
    for species, group in d.groupby("species"):
        s = survival[survival.species.eq(species)]
        complete = (len(group) == 54
                    and set(group.mesocosm) == set(allocation.mesocosm)
                    and np.isfinite(group.total_g).all() and group.total_g.gt(0).all())
        agrees = (len(s) == 9 and s.status.eq("agrees").all()
                  and s.surviving_mesocosms_S1.eq(6).all()
                  and set(zip(s.trophic, s.light)) == set(cell_sizes.index))
        if complete and agrees and not group.component_inconsistency.any():
            selected.append(species)
    require(selected == SPECIES, f"Eligible species changed: {selected}")
    units = d[d.species.isin(selected)].copy()
    require(len(units) == 216, "Expected 216 complete-species source rows")
    joined = units.merge(allocation, on="mesocosm", validate="many_to_one",
                         suffixes=("", "_allocation"), how="left", indicator=True)
    require(joined["_merge"].eq("both").all()
            and joined.light.eq(joined.light_allocation).all()
            and joined.trophic.eq(joined.trophic_allocation).all(),
            "Source and allocation treatment labels disagree")
    selected_saved = pd.read_csv(analysis_dir / "selected_source_rows.csv")
    require(len(selected_saved) == 216
            and not selected_saved.duplicated(["species", "mesocosm"]).any(),
            "Saved selection has missing or duplicate rows")
    original_index = units.set_index(["species", "mesocosm"]).sort_index()
    saved_index = selected_saved.set_index(["species", "mesocosm"]).sort_index()
    require(original_index.index.equals(saved_index.index), "Saved selection keys differ")
    selected_row_error = check_numeric(original_index.total_g, saved_index.total_g,
                                       "Selected source values")

    # Independent route: 36 arithmetic arm means and nine separate 4x4
    # covariance blocks. No Gram/influence-score code is used here.
    columns = [(tr, li, species) for tr in TROPHIC for li in LIGHTS for species in SPECIES]
    means = []
    arm_covariance = np.zeros((36, 36))
    for block, (tr, li) in enumerate((tr, li) for tr in TROPHIC for li in LIGHTS):
        wide = units[units.trophic.eq(tr) & units.light.eq(li)].pivot(
            index="mesocosm", columns="species", values="total_g")[SPECIES]
        require(wide.shape == (6, 4) and wide.notna().all().all(),
                f"Incomplete source block {tr}/{li}")
        means.extend(wide.mean().tolist())
        start = block * 4
        arm_covariance[start:start + 4, start:start + 4] = wide.cov().to_numpy() / 6
    means = np.asarray(means)
    effects = pd.read_csv(analysis_dir / "species_effects.csv")
    require(len(effects) == 24 and effects.effect_id.is_unique,
            "Expected 24 unique species effects")
    expected_keys = {(sp, tr, li) for sp in SPECIES for tr in TROPHIC for li in LIGHTS[1:]}
    require(set(zip(effects.species, effects.trophic, effects.light)) == expected_keys,
            "Unexpected species-effect support")
    log_transform = np.zeros((24, 36))
    for i, row in enumerate(effects.itertuples()):
        log_transform[i, columns.index((row.trophic, row.light, row.species))] = 1
        log_transform[i, columns.index((row.trophic, "L1", row.species))] = -1
    theta = log_transform @ np.log(means)
    jacobian = log_transform / means
    covariance = jacobian @ arm_covariance @ jacobian.T
    stored_covariance = read_covariance(analysis_dir / "species_covariance.csv",
                                        effects.effect_id.tolist())
    errors = {
        "selected_source_biomass": selected_row_error,
        "species_log_effects": check_numeric(theta, effects.lnRR, "Species lnRR"),
        "species_covariance": check_numeric(covariance, stored_covariance, "Species covariance"),
        "species_variance": check_numeric(np.diag(covariance), effects.sampling_variance,
                                           "Species variances"),
    }

    comparators = [species for species in SPECIES if species != FOCAL]

    def cell_vector(trophic, light, compare=comparators):
        vector = np.zeros(len(effects))
        for i, row in enumerate(effects.itertuples()):
            if row.trophic == trophic and row.light == light:
                vector[i] = (1.0 if row.species == FOCAL else
                             -1 / len(compare) if row.species in compare else 0.0)
        return vector

    derived = pd.read_csv(analysis_dir / "composition_and_trophic_contrasts.csv")
    require(len(derived) == 12 and derived.contrast_id.is_unique,
            "Expected 12 unique derived contrasts")
    require(derived.kind.value_counts().to_dict() == {"cell": 6, "modifier": 4, "marginal": 2},
            "Unexpected derived contrast counts")
    contrast_matrix = np.zeros((12, 24))
    expected_families = {"cell": 6, "marginal": 2, "modifier": 4}
    for i, row in enumerate(derived.itertuples()):
        require(row.family_size == expected_families[row.kind], "Unexpected correction family")
        if row.kind == "cell":
            contrast_matrix[i] = cell_vector(row.trophic, row.light)
        elif row.kind == "marginal":
            contrast_matrix[i] = sum(cell_vector(tr, row.light) for tr in TROPHIC) / 3
        else:
            high, low = row.trophic.split("_minus_")
            require((high, low) in [("T2", "T1"), ("T3", "T2")], "Unexpected modifier")
            contrast_matrix[i] = cell_vector(high, row.light) - cell_vector(low, row.light)
    stored_weights = pd.read_csv(analysis_dir / "contrast_weights.csv", index_col=0)
    stored_weights = stored_weights.loc[derived.contrast_id, effects.effect_id].to_numpy()
    errors["contrast_weights"] = check_numeric(contrast_matrix, stored_weights, "Contrast weights")
    estimates = contrast_matrix @ theta
    derived_covariance = contrast_matrix @ covariance @ contrast_matrix.T
    stored_derived_cov = read_covariance(analysis_dir / "contrast_covariance.csv",
                                         derived.contrast_id.tolist())
    errors["derived_log_effects"] = check_numeric(estimates, derived.estimate, "Derived estimates")
    errors["derived_covariance"] = check_numeric(derived_covariance, stored_derived_cov,
                                                 "Derived covariance")
    errors["derived_variance"] = check_numeric(np.diag(derived_covariance),
                                               derived.sampling_variance, "Derived variance")
    standard_error = np.sqrt(np.diag(derived_covariance))
    radius = t.ppf(0.975, 5) * standard_error
    corrected_radius = t.ppf(1 - 0.05 / (2 * derived.family_size.to_numpy()), 5) * standard_error
    for actual, label in [(np.exp(estimates), "ratio"),
                          (np.exp(estimates - radius), "lower_95_t5"),
                          (np.exp(estimates + radius), "upper_95_t5"),
                          (np.exp(estimates - corrected_radius), "lower_95_family_t5"),
                          (np.exp(estimates + corrected_radius), "upper_95_family_t5")]:
        errors[label] = check_numeric(actual, derived[label], label)

    # These are existing requested sensitivities, not additional analyses.
    log_means = units.assign(logmass=np.log(units.total_g)).groupby(
        ["trophic", "light", "species"]).logmass.mean()
    mean_log = []
    for tr in TROPHIC:
        for li in LIGHTS[1:]:
            value = (log_means[tr, li, FOCAL] - log_means[tr, "L1", FOCAL]
                     - sum(log_means[tr, li, sp] - log_means[tr, "L1", sp]
                           for sp in comparators) / 3)
            mean_log.append(dict(trophic=tr, light=li, estimate=float(value)))
    for li in LIGHTS[1:]:
        value = np.mean([r["estimate"] for r in mean_log if r["light"] == li])
        mean_log.append(dict(trophic="equal_weight_T1_T2_T3", light=li, estimate=float(value)))
    mean_log = pd.DataFrame(mean_log).set_index(["trophic", "light"])
    mean_log_stored = pd.read_csv(analysis_dir / "mean_log_estimand_sensitivity.csv").set_index(
        ["trophic", "light"])
    errors["mean_log_estimand"] = check_numeric(mean_log.estimate,
        mean_log_stored.loc[mean_log.index, "estimate"], "Mean-log alternative estimand")
    comparator_stored = pd.read_csv(analysis_dir / "comparator_sensitivity.csv")
    single = []
    single_errors = []
    for sp in comparators:
        for li in LIGHTS[1:]:
            vector = sum(cell_vector(tr, li, [sp]) for tr in TROPHIC) / 3
            estimate = float(vector @ theta)
            se = float(np.sqrt(vector @ covariance @ vector))
            stored = comparator_stored[(comparator_stored["mode"] == "single_comparator")
                & comparator_stored.named_comparator.eq(sp)
                & comparator_stored.kind.eq("marginal") & comparator_stored.light.eq(li)]
            require(len(stored) == 1, "Missing or duplicated single-comparator sensitivity")
            single_errors.append(check_numeric([estimate, se],
                [stored.iloc[0].estimate, stored.iloc[0].se], "Single-comparator sensitivity"))
            single.append(dict(comparator=sp, light=li, ratio=float(np.exp(estimate)),
                lower_95_t5=float(np.exp(estimate - t.ppf(.975, 5) * se)),
                upper_95_t5=float(np.exp(estimate + t.ppf(.975, 5) * se))))
    errors["single_comparator_marginals"] = max(single_errors)

    source_hashes = {name: sha256(base / name) for name in source_names}
    saved_audit = json.loads((analysis_dir / "analysis_summary.json").read_text())
    require(saved_audit["input_sha256"] == source_hashes, "Recorded source hashes differ")
    marginal = derived[derived.kind.eq("marginal")]
    return {
        "schema": "penryn_independent_jacobian_verification_v1",
        "status": "passed",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_dir": str(analysis_dir),
        "method": "36 arithmetic arm means; nine 4x4 sample covariance blocks divided by 6; 24x36 log-ratio Jacobian; independent 12x24 contrast matrix",
        "analysis_implementation_imported_or_executed": False,
        "verification_agent": "AI independent calculation review",
        "human_duplicate_extraction": False,
        "scope": "Repository source tables and saved calculations; this script does not re-extract the original PDF or verify trait values",
        "selected_species": selected,
        "source_rows": len(units), "physical_mesocosms": int(units.mesocosm.nunique()),
        "independent_programmes": 1,
        "species_effects_checked": len(effects), "species_covariance_cells_checked": int(covariance.size),
        "derived_contrasts_checked": len(derived), "derived_covariance_cells_checked": int(derived_covariance.size),
        "species_covariance_rank": int(np.linalg.matrix_rank(covariance)),
        "species_covariance_min_eigenvalue": float(np.linalg.eigvalsh(covariance)[0]),
        "absolute_tolerance": ATOL, "relative_tolerance": RTOL,
        "max_absolute_errors": errors,
        "source_sha256": source_hashes,
        "checked_output_sha256": {name: sha256(analysis_dir / name) for name in output_names},
        "verifier_sha256": sha256(Path(__file__)),
        "marginal_results": marginal[["light", "ratio", "lower_95_t5", "upper_95_t5",
            "lower_95_family_t5", "upper_95_family_t5"]].to_dict("records"),
        "all_four_modifier_nominal_intervals_include_one": bool(
            ((derived.loc[derived.kind.eq("modifier"), "lower_95_t5"] <= 1)
             & (derived.loc[derived.kind.eq("modifier"), "upper_95_t5"] >= 1)).all()),
        "mean_log_alternative_marginal_ratios": {
            li: float(np.exp(mean_log.loc[("equal_weight_T1_T2_T3", li), "estimate"]))
            for li in LIGHTS[1:]},
        "single_comparator_marginals": single,
        "limitations": ["First-order delta and t5 intervals remain approximate with six units per cell",
            "Post hoc complete-survivor species subset; not a general moisture effect or new independent study",
            "Sensitivity point directions do not establish significant individual comparator differences"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).with_name("independent_verification.json"))
    args = parser.parse_args()
    result = verify(args.analysis_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "output": str(args.out.resolve()),
        "species_effects": result["species_effects_checked"],
        "species_covariance_cells": result["species_covariance_cells_checked"],
        "derived_contrasts": result["derived_contrasts_checked"],
        "derived_covariance_cells": result["derived_covariance_cells_checked"],
        "max_absolute_errors": result["max_absolute_errors"]}, indent=2))


if __name__ == "__main__":
    main()
