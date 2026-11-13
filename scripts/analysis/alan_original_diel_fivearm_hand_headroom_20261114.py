#!/usr/bin/env python3
"""Revisit real four-arm diel reproduction with original HAND fifth arms.

Original 1652 author treatments, original 301 site/year pairs, 153 four-arm
outcomes, authenticated prior-light source. This is a *descriptive* protocol
sensitivity: hand pollination is not a uniformly saturated pollen treatment.
Never infer ALAN effects, moth-mediated compensation, or saturation causality.
Only aggregate source results leave CI; no original plant/site rows exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

from scripts.analysis.alan_original_DvN_reproductive_vs_PRIOR_PCNL_20261113 import (
    load_source, require,
)
from scripts.analysis.alan_original_diel_fourarm_reproductive_function_20261114 import (
    exact_source_arms, qualify, descriptive_study_sign, AUTHOR_STUDY_BLOB,
)


def study_summary(frame: pd.DataFrame, infer: bool = True) -> dict:
    """Independent publication medians, no individual taxon/site output."""
    paper = frame.groupby("study_key", observed=True).agg(
        I=("two_window_nonadditivity", "median"),
        H=("hand_headroom_proxy", "median"),
    )
    n = len(paper)
    out = {
        "papers": int(n),
        "pairs": int(len(frame)),
        "plants": int(frame.species.nunique()),
        "median_study_I": float(paper.I.median()) if n else None,
        "median_study_H": float(paper.H.median()) if n else None,
        "negative_I_papers": int(paper.I.lt(0).sum()),
        "positive_I_papers": int(paper.I.gt(0).sum()),
        "positive_H_papers": int(paper.H.gt(0).sum()),
        "negative_H_papers": int(paper.H.lt(0).sum()),
        "zero_H_papers": int(paper.H.eq(0).sum()),
    }
    if infer and n >= 25 and frame.species.nunique() >= 30:
        v = paper.I.to_numpy(float)
        nn = int(np.sum(v < 0))
        np_ = int(np.sum(v > 0))
        pos_h = int(np.sum(paper.H.to_numpy(float) > 0))
        neg_h = int(np.sum(paper.H.to_numpy(float) < 0))
        out["study_sign_p_nominal"] = float(binomtest(nn, nn + np_, 0.5).pvalue)
        out["hand_headroom_sign_p_nominal"] = float(
            binomtest(pos_h, pos_h + neg_h, 0.5).pvalue)
    elif infer:
        out["study_sign_p_nominal"] = None
        out["hand_headroom_sign_p_nominal"] = None
        out["inferential_gate"] = "NOT_ESTIMABLE_FEW_INDEPENDENT_PAPERS_OR_TAXA"
    else:
        out["inference_not_run_on_outcome_selected_subgroup"] = True
    return out


def describe_hand_subset(four: pd.DataFrame) -> dict:
    """The fifth arm is genuine author hand treatment within exact same keys."""
    h = four.loc[four.hand.notna()].copy()
    require(np.isfinite(h.hand).all() and (h.hand >= 0).all(),
            "Original HAND reproductive means must be numeric nonnegative")
    den = h.hand + h.open
    undefined = int(den.eq(0).sum())
    h = h.loc[den.gt(0)].copy()
    h["hand_headroom_proxy"] = (h.hand - h.open) / (h.hand + h.open)
    require(h.hand_headroom_proxy.between(-1, 1).all(),
            "Hand/open bounded difference corrupted")
    require((len(h), h.study_key.nunique(), h.species.nunique()) == (63, 32, 39),
            "Same-site original HAND fifth-arm 63/32/39 source invariant changed")
    base = descriptive_study_sign(four, "two_window_nonadditivity")
    main = descriptive_study_sign(h, "two_window_nonadditivity")
    require(base["source_original_comparisons"] == 153 and
            base["source_studies"] == 79 and
            base["source_species"] == 85 and
            abs(base["median_original_independent_study_effects"] +
                0.057470648581547865) < 1e-10,
            "Four-arm parent original study baseline changed")
    require(main["source_original_comparisons"] == 63 and
            main["source_studies"] == 32 and
            main["source_species"] == 39,
            "Fifth arm study support changed")
    paper = h.groupby("study_key", observed=True).agg(
        I=("two_window_nonadditivity", "median"),
        H=("hand_headroom_proxy", "median"),
    )
    species = h.groupby("species", observed=True).hand_headroom_proxy.median()
    species_pos = int(species.gt(0).sum())
    species_neg = int(species.lt(0).sum())
    require((len(species), species_pos, species_neg) == (39, 31, 8),
            "Authenticated plant-grain HAND>OPEN source support changed")
    rho, nominal_p = spearmanr(paper.H.to_numpy(), paper.I.to_numpy())
    require(np.isfinite(rho), "Non-identifiable hand/headroom association")
    # Paper-resampled interval is exploratory; shared OPEN in I and H
    # induces mathematical coupling even when visitor biology is unrelated.
    rng = np.random.default_rng(20261114)
    values = paper[["H", "I"]].to_numpy(float)
    boot = []
    for _ in range(3000):
        i = rng.integers(0, len(paper), len(paper))
        rr, _ = spearmanr(values[i, 0], values[i, 1])
        if np.isfinite(rr):
            boot.append(float(rr))
    require(len(boot) >= 2900, "Study bootstrap correlation not supported")
    positive = h.loc[h.hand_headroom_proxy.gt(0)]
    nonpositive = h.loc[h.hand_headroom_proxy.le(0)]
    return {
        "schema": "alan_real301_original_fivearm_hand_reproductive_headroom_v1",
        "status": "SOURCE_AUTHENTIC_EXPLORATORY_FIVE_ARM_REPRODUCTIVE_SENSITIVITY",
        "original_author_treatment_gitblob": AUTHOR_STUDY_BLOB,
        "original_fourarm_reconstructed": {
            "pairs": len(four),
            "studies": int(four.study_key.nunique()),
            "plant_species": int(four.species.nunique()),
            "study_median_I": base["median_original_independent_study_effects"],
            "study_sign_p_nominal": base["source_study_sign_test_two_sided_p_nominal"],
        },
        "fifth_hand_arm_support": {
            "pairs": len(h),
            "studies": int(h.study_key.nunique()),
            "plant_species": int(h.species.nunique()),
            "both_HAND_and_OPEN_zero_not_imputed": undefined,
            "positive_hand_minus_open_pairs": len(positive),
            "nonpositive_hand_minus_open_pairs": len(nonpositive),
            "median_all_pair_headroom_proxy": float(h.hand_headroom_proxy.median()),
            "median_original_study_headroom_proxy": float(paper.H.median()),
            "original_plant_median_headroom_positive": species_pos,
            "original_plant_median_headroom_negative": species_neg,
            "original_plant_median_headroom": float(species.median()),
            "plant_grain_sign_p_nominal_exploratory": float(
                binomtest(species_pos, species_pos + species_neg, .5).pvalue),
            "source_species_phylogenetically_nonindependent": True,
        },
        "new_fivearm_study_grain": {
            "all_original_fifth_arm": study_summary(h),
            "HAND_better_than_OPEN_outcome_condition_descriptive": study_summary(positive, infer=False),
            "HAND_not_better_than_OPEN_outcome_condition_descriptive": study_summary(nonpositive, infer=False),
            "original_I_study_sign_test_nominal": main["source_study_sign_test_two_sided_p_nominal"],
            "study_equal_spearman_headroom_I": float(rho),
            "study_equal_spearman_p_nominal_exploratory": float(nominal_p),
            "study_bootstrap_spearman_95pct_DESCRIPTIVE": [
                float(np.quantile(boot, .025)), float(np.quantile(boot, .975))
            ],
            "study_resamples": 3000,
        },
        "limits": {
            "HAND_is_not_uniform_standardised_pollen_saturation": True,
            "HAND_minus_OPEN_not_unbiased_pollen_limitation_causal_effect": True,
            "I_and_H_both_include_OPEN_outcome_mechanically": True,
            "positive_vs_nonpositive_HAND_subsets_are_outcome_selected": True,
            "source_publications_may_contribute_to_both_subsets": True,
            "original_direct_experiments_not_randomised_by_satellite_ALAN": True,
            "reproductive_fourarm_novelty_is_exploratory_extension_of_Kendall_Nicholson": True,
            "no_local_moth_or_diurnal_pollinator_taxon_assigned": True,
            "only_aggregate_no_identifying_study_species_site_rows_exported": True,
        },
    }


def run(joined: Path, receipt: Path, author: Path, out: Path) -> dict:
    require(not out.exists(), "Never overwrite original real-source results")
    source, prior_receipt = load_source(joined, receipt)
    require((len(source), source.study_key.nunique(), source.species.nunique()) ==
            (301, 118, 123), "Original DvN reproductive source changed")
    all_source, coverage = exact_source_arms(author, source)
    four, eligible = qualify(all_source)
    require((len(four), four.study_key.nunique(), four.species.nunique()) ==
            (153, 79, 85), "True source four-arm support must match parent")
    r = describe_hand_subset(four)
    r["source"] = {
        "paper_doi": "10.1111/ele.70036",
        "original_source_sha256": prior_receipt["joined_file_sha256"],
        "annual_prior_PCNL_MD5_verified": prior_receipt["annual_rasters_MD5_source_verified"],
        "original_true_site_year_DvN_pairs": 301,
        "source_total_original_treatments": coverage["original_source_treatment_rows"],
        "original_4arm_source_eligible": eligible["valid_four_arm_pairs"],
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print("ALAN_ORIGINAL_FIVEARM_HAND_HEADROOM", json.dumps(r, ensure_ascii=False), flush=True)
    return r


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("joined", "receipt", "author", "out"):
        p.add_argument("--" + name, required=True, type=Path)
    z = p.parse_args()
    run(z.joined, z.receipt, z.author, z.out)
