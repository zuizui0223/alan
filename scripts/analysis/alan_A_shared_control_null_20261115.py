#!/usr/bin/env python3
"""Pre-specified *synthetic null* diagnostic for shared dark/open controls.

Never use simulated slopes as observed ALAN results. This shows how reusing
dark-open data to construct both ALAN lnRR and night-access dependence V can
spuriously favor the preselected negative coefficient.
"""
from __future__ import annotations
import json
import math
import pathlib

import numpy as np

SEED = 20261115
V_TRUE = np.linspace(0.15, 0.75, 8)
P_DARK_OPEN = .60
REPLICATES = 50000


def delta_shared_covariance(p_day: float, p_open: float, n_open: int) -> float:
    """delta-method Cov(1 - D/O, log(L/O)) assuming between-arm independence."""
    if not (0 < p_open <= 1 and 0 <= p_day <= 1 and n_open > 0):
        raise ValueError("Probability or n invalid")
    var_o = p_open * (1 - p_open) / n_open
    return -p_day * var_o / p_open ** 3


def slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    xm=x.mean(axis=1, keepdims=True)
    ym=y.mean(axis=1, keepdims=True)
    denom=((x-xm)**2).sum(axis=1)
    return ((x-xm)*(y-ym)).sum(axis=1)/denom


def simulate_one(n: int, seed: int, replicates: int = REPLICATES) -> dict:
    """True ALAN effect is zero for every species. Trials are independent by cell."""
    rng = np.random.default_rng(seed)
    pday = P_DARK_OPEN * (1 - V_TRUE)
    open_dark = rng.binomial(n, P_DARK_OPEN, (replicates, len(V_TRUE)))/n
    day_dark = rng.binomial(n, pday, (replicates, len(V_TRUE)))/n
    open_lit = rng.binomial(n, P_DARK_OPEN, (replicates, len(V_TRUE)))/n
    # Control measurement independent of all three other cells in this arm.
    independent_open_for_v = rng.binomial(n, P_DARK_OPEN, (replicates, len(V_TRUE)))/n
    valid=(open_dark>0).all(axis=1) & (open_lit>0).all(axis=1) & (independent_open_for_v>0).all(axis=1)
    od=open_dark[valid]; dl=day_dark[valid]; ol=open_lit[valid]; other=independent_open_for_v[valid]
    y=np.log(ol/od)
    shared_v=1-dl/od
    independent_v=1-dl/other
    betas_shared=slopes(shared_v,y)
    betas_independent=slopes(independent_v,y)
    def digest(b):
        return {
            "mean_beta": float(np.mean(b)), "median_beta": float(np.median(b)),
            "fraction_negative": float(np.mean(b<0)),
            "mc_interval_2.5_97.5": [float(q) for q in np.quantile(b,[.025,.975])]
        }
    return {
        "n_per_arm_per_species": n, "n_species": len(V_TRUE),
        "replicates_attempted": replicates, "replicates_usable": int(valid.sum()),
        "shared_dark_open": digest(betas_shared),
        "independent_moderator_open": digest(betas_independent),
        "population_p_open": P_DARK_OPEN,
        "true_ALAN_log_RR": 0.0,
        "true_night_access_dependence_range": [float(V_TRUE.min()), float(V_TRUE.max())]
    }


def main() -> None:
    out = {
      "schema": "ALAN_A_shared_control_spurious_slope_null_v1",
      "run_type": "SYNTHETIC_SIMULATION_NOT_OBSERVED_ALAN_DATA",
      "fixed_seed": SEED, "numpy_version": np.__version__,
      "estimands": {
        "y": "ln(lit_open/dark_open)",
        "V_shared": "1-dark_day_only/dark_open",
        "V_independent": "1-dark_day_only/independently_sampled_dark_open",
        "covariance_delta_under_independent_arms": "-p_day*Var(p_dark_open)/(p_dark_open**3)",
        "within_experiment_direct_interaction_DID":
            "ln(lit_open/dark_open)-ln(lit_day_only/dark_day_only), not equivalent to beta_V"
      },
      "synthetic_scenarios": [simulate_one(n, SEED+n) for n in (20, 60, 200)],
      "notes": [
        "The artificial negative slope arises by construction from shared measurement error in the denominator, even though each species has y_true=0.",
        "Independent outside-study V eliminates direct denominator coupling but may introduce transported biological mismatch and measurement error.",
        "Uncertainty must respect the original light-assignment cluster, plot, period, plant and shared contrasts.",
        "Do not estimate meta beta_V until comparable reproduction outcomes, true site-level covariance and V data have been recovered."
      ]
    }
    result=pathlib.Path("results/alan_A_shared_control_null_simulation_20261115.json")
    result.parent.mkdir(parents=True,exist_ok=True)
    result.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"run_type":out["run_type"], "scenarios":out["synthetic_scenarios"]},ensure_ascii=False))


if __name__=="__main__":
    main()
