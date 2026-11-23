#!/usr/bin/env python3
"""Exploratory colour-balanced matched-block contrasts for Sullivan et al. (2019).

Input is the already authenticated *derived* programme summary JSON, not a
replacement for the published raw original data. The three physical habitat
blocks, not dates, station observations, colours, or taxa, are the n=3 units.

No wet/dry, predation-rate, species-fitness, or causal-mediation conclusion.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
from pathlib import Path

LIGHTS = ("LEDI", "LEDH", "HPSI", "HPSH")
ENDPOINTS = ("AquaticEmergent", "TerrestrialArthropod", "WebCount")
T2_CRITICAL_95 = math.sqrt(2 * 0.95 ** 2 / (1 - 0.95 ** 2))


def endpoint_blocks(source: dict, endpoint: str) -> dict[str, list[float]]:
    assert source["provisional_block_contrasts"]["same_three_groups_all_treatments"] == 3
    entry = source["provisional_block_contrasts"]["results"][endpoint]
    assert entry["independent_physical_block_n"] == 3
    data = {}
    for light in LIGHTS:
        item = entry["estimands"][light]
        vals = [float(v) for v in item["raw_block_lnRR"]]
        if len(vals) != 3 or not all(math.isfinite(v) for v in vals):
            raise ValueError(f"invalid n=3 original block values: {endpoint}/{light}")
        mean = statistics.mean(vals)
        se = statistics.stdev(vals) / math.sqrt(3)
        if not math.isclose(mean, float(item["lnRR"]), abs_tol=1e-10):
            raise ValueError(f"mean source mismatch: {endpoint}/{light}")
        if not math.isclose(se, float(item["SE_across_physical_blocks"]), abs_tol=1e-10):
            raise ValueError(f"SE source mismatch: {endpoint}/{light}")
        data[light] = vals
    return data


def summarize_blocks(values: list[float]) -> dict:
    if len(values) != 3:
        raise ValueError("Only three independent physical blocks supported")
    avg = statistics.mean(values)
    se = statistics.stdev(values) / math.sqrt(len(values))
    t_stat = avg / se if se else math.copysign(math.inf, avg)
    # Exact df=2 Student-t identity: 2-sided p = 1-|t|/sqrt(t^2+2)
    t_p = (1 - abs(t_stat) / math.sqrt(t_stat * t_stat + 2)) if math.isfinite(t_stat) else 0
    t_ci = [math.exp(avg - T2_CRITICAL_95 * se), math.exp(avg + T2_CRITICAL_95 * se)]
    # Exploratory sign-flip reference, not a known-verified randomization test.
    null = [abs(sum(s * v for s, v in zip(signs, values)) / 3)
            for signs in itertools.product((-1, 1), repeat=3)]
    exact_p = sum(v >= abs(avg) - 1e-12 for v in null) / len(null)
    return {
        "block_ln_differences": values,
        "block_ratios": [math.exp(x) for x in values],
        "n_independent_habitat_blocks": 3,
        "average_ln_ratio_of_RRs": avg,
        "colour_balanced_ratio_of_RRs": math.exp(avg),
        "se_across_blocks": se,
        "approx_t_df": 2,
        "approx_t95_ratio_interval": t_ci,
        "approx_t_two_sided_p": t_p,
        "block_signs_positive": sum(v > 0 for v in values),
        "block_signflip_two_sided_p": exact_p,
    }


def block_comparison(lhs: dict[str, list[float]], rhs: dict[str, list[float]], colours: tuple[str, ...]) -> dict:
    # A single log contrast per habitat block, averaged across lamp colours.
    # Shared dark stations therefore enter the mean only once per endpoint.
    for c in colours:
        if c not in LIGHTS:
            raise ValueError(f"unknown colour {c}")
    values = [sum(lhs[c][i] - rhs[c][i] for c in colours) / len(colours)
              for i in range(3)]
    return summarize_blocks(values)


def holm2(p1: float, p2: float) -> tuple[float, float]:
    """Descriptive Holm sensitivity for exactly two selected contrasts."""
    a = [p1, p2]
    order = sorted(range(2), key=lambda i: a[i])
    out = [0.0, 0.0]
    last = 0.0
    for rank, i in enumerate(order):
        last = max(last, min(1., a[i] * (2 - rank)))
        out[i] = last
    return tuple(out)


def analyze(source: dict) -> dict:
    if source.get("publication") != "10.1002/eap.1821":
        raise ValueError("wrong study DOI")
    ep = {name: endpoint_blocks(source, name) for name in ENDPOINTS}
    contrasts = {}
    for denominator in ("TerrestrialArthropod", "WebCount"):
        main = block_comparison(ep["AquaticEmergent"], ep[denominator], LIGHTS)
        main["leave_one_colour_out"] = {
            leave: block_comparison(ep["AquaticEmergent"], ep[denominator],
                                    tuple(c for c in LIGHTS if c != leave))
            for leave in LIGHTS
        }
        contrasts[denominator] = main
    a, b = contrasts.values()
    adj = holm2(a["approx_t_two_sided_p"], b["approx_t_two_sided_p"])
    for obj, p in zip((a, b), adj):
        obj["exploratory_two_contrast_Holm_p"] = p
    return {
        "programme_id": "Sullivan2013_Wetland15",
        "publication_doi": "10.1002/eap.1821",
        "original_data_doi": "10.5061/dryad.gb9534p",
        "source_raw_tsv_sha256": "385444ab318ec4b15edd3161c740aa7f4b2d9b973fe4c2ddf6c29cc56cdec484",
        "analysis_date": "2026-11-24",
        "methods": "Log of ratio of lighting/dark response ratios, averaged equally over four colour treatments within each of three independent habitat blocks; SE between the three block means; descriptive t2 95% interval; two-sided 8-way sign flips of block contrasts; 2 exploratory t-test p-values Holm adjusted.",
        "interpretation_boundary": "Post hoc, only one wetland study with three spatial blocks; placement randomization and light leakage unresolved; t2 is not a design-based randomization test; sign flips are a conservative sensitivity not a validated physical randomization distribution. Aquatic emergence trap catch, terrestrial arthropod catch, and spider web counts are distinct endpoints. Not population production, predation rate, wet-dry moderator, or plant-fitness mediation. Colours share the same dark stations and cannot be treated as independent experiments.",
        "comparisons": contrasts,
    }


def main() -> None:
    path = Path(__file__).with_name("dryad_provenance_and_frozen_block_results.json")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    result = analyze(json.loads(args.source.read_text(encoding="utf-8")))
    contents = json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(contents, encoding="utf-8")
    else:
        print(contents, end="")


if __name__ == "__main__":
    main()
