"""Re-extract Henn et al. 2014 Table 1 paired counts with exact river units.

This is NOT a drift-density effect, and cannot be combined directly with the
EERES 2021 flow-regime experiment or Perkin 2014 modeled drift without
homogenizing denominators and verified independence.
"""
import csv
import itertools
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "henn2014_table1_five_stream_counts.csv"
# t(.975, df=4), from scipy.stats.t.ppf(0.975, 4) reference
T95_DF4 = 2.7764451051977987


def aggregate(rows):
    if len(rows) != 5:
        raise ValueError("exactly five original paired streams expected")
    names = [r["stream"] for r in rows]
    if len(set(names)) != 5:
        raise ValueError("duplicate or non-paired stream")
    if sum(r["stream_size"] == "large" for r in rows) != 2:
        raise ValueError("two large-stream physical units expected")
    if sum(r["stream_size"] == "small" for r in rows) != 3:
        raise ValueError("three small-stream physical units expected")
    if any(r["doi"] != "10.1080/02705060.2014.900654" for r in rows):
        raise ValueError("mixed publication provenance")
    arms = []
    for r in rows:
        l = int(r["lit_count"])
        d = int(r["dark_count"])
        if min(l, d) <= 0:
            raise ValueError("zero counts require separate model")
        arms.append({
            "stream": r["stream"], "size": r["stream_size"],
            "lit_count": l, "dark_count": d,
            "lit_over_dark_count_ratio": l / d,
            "ln_ratio": math.log(l / d),
        })

    lit_sum = sum(x["lit_count"] for x in arms)
    dark_sum = sum(x["dark_count"] for x in arms)
    if (lit_sum, dark_sum) != (1229, 1961):
        raise ValueError("published Table 1 totals do not match")
    logs = [x["ln_ratio"] for x in arms]
    mean_ln = statistics.mean(logs)
    se = statistics.stdev(logs) / math.sqrt(5)
    t_stat = mean_ln / se
    # Exact df=4 t p-value needs scipy for formal assessment; avoid
    # introducing dependency merely for an observational sensitivity.
    flip_means = [
        abs(sum(sign * x for sign, x in zip(logs, signs)) / 5)
        for signs in itertools.product((-1, 1), repeat=5)
    ]
    flip_p = sum(x >= abs(mean_ln) - 1e-12 for x in flip_means) / len(flip_means)
    grouped = {}
    for size in ("large", "small"):
        subset = [x for x in arms if x["size"] == size]
        grouped[size] = {
            "n_streams": len(subset),
            "geometric_mean_stream_RR": math.exp(statistics.mean(x["ln_ratio"] for x in subset)),
            "grand_total_count_RR": sum(x["lit_count"] for x in subset) / sum(x["dark_count"] for x in subset),
        }
    return {
        "programme_id": "Henn2014_TexasFive",
        "doi": "10.1080/02705060.2014.900654",
        "authoritative_extraction": "Published author Table 1, original five paired streams, not text-only grand total",
        "publication_text_counts": {"lit": 1212, "dark": 1979, "combined": 3190},
        "table1_count_sums": {"lit": lit_sum, "dark": dark_sum, "combined": lit_sum + dark_sum},
        "table_vs_text_warning": "Table1 sums differ by +17 lit and -18 dark versus prose; prose lit+dark=3191, Table1 combined=3190. Preserve, do not silently reconcile.",
        "stream_level_data": arms,
        "grand_total_lit_over_dark": lit_sum / dark_sum,
        "stream_equal_geometric_mean_RR": math.exp(mean_ln),
        "stream_equal_lnRR": mean_ln,
        "stream_equal_lnRR_se_between_five_streams": se,
        "descriptive_t4_interval_stream_RR": [math.exp(mean_ln-T95_DF4*se), math.exp(mean_ln+T95_DF4*se)],
        "unadjusted_32_signflip_two_sided_p_sensitivity": flip_p,
        "leave_one_stream_out_equal_weight_ratio": {
            r["stream"]: math.exp(statistics.mean(x["ln_ratio"] for x in arms if x["stream"] != r["stream"]))
            for r in arms
        },
        "size_strata_exploratory": grouped,
        "meta_adopted": False,
        "independent_programme_increment": 0,
        "limitations": [
            "Counts per net and two-hour sampling night are not independently flow-volume-standardized drifting densities.",
            "Lit and dark treatments were sampled on different randomized nights in each stream; hydrology and night conditions may differ.",
            "The source paper reports the all-stream paired drift density difference as not significant (p>0.05).",
            "The extreme approximately 1482 lux lamp arrangement is unlike the 0.4-1.2 lux Perkin forest-stream experiment.",
            "Five streams, not ten independent samples or 3190 individual lighting replicates.",
            "These two summary estimands weight streams differently and can diverge even with accurate counts.",
        ],
    }


def main():
    with INPUT.open(encoding="utf-8", newline="") as src:
        result = aggregate(list(csv.DictReader(src)))
    output = ROOT / "henn2014_paired_count_result_20261124.json"
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({
        "status": "PASS",
        "paired_n": len(result["stream_level_data"]),
        "table_rr": result["grand_total_lit_over_dark"],
        "stream_equal_rr": result["stream_equal_geometric_mean_RR"],
        "t4_interval": result["descriptive_t4_interval_stream_RR"],
        "signflip_p": result["unadjusted_32_signflip_two_sided_p_sensitivity"],
        "table_vs_text": result["table_vs_text_warning"],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
