#!/usr/bin/env python3
"""Henn et al. 2014: paired stream E+T vs other-order composition contrast.

Audit input: *derived* six-order summaries manually summed from the original
Table 1 family rows, not the total row (which disagrees in San Marcos light).
This is one experiment across five physical streams; not six independent
order-level studies and not a flow-standardised drift-density meta effect.
"""
import csv
import itertools
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ORDER_CSV = ROOT / "henn2014_derived_order_by_stream.csv"
STREAM_CSV = ROOT / "henn2014_table1_five_stream_counts.csv"
ORDERS = ("Coleoptera", "Diptera", "Ephemeroptera", "Lepidoptera", "Odonata", "Trichoptera")
RIVERS = (("Barton", "Barton Creek"), ("Bull", "Bull Creek"),
          ("Onion", "Onion Creek"), ("Comal", "Comal River"),
          ("SanMarcos", "San Marcos River"))
EPT = ("Ephemeroptera", "Trichoptera")  # Plecoptera not observed in Table 1
T4_975 = 2.7764451051977987
PUBLISHED_TABLE_TOTALS = {
    "Barton": (57, 7), "Bull": (49, 59),
    "Onion": (46, 98), "Comal": (286, 601), "SanMarcos": (791, 1196),
}
SOURCE_DERIVED_TOTALS = {
    "Barton": (57, 7), "Bull": (49, 59),
    "Onion": (46, 98), "Comal": (286, 601), "SanMarcos": (793, 1196),
}


def parse_order_rows(rows):
    if len(rows) != 6 or set(r["order"] for r in rows) != set(ORDERS):
        raise ValueError("expected exactly six unique original order groups")
    by_order = {r["order"]: {} for r in rows}
    for r in rows:
        order = r["order"]
        for short, _ in RIVERS:
            by_order[order][short] = {}
            for arm in ("lit", "dark"):
                val = int(r[f"{short}_{arm}"])
                if val < 0:
                    raise ValueError("negative family counts invalid")
                by_order[order][short][arm] = val
    return by_order


def get_table_total_counts():
    with STREAM_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    by_name = {r["stream"]: r for r in rows}
    if len(by_name) != 5 or len(rows) != 5:
        raise ValueError("duplicate original stream or not five streams")
    for short, full in RIVERS:
        r = by_name[full]
        if (int(r["lit_count"]), int(r["dark_count"])) != PUBLISHED_TABLE_TOTALS[short]:
            raise ValueError(f"original Table 1 total row mismatch: {full}")
    return by_name


def one_river(ept_lit, ept_dark, other_lit, other_dark):
    if min(ept_lit, ept_dark, other_lit, other_dark) <= 0:
        raise ValueError("this log ratio requires all four original grouped counts positive")
    return math.log(ept_lit / other_lit) - math.log(ept_dark / other_dark)


def exact_sign_flips(block_log_contrasts):
    mu = abs(statistics.mean(block_log_contrasts))
    values = [abs(statistics.mean(s*v for s,v in zip(signs, block_log_contrasts)))
              for signs in itertools.product((-1,1), repeat=len(block_log_contrasts))]
    return sum(v >= mu - 1e-12 for v in values)/len(values)


def analyze(rows):
    by_order = parse_order_rows(rows)
    get_table_total_counts()
    by_stream = []
    logs = []
    for short, full in RIVERS:
        a = sum(by_order[o][short]["lit"] for o in EPT)
        b = sum(by_order[o][short]["dark"] for o in EPT)
        lit = sum(by_order[o][short]["lit"] for o in ORDERS)
        dark = sum(by_order[o][short]["dark"] for o in ORDERS)
        if (lit, dark) != SOURCE_DERIVED_TOTALS[short]:
            raise ValueError(f"order-sum mismatch vs original 36 family rows: {short}")
        other_lit = lit-a
        other_dark = dark-b
        log_ratio = one_river(a,b,other_lit,other_dark)
        by_stream.append({
            "stream":full, "size": "large" if short in ("Comal","SanMarcos") else "small",
            "ET_lit":a, "ET_dark":b,
            "other_lit":other_lit, "other_dark":other_dark,
            "ET_fraction_lit":a/lit, "ET_fraction_dark":b/dark,
            "composition_ratio_of_RRs":math.exp(log_ratio),
            "source_total_lit":lit, "source_total_dark":dark,
            "source_table_total_lit":PUBLISHED_TABLE_TOTALS[short][0],
            "source_table_total_dark":PUBLISHED_TABLE_TOTALS[short][1],
        })
        logs.append(log_ratio)
    mu = statistics.mean(logs)
    se = statistics.stdev(logs)/math.sqrt(5)
    ept_gm = math.exp(mu)
    low,high = math.exp(mu-T4_975*se),math.exp(mu+T4_975*se)

    # The published original family row sum for San Marcos lit is 793,
    # while the published Table 1 total row says 791. We do not alter source.
    # Bound the estimand by assigning exactly two excess individuals entirely
    # to either E+T or other orders; do not call this a corrected source.
    ref = by_stream[-1]
    lows = []
    for remove_ept, remove_other in ((2,0),(0,2)):
        last_log = one_river(ref["ET_lit"]-remove_ept, ref["ET_dark"],
                             ref["other_lit"]-remove_other, ref["other_dark"])
        lows.append(math.exp(statistics.mean(logs[:4]+[last_log])))

    # Four numerically dominant families from author Table 1:
    # Baetidae 209/400, Leptohyphidae 449/587, Chironomidae 174/284,
    # Simuliidae 77/188 (lit/dark respectively). These family totals
    # are used only for a descriptive deficit contribution audit.
    focal_lit = 209+449+174+77
    focal_dark = 400+587+284+188
    all_lit = sum(x["source_total_lit"] for x in by_stream)
    all_dark = sum(x["source_total_dark"] for x in by_stream)
    return {
        "study":"Henn et al. 2014",
        "doi":"10.1080/02705060.2014.900654",
        "original_table":"Table 1",
        "data_level":"36 original family entries manually aggregated to 6 orders by 5 physical streams",
        "ambiguous_source_family_label":"Tipulidae occurs under Coleoptera and Diptera in printed table; original grouping retained without speculative taxonomic repair",
        "authors_reported_family_count":36,
        "printed_family_rows":36,
        "distinct_family_labels":35,
        "exposure_unit":"5 streams, each sampled on distinct light versus dark nights",
        "ET_definition":"Ephemeroptera + Trichoptera; no Plecoptera in author source",
        "stream_results":by_stream,
        "n_independent_stream_pairs":5,
        "estimand":"mean across streams of log[(ETlit/Otherlit)/(ETdark/Otherdark)]",
        "equal_stream_compositional_ratio_of_RRs":ept_gm,
        "descriptive_t4_95_interval":[low,high],
        "se_between_streams_ln_scale":se,
        "posthoc_two_sided_signflip_p":exact_sign_flips(logs),
        "all_five_ratios_above_one":all(v>0 for v in logs),
        "leave_one_stream_out":{
            RIVERS[i][1]:math.exp(statistics.mean(v for j,v in enumerate(logs) if i!=j))
            for i in range(5)
        },
        "sum_source_family_group_counts":{"lit":all_lit,"dark":all_dark,"total":all_lit+all_dark},
        "sum_original_table_total_row":{"lit":1229,"dark":1961,"total":3190},
        "prose_author_total":{"lit":1212,"dark":1979,"claimed_grand_total":3190},
        "family_minus_table_total":{"lit":all_lit-1229,"dark":all_dark-1961},
        "hypothetical_san_marcos_two_count_allocation_bound":[min(lows),max(lows)],
        "four_abundant_families":{
            "names":["Baetidae","Leptohyphidae","Chironomidae","Simuliidae"],
            "lit":focal_lit,"dark":focal_dark,
            "share_of_source_aggregate_count_deficit":(focal_dark-focal_lit)/(all_dark-all_lit),
            "share_of_all_captured_individuals":(focal_dark+focal_lit)/(all_dark+all_lit),
            "warning":"High contribution to the overall deficit largely tracks these families' original abundance, not proof of greater response sensitivity.",
        },
        "meta_programme_increment":0,
        "limitations":[
            "Exploratory post hoc compositional contrast, unadjusted for all analyses tried in this programme",
            "Student-t df4 interval is descriptive; two-sided sign flips assume symmetric exchangeability of five stream contrasts, not confirmed light-allocation randomization",
            "Different lit/dark sampling nights confound within-stream community turnover and abiotic conditions",
            "Within-sample E+T/other odds cancel equal effort multipliers for all taxa but not taxon-dependent detection or drift selectivity",
            "Source mismatch of two San Marcos individuals retained and bounded rather than repaired",
            "Relative-composition change is not evidence of greater E+T absolute abundance, functional feeding, transfer across land-water, or plant reproduction",
        ]
    }


def main():
    with ORDER_CSV.open(newline="",encoding="utf-8") as handle:
        rows=list(csv.DictReader(handle))
    result=analyze(rows)
    path=ROOT/"henn2014_order_composition_20261124.json"
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":"REPLAY_OK",
        "n_physical_streams":result["n_independent_stream_pairs"],
        "ratio_ET_vs_other":result["equal_stream_compositional_ratio_of_RRs"],
        "t4_95":result["descriptive_t4_95_interval"],
        "signflip_p":result["posthoc_two_sided_signflip_p"],
        "source_table_mismatch":result["family_minus_table_total"],
        "meta_increment":result["meta_programme_increment"],
    },ensure_ascii=False))


if __name__=="__main__":
    main()
