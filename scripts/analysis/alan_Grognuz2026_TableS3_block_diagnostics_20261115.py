#!/usr/bin/env python3
"""Non-raw derived block-level reanalysis of Grognuz & Knop 2026 Table S3.

Source: Grognuz & Knop, March 2026 PUBLIC PREPRINT DOI
10.13140/RG.2.2.36822.10563. Author Table S3 contains 28 choice blocks,
16 target pots/block, visit counts and parenthesized seed counts.

The 28x4 block means below were independently transcribed from that table and
cross-checked against the author-published 8 treatment-level seed and visitor
means. The individual 448 source observations are NOT redistributed.

LIMITATION: same live potted plants were reused across calendar weeks;
lighting was installed at the same field site for repeated exposure groups.
28 pollinator-choice blocks are NOT independent light-installation replicates.
Block bootstrap is descriptive uncertainty ONLY, not an intervention-cluster
95% CI or newly certified causal p value. Published fitted β is preserved as
an author-reported parameter, not replaced by the aggregates.
"""
from __future__ import annotations
import csv
import json
import math
import pathlib
import random
from collections import defaultdict

ROOT=pathlib.Path(__file__).resolve().parents[2]
SOURCE=ROOT/"results/alan_Grognuz2026_tableS3_block_derived_20261115.csv"
SEED=20261115
REPS=30000
PUBLISHED_POT_MEANS={
    "Trifolium":(4.56,5.34,3.66,3.39),
    "Centaurea":(1.71,1.53,.41,1.85)
}
LABELS=("dark_no_slug","dark_slug","lit_no_slug","lit_slug")
GROUPS={"Trifolium":20,"Centaurea":8}


def load_source(path=SOURCE):
    result=defaultdict(list)
    with path.open(newline="",encoding="utf-8") as handle:
        rows=list(csv.DictReader(handle))
    for row in rows:
        block=row["block_id"]
        species="Trifolium" if block.startswith("T") else ("Centaurea" if block.startswith("C") else None)
        if not species: raise ValueError("Unrecognized block id")
        means=[float(row["mean_"+s]) for s in LABELS]
        ns=[int(row["n_"+s]) for s in LABELS]
        if not all(math.isfinite(x) and x>=0 for x in means):raise ValueError("Invalid seed mean")
        if not all(1<=n<=4 for n in ns):raise ValueError("Invalid number of measured pots")
        result[species].append({"block":block,"means":means,"n":ns})
    for name,n in GROUPS.items():
        if len(result[name])!=n:raise AssertionError("Unexpected published block count")
    all_obs=sum(sum(r["n"]) for group in result.values() for r in group)
    if all_obs!=425:raise AssertionError("Unexpected seed-observed pot count")
    for name,group in result.items():
        for a,expected in enumerate(PUBLISHED_POT_MEANS[name]):
            actual=sum(r["means"][a]*r["n"][a] for r in group)/sum(r["n"][a] for r in group)
            if abs(actual-expected)>.011:
                raise AssertionError("Treatment mean failed to reproduce author published Table S3: "+name+" "+str(a))
    return result


def extract(group):
    means=[sum(r["means"][a] for r in group)/len(group) for a in range(4)]
    obs=[sum(r["n"][a] for r in group) for a in range(4)]
    return means,obs


def effects(means):
    dark=(means[0]+means[1])/2
    lit=(means[2]+means[3])/2
    return {
       "mean_dark":dark,"mean_lit":lit,
       "delta_lit_minus_dark":lit-dark,
       "light_dark_ratio":None if dark==0 else lit/dark,
       "no_slug_ratio":None if means[0]==0 else means[2]/means[0],
       "slug_ratio":None if means[1]==0 else means[3]/means[1],
       "slug_x_light_difference_of_differences":(means[3]-means[1])-(means[2]-means[0])
    }


def percentile(sorted_values,p):
    k=(len(sorted_values)-1)*p
    low=int(k);high=min(low+1,len(sorted_values)-1)
    return sorted_values[low]*(high-k)+sorted_values[high]*(k-low)


def block_boot(group,seed=SEED,reps=REPS):
    rng=random.Random(seed)
    names=("delta_lit_minus_dark","light_dark_ratio","slug_x_light_difference_of_differences")
    vals={name:[] for name in names}
    for _ in range(reps):
        sample=[group[rng.randrange(len(group))] for __ in range(len(group))]
        e=effects(extract(sample)[0])
        for name in names:
            if e[name] is not None and math.isfinite(e[name]):vals[name].append(e[name])
    out={}
    for name,x in vals.items():
        x.sort()
        out[name]={
           "n_valid":len(x),
           "percentile_2p5_50_97p5":None if not x else [percentile(x,p) for p in (.025,.5,.975)]
        }
    out["fraction_delta_negative"]=sum(x<0 for x in vals["delta_lit_minus_dark"])/len(vals["delta_lit_minus_dark"])
    return out


def one_species(name,group,reps):
    means,nobs=extract(group)
    point=effects(means)
    loo=[effects(extract([g for j,g in enumerate(group) if i!=j])[0])["delta_lit_minus_dark"] for i in range(len(group))]
    missing=[len(group)*4-x for x in nobs]
    # Exact pot totals require multiplying each BLOCK cell mean by that
    # cell's observed n first. Grand block-mean * grand n is WRONG when
    # missingness causes heterogeneous n across blocks.
    dark_total=sum(r["means"][a]*r["n"][a] for r in group for a in (0,1))
    lit_total=sum(r["means"][a]*r["n"][a] for r in group for a in (2,3))
    n_dark_miss=sum(missing[a] for a in (0,1))
    n_lit_miss=sum(missing[a] for a in (2,3))
    # All missing dark seeds=0, missing lit seeds=x:
    # both arms originally have 4 pots per factorial cell per block, hence
    # balanced total 8*N pots in lit and dark sides; x at sign reversal:
    x_tip=None if n_lit_miss==0 else (dark_total-lit_total)/n_lit_miss
    return {
       "species":name,"n_choice_blocks":len(group),
       "n_target_pots_total":len(group)*16,
       "n_target_pots_seed_observed":sum(nobs),
       "missing_seeds_by_cell":dict(zip(LABELS,missing)),
       "n_seed_observed_by_cell":dict(zip(LABELS,nobs)),
       "seed_per_pot_cell_means_equal_block":dict(zip(LABELS,means)),
       "effect":point,
       "leave_one_choice_block_delta_range":[min(loo),max(loo)],
       "descriptive_block_resample":block_boot(group,SEED+len(group),reps),
       "tipping_impute_missing_lit_pots_seeds_each_if_dark_missing_assumed_zero":x_tip,
       "tip_assumption":"This is a deliberately favorable-to-light missing-data stress test, NOT a plausible missing-at-random model. No missing seed count is silently set to zero.",
       "actual_independent_lamp_installations":None,
       "source_has_plant_reuse_across_weeks":True,
       "interpretation_scope":"Do not report resample percentile ranges as causal ALAN treatment CIs. The published original negative seed coefficient was a model parameter, not the block mean ratio."
    }


def run():
    source=load_source()
    results=[one_species(name,source[name],REPS) for name in ("Trifolium","Centaurea")]
    out={
       "schema":"ALAN_Grognuz2026_public_TableS3_block_derived_audit_v1",
       "source":"https://www.researchgate.net/publication/401424487_Artificial_light_at_night_restructures_daytime_plant-pollinator_interactions_through_direct_and_herbivore-mediated_pathways",
       "source_public_preprint_doi":"10.13140/RG.2.2.36822.10563",
       "source_author_open_data_doi":"10.17632/drhzznh4ps.1",
       "status":"AUTHOR_TABLE_S3_SECONDARY_MANUAL_DATA_EXTRACTION_VERIFIED_AGAINST_PUBLISHED_TREATMENT_MEANS_NOT_GENUINE_SOURCE_CSV",
       "n_blocks":28,"n_pots_total":448,"n_pots_seed_observed":425,"n_seed_missing":23,
       "four_factorial_cells":LABELS,
       "source_author_seed_model_beta":-.40,
       "source_author_seed_beta_ci95":[-.73,-.08],
       "source_author_seed_beta_p":.014,
       "seed_model_beta_unit":"AUTHOR log-linked model coefficient; not a log-response-ratio from raw pot means",
       "bootstrap_resamples":REPS,
       "bootstrap_warning":"Repeated plants and shared lighting site mean independent treatment replication is smaller than 28 blocks; descriptive bootstrap alone is NOT a valid causal confidence interval.",
       "results":results
    }
    target=ROOT/"results/alan_Grognuz2026_tableS3_block_derived_diagnostics_20261115.json"
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"n_blocks":out["n_blocks"],"n_pots_seed_observed":out["n_pots_seed_observed"],
                      "n_missing":out["n_seed_missing"],
                      "results":[{"species":r["species"],"ratio":r["effect"]["light_dark_ratio"],
                                  "delta":r["effect"]["delta_lit_minus_dark"],
                                  "ci_descriptive":r["descriptive_block_resample"]["delta_lit_minus_dark"]["percentile_2p5_50_97p5"],
                                  "interaction":r["effect"]["slug_x_light_difference_of_differences"],
                                  "tipping_missing":r["tipping_impute_missing_lit_pots_seeds_each_if_dark_missing_assumed_zero"]}
                                 for r in results]},ensure_ascii=False))


if __name__=="__main__":
    run()
