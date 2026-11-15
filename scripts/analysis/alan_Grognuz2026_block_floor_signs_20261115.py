#!/usr/bin/env python3
"""Grognuz & Knop 2026 supplementary table S3 block-aggregated ALAN seed counts.

28 choice blocks are NOT 28 independent light installations; the pots and
lighting setups are repeatedly used. The source was transcribed/checked
against eight published means by earlier ALAN study; this is a new
descriptive sign/floor analysis of the fixed derived block table.
"""
from __future__ import annotations
import csv,json
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/"results/alan_Grognuz2026_tableS3_block_derived_20261115.csv"
OUT=ROOT/"results/alan_Grognuz2026_block_floor_and_signs_20261115.json"
COLS=("mean_dark_no_slug","mean_dark_slug","mean_lit_no_slug","mean_lit_slug")
N=("n_dark_no_slug","n_dark_slug","n_lit_no_slug","n_lit_slug")
def load(path=SRC):
    with path.open(newline="",encoding="utf-8") as f:d=list(csv.DictReader(f))
    if len(d)!=28 or Counter(x["block_id"][0] for x in d)!=Counter({"T":20,"C":8}):
        raise ValueError("Previous author-derived block source drift")
    for x in d:
        for col in COLS+N:x[col]=float(x[col])
        if min(x[c] for c in COLS)<0 or min(x[c] for c in N)<=0:
            raise ValueError("Seeds nonnegative and each group has >0 observed source pots")
    return d
def summary(rows):
    vals=[]
    for x in rows:
        dark=(x["mean_dark_no_slug"]+x["mean_dark_slug"])/2
        light=(x["mean_lit_no_slug"]+x["mean_lit_slug"])/2
        vals.append({"block":x["block_id"],"dark":dark,"lit":light,
                     "delta":light-dark,"dark_zero":dark==0,
                     "lit_zero":light==0})
    if not vals:return {"n":0}
    n=len(vals)
    out={"n_blocks":n,
         "dark_zero_blocks":[z["block"] for z in vals if z["dark_zero"]],
         "light_zero_blocks":[z["block"] for z in vals if z["lit_zero"]],
         "all_block_signs":{"lit_positive":sum(z["delta"]>0 for z in vals),
                            "lit_negative":sum(z["delta"]<0 for z in vals),
                            "zero":sum(z["delta"]==0 for z in vals)},
         "mean_lit_minus_dark_seed_per_target":
             sum(z["delta"] for z in vals)/n,
         "mean_dark_seed_per_target":sum(z["dark"] for z in vals)/n,
         "mean_lit_seed_per_target":sum(z["lit"] for z in vals)/n,
         "control_outcome_selected_not_causal":True}
    nz=[z for z in vals if z["dark"]>0]
    out["dark_nonzero_only_posthoc"]={
        "n_blocks":len(nz),
        "lit_positive":sum(z["delta"]>0 for z in nz),
        "lit_negative":sum(z["delta"]<0 for z in nz),
        "zero":sum(z["delta"]==0 for z in nz),
        "mean_lit_minus_dark_seed_per_target":
           sum(z["delta"] for z in nz)/len(nz) if nz else None}
    # If the source nonnegative dark counts average 0, all observed dark pots
    # in both slug levels must have zero. Then lit - dark >= 0 is a floor fact.
    if any(z["dark_zero"] and z["delta"]<0 for z in vals):
        raise AssertionError("Impossible negative light-minus-dark when measured dark zero")
    return out
def analyze(data):
    out={"schema":"ALAN_Grognuz2026_block_descriptive_floor_signs_v1",
         "source":"Grognuz & Knop 2026 preprint, supplementary Table S3",
         "source_doi":"10.13140/RG.2.2.36822.10563",
         "source_is_derived_block_means_NOT_original_pot_rows":True,
         "total_blocks":len(data),
         "Trifolium_pratense":summary([x for x in data if x["block_id"].startswith("T")]),
         "Centaurea_jacea":summary([x for x in data if x["block_id"].startswith("C")]),
         "limitations":[
           "Same light installation shared across blocks; no valid experiment-level independence",
           "Dark-zero block sign is algebraically nonnegative and cannot prove rescue",
           "Selecting dark-nonzero blocks conditions on observed outcome and can bias",
           "Only original aggregate 4-pot treatment means; individual zero counts unobserved",
           "No nighttime target visits: plants bagged during night-light exposure",
           "The 2 species are in ONE preprint, not two independent studies",
           "Not a primary ALAN×same-population nocturnal dependence meta regression"]
    }
    return out
def main():
    out=analyze(load())
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("ALAN_GROGNUZ_SOURCE_FLOOR="+json.dumps(out,ensure_ascii=False))
if __name__=="__main__":main()
