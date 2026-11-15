#!/usr/bin/env python3
"""Source-locked diel function vs floral syndrome (post hoc exploratory B, NOT ALAN A).

Kendall & Nicholson (2025), Ecology Letters doi 10.1111/ele.70036,
source github liamkendall/bumble-bat frozen commit
1945b7c114ad87000608ce6018d0280001c04ef1.

Full 1,652 source rows with matching dark day-only/open, not illumination.
340 genuine reproductive comparison panels from 113 source study IDs.
Original trait table 142 rows, bloom time + tube/open etc.

No trait completeness filter is placed on the original 340-contrast census.
Within each denominator-specific endpoint, pair bloom timing and flower tube
shape, calculate study-equal descriptive medians and paper-bootstrap intervals.
Study bootstrap ignores phylogeny, common authors and original within-study
measurement error and IS NOT causal syndrome effect or ALAN light sensitivity.

Crucial: bloom time, floral tubes, color, odour are only TRAITS and do not
identify effective guilds or light-sensitive moth pollen transfer.
"""
from __future__ import annotations
import collections
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import urllib.request

import numpy as np
import pandas as pd
from scripts.analysis.alan_functional_dependence_20261115 import (
  SOURCE_SHA,TRAIT_SHA,match_panels,join_traits,endpoint_subtype
)

COMMIT="1945b7c114ad87000608ce6018d0280001c04ef1"
UPSTREAM="https://raw.githubusercontent.com/liamkendall/bumble-bat/"+COMMIT+"/data/"
FILES={
 "source":("DvN_LiteratureDataClean_2024-09-13.csv",SOURCE_SHA),
 "traits":("SppTraitData_clean.csv",TRAIT_SHA),
}
SEED=20261115
BOOT=12000
ENDPOINTS=("fruit_proportion","seed_per_fruit","seed_per_ovule")

def verified_bytes(filename,sha):
    req=urllib.request.Request(UPSTREAM+filename,headers={
      "User-Agent":"ALAN-P1-syndrome-functional-dependence-2026-11-15 (+https://github.com/zuizui0223/alan)"
    })
    with urllib.request.urlopen(req,timeout=50) as r:
        b=r.read(2_500_000)
    if hashlib.sha256(b).hexdigest()!=sha:
        raise ValueError("Upstream source SHA mismatch "+filename)
    return b

def classify_bloom(s:str)->str:
    if s=="night":return "night"
    if s=="day":return "day"
    if s=="both":return "both"
    return "unknown"

def classify_shape(s:str)->str:
    if s=="tube":return "tube"
    if s and s.lower() not in ("nan","<na>","unknown",""):
        return "non_tube"
    return "unknown"

def prep(raw,traits):
    p,audit=match_panels(raw)
    q=p[p.treatment_effectiveness_metric.isin(("fruit set","seed set"))].copy()
    if not (len(raw)==1652 and len(q)==340 and q.study_ID.nunique()==113):
        raise AssertionError("Authentic Kendall source count/selection drift")
    j=join_traits(q,traits)
    j["endpoint_subtype"]=[
      endpoint_subtype(m,u) for m,u in zip(j.treatment_effectiveness_metric,j.effectiveness_units)
    ]
    j["bloom"]=j.bloom_period_simple.map(classify_bloom)
    j["shape_bin"]=j.flower_shape.map(classify_shape)
    j["trait_profile"]=[
       (b+"_"+t if b in ("night","day") and t!="unknown" else "other_or_unknown")
       for b,t in zip(j.bloom,j.shape_bin)
    ]
    j["night_loss"]=pd.to_numeric(j.night_loss,errors="raise")
    if not np.isfinite(j.night_loss).all():
        raise AssertionError("Nonfinite physiological night access contrast")
    return j,audit

def per_study(values):
    return values.groupby("study_ID").night_loss.median().to_numpy(float)

def study_summary(q):
    study=per_study(q)
    return {
      "panels":len(q),
      "study_IDs":int(q.study_ID.nunique()),
      "original_taxa":int(q.plant_species.nunique()),
      "study_equal_median":float(np.median(study)) if len(study) else None,
      "study_equal_mean":float(np.mean(study)) if len(study) else None,
      "source_panels_above_20pct":int((q.night_loss>.2).sum()),
      "source_negative_night_losses":int((q.night_loss<0).sum())
    }

def diag_boot(q,idx=0):
    if q.empty:return {"status":"NO_SOURCE_PANEL"}
    keys=sorted(q.study_ID.unique())
    grouped=[q[q.study_ID==k].night_loss.median() for k in keys]
    x=np.array(grouped,dtype=float)
    if len(x)<3:
        return {"status":"SUPPRESSED_LT3_SOURCE_PAPERS",
                "n_studies":len(x)}
    rng=np.random.default_rng(SEED+idx)
    draws=rng.choice(x,(BOOT,len(x)),replace=True)
    qs=np.quantile(np.median(draws,axis=1),[.025,.5,.975])
    return {"status":"DESCRIPTIVE_RESEARCH_PAPER_BOOTSTRAP",
            "n_studies":len(x),
            "interval_2p5_50_97p5":qs.tolist(),
            "not_phylogenetic_or_shared_author_CI":True,
            "not_controlling_plant_sampling_variance":True}

def contrast(q,a,b,index=0):
    x=q[q.trait_profile==a]
    y=q[q.trait_profile==b]
    xx=per_study(x); yy=per_study(y)
    out={
        "comparison":a+"_minus_"+b,
        "a":study_summary(x),"b":study_summary(y),
        "shared_study_count":len(set(x.study_ID)&set(y.study_ID)),
        "comparison_is_mostly_cross_study":len(set(x.study_ID)&set(y.study_ID))<3,
        "median_study_equal_difference":None,
        "paper_bootstrap_95_descriptive":None,
        "status":"UNDERPOWERED_OR_MISSING"
    }
    if not len(xx) or not len(yy):return out
    out["median_study_equal_difference"]=float(np.median(xx)-np.median(yy))
    if min(len(xx),len(yy))<3:return out
    rng=np.random.default_rng(SEED+index)
    aa=rng.choice(xx,(BOOT,len(xx)),replace=True)
    bb=rng.choice(yy,(BOOT,len(yy)),replace=True)
    dif=np.median(aa,axis=1)-np.median(bb,axis=1)
    out["paper_bootstrap_95_descriptive"]=np.quantile(dif,[.025,.975]).tolist()
    out["status"]="DESCRIPTIVE_UNPAIRED_STUDY_BOOTSTRAP"
    return out

def analyze(raw,traits):
    df,audit=prep(raw,traits)
    data={
      "schema":"ALAN_B_20261115_original_flowering_shape_x_night_dependence_v1",
      "source_doi":"10.1111/ele.70036",
      "source_commit":COMMIT,
      "source_SHA256":SOURCE_SHA,
      "traits_SHA256":TRAIT_SHA,
      "original_raw_treatment_rows":len(raw),
      "source_reproductive_panels":len(df),
      "source_studies":int(df.study_ID.nunique()),
      "source_trait_species_rows":len(traits),
      "source_paper_matching_audit":audit,
      "bloom_counts_all_repro":df.bloom.value_counts().to_dict(),
      "shape_counts_all_repro":df.shape_bin.value_counts().to_dict(),
      "floral_syndrome_profile_counts_all_repro":df.trait_profile.value_counts().to_dict(),
      "study_outcomes_by_endpoint":{},
      "global_illumination_reproduction_effect":None,
      "global_ALAN_night_dependence_beta_V":None,
      "warning":"The source is a systematic day/night exclosure synthesis, NOT actual lit/dark ALAN experiments. Differences are between original papers and different taxa; no causal syndrome effect inferred."
    }
    for ix,ep in enumerate(ENDPOINTS):
        q=df[df.endpoint_subtype==ep]
        cats={k:study_summary(q[q.trait_profile==k]) for k in ("night_tube","night_non_tube","day_tube","day_non_tube","other_or_unknown")}
        d={
          "endpoint":ep,"coverage":study_summary(q),
          "profiles":cats,
          "night_tube_vs_day_non_tube":contrast(q,"night_tube","day_non_tube",10*ix),
          "night_tube_vs_day_tube":contrast(q,"night_tube","day_tube",10*ix+1),
          "night_nontube_vs_day_nontube":contrast(q,"night_non_tube","day_non_tube",10*ix+2),
          "study_bootstrap_by_profile":{
             group:diag_boot(q[q.trait_profile==group],100*ix+n)
             for n,group in enumerate(("night_tube","night_non_tube","day_tube","day_non_tube"))
          },
          "all_traits_unknown_kept_in_denominator":True
        }
        data["study_outcomes_by_endpoint"][ep]=d
    if not (len(df[df.endpoint_subtype=="fruit_proportion"])==79 and
            len(df[df.endpoint_subtype=="seed_per_fruit"])==105):
        raise AssertionError("Source-specific denominator stage changed")
    return data

def main():
    with tempfile.TemporaryDirectory() as root:
        root=Path(root)
        for alias,(name,sha) in FILES.items():
            (root/name).write_bytes(verified_bytes(name,sha))
        raw=pd.read_csv(root/FILES["source"][0])
        traits=pd.read_csv(root/FILES["traits"][0])
        out=analyze(raw,traits)
    dest=Path("results/alan_B_P1_global_floral_syndrome_functional_v_20261115.json")
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
      "source":out["source_reproductive_panels"],"papers":out["source_studies"],
      "bloom_counts":out["bloom_counts_all_repro"],
      "shape_counts":out["shape_counts_all_repro"],
      "results":{k:{"profiles":v["profiles"],
                    "night_tube_vs_day_non_tube":v["night_tube_vs_day_non_tube"]}
                for k,v in out["study_outcomes_by_endpoint"].items()},
      "beta_ALAN":None
    },ensure_ascii=False))

if __name__=="__main__":
    main()
