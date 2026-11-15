#!/usr/bin/env python3
"""Independent author-source comparison of nighttime access at fruit vs seed stages.

Kendall & Nicholson, Ecology Letters (2025), doi:10.1111/ele.70036.
Same study + plant + site + year + month + supplementary treatment required.
All outputs exploratory; NOT an artificial-light experiment or a causal meta effect.
"""
from __future__ import annotations
import hashlib, json, math, urllib.request
from collections import defaultdict
from io import BytesIO
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_functional_dependence_20261115 import (
    SOURCE_SHA, KEYS, match_panels, endpoint_subtype
)
ROOT=Path(__file__).resolve().parents[2]
COMMIT="1945b7c114ad87000608ce6018d0280001c04ef1"
URL="https://raw.githubusercontent.com/liamkendall/bumble-bat/"+COMMIT+"/data/DvN_LiteratureDataClean_2024-09-13.csv"
OUT=ROOT/"results/alan_B_global340_original_paired_reproductive_stages_20261115.json"
EXACT=[x for x in KEYS if x not in ("treatment_effectiveness_metric","effectiveness_units","treatmenteffectiveness_metric_comments")]
COARSE=["study_ID","plant_species","country","locality","year_start","additional_treatment"]
FRUIT="fruit_proportion"
SEED=("seed_per_fruit","seed_per_ovule","seed_per_flower","seed_per_plant","seed_formation_probability")

def source():
    req=urllib.request.Request(URL,headers={"User-Agent":"ALAN-20261115-source-stage-gate"})
    with urllib.request.urlopen(req,timeout=75) as r:b=r.read(2500000)
    if hashlib.sha256(b).hexdigest()!=SOURCE_SHA:raise ValueError("Author SHA does not match")
    return pd.read_csv(BytesIO(b))

def get_panels(raw):
    allp,audit=match_panels(raw)
    p=allp[allp.treatment_effectiveness_metric.isin(["fruit set","seed set"])].copy()
    if (len(raw),len(p),p.study_ID.nunique())!=(1652,340,113):
        raise ValueError("Original 1652/340/113 source drift")
    p["stage"]=[endpoint_subtype(m,u) for m,u in zip(p.treatment_effectiveness_metric,p.effectiveness_units)]
    return p,audit

def exact_pairs(p,keys=EXACT):
    needed=set(keys+["study_ID","plant_species","stage","day_mean","open_mean","effectiveness_units"])
    if not needed<=set(p.columns):raise AssertionError("Missing source columns")
    subset=p[p.stage.isin((FRUIT,)+SEED)]
    pairs=[];bad=[];amb=[];candidate=0
    for group,g in subset.groupby(keys,dropna=False,sort=True):
        if not (g.stage==FRUIT).any() or not g.stage.isin(SEED).any():continue
        candidate+=1
        ff=g[g.stage==FRUIT]
        if len(ff)!=1:
            amb.append({"study":str(g.iloc[0].study_ID),"reason":"MULTIPLE_FRUIT_PER_KEY"});continue
        f=ff.iloc[0]
        for stage in SEED:
            ss=g[g.stage==stage]
            if len(ss)>1:
                amb.append({"study":str(f.study_ID),"reason":"MULTIPLE_SEED_PER_KEY","stage":stage})
                continue
            if ss.empty:continue
            s=ss.iloc[0]
            v=[float(f.day_mean),float(f.open_mean),float(s.day_mean),float(s.open_mean)]
            if any(not math.isfinite(t) for t in v) or min(v)<0 or v[1]<=0 or v[3]<=0:
                bad.append({"study":str(f.study_ID),"reason":"INVALID_ARM_MEAN","stage":stage});continue
            if max(v[:2])>1:
                bad.append({"study":str(f.study_ID),"reason":"FRUIT_PROPORTION_GT_ONE","stage":stage});continue
            vf=1-v[0]/v[1]
            vs=1-v[2]/v[3]
            pairs.append({
                "study_id":str(f.study_ID),"plant_species":str(f.plant_species),
                "seed_stage":stage,"V_fruit":vf,"V_seed":vs,
                "delta_seed_minus_fruit":vs-vf,
                "log_DoverO_seed_minus_fruit":math.log(v[2]/v[3])-math.log(v[0]/v[1])
                    if min(v[0],v[2])>0 else None,
                "known_lighting_pooled":str(f.study_ID)=="Macgregor_2019",
                "same_plants_proven":False,"plant_level_SE_available":False,
                "units_fruit":str(f.effectiveness_units),"units_seed":str(s.effectiveness_units)
            })
    return pairs,{"candidate_group_count":candidate,"ambiguities":amb,"invalid":bad,"stage_rows":len(subset)}

def summarize(pairs):
    out={}
    for i,stage in enumerate(SEED):
        qs=[p for p in pairs if p["seed_stage"]==stage and not p["known_lighting_pooled"]]
        groups=defaultdict(list)
        for p in qs:groups[p["study_id"]].append(p["delta_seed_minus_fruit"])
        study={k:float(np.median(v)) for k,v in sorted(groups.items())}
        a=np.array(list(study.values()),float)
        positive=int((a>0).sum());negative=int((a<0).sum())
        # Directional probability is posthoc descriptive only, not a fitted effect/SE.
        n=positive+negative
        if n:
            from math import comb
            exact_p=min(1,2*sum(comb(n,k) for k in range(min(positive,negative)+1))/2**n)
        else:exact_p=None
        if n>1:
            rng=np.random.default_rng(20261115+i)
            draws=rng.choice(a,size=(20000,len(a)),replace=True)
            boot=np.quantile(np.median(draws,axis=1),[.025,.975]).tolist()
        else:boot=None
        out[stage]={"pairs":len(qs),"papers":len(study),"study_median_deltas":study,
                    "study_equal_median":float(np.median(a)) if len(a) else None,
                    "study_equal_mean":float(np.mean(a)) if len(a) else None,
                    "papers_positive":positive,"papers_negative":negative,
                    "descriptive_paper_bootstrap95":boot,
                    "posthoc_two_sided_sign_probability":exact_p,
                    "no_original_within_paper_sampling_covariance":True}
    return out

def compute(raw):
    p,source_audit=get_panels(raw)
    true,meta=exact_pairs(p)
    loose,loose_meta=exact_pairs(p,COARSE)
    return {
        "schema":"alan_B_global340_original_exact_fruit_seed_stage_v1",
        "author_source_sha256":SOURCE_SHA,"author_commit":COMMIT,
        "original_rows":len(raw),"reproductive_panels":len(p),
        "papers":p.study_ID.nunique(),"exact_keys":EXACT,
        "exact_pairs":len(true),"coarse_pairs_audit_only":len(loose),
        "paired_endpoint_summaries":summarize(true),
        "exact_source_quality":meta,
        "coarse_source_quality_audit_only":loose_meta,
        "source_exclusion_audit":source_audit,"source_pairs":true,
        "limits":[
          "Endpoint pairing by metadata does not prove identical initial flower or individual plants",
          "Seed/formed fruit is conditioned on fruit set; no reconstructed total seed per original flower",
          "Some comparisons can be observational or dark/illumination pooled",
          "No valid plant/experiment-level covariance or meta-analysis SE",
          "Study equal bootstrap and sign probability are posthoc descriptive, not confirmatory",
          "No light manipulation in this pooled diel access source; ALAN beta_V not fitted"
        ],
        "ALAN_reproduction_moderator_beta_V":None
    }

def main():
    out=compute(source())
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("ALAN_GLOBAL_AUTHOR_SOURCE_RESULT="+json.dumps({
       k:v for k,v in out.items() if k not in ("source_pairs","source_exclusion_audit",
           "exact_source_quality","coarse_source_quality_audit_only")
    },ensure_ascii=False,allow_nan=False))
if __name__=="__main__":main()
