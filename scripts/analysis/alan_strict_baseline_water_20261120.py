#!/usr/bin/env python3
"""Pre-outcome stricter original-water baseline matching on the frozen ALAN panel.

This is a falsification of the source-overlap gate, not an ALAN causal estimate.
The previous 187 pairs have baseline water SMD .72/.74 (SOURCE_SUPPORT_STOP).
Prespecified 0.5 and 0.25 IQR area0/freq0 calipers are compared without
looking at outcomes or claiming reduced-sample sensitivity as independent proof.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from alan_light_development_overlap_20261119 import (
    REPO, SOURCE_COMMIT, SOURCE_PATH, BASE, WATER_EXTRA,
    robust_scale, maxgap, read_source, pair_sites, summarize_pairs,
    same_cohort_water_adjustment
)

SEED=20261120
BASELINE_WATER=("area0","freq0")
CORE_BALANCE=("built","built0","light0","area0","freq0",
              "temp_trend","rain_trend")
PRESPECIFIED_FRACS=(.50,.25)

def water_caliper(data, fraction):
    """Nondegenerate scale if a water variable is zero-inflated (IQR=0).

    Fraction must scale the fallback itself, otherwise 0.5 and 0.25
    IQR produce IDENTICAL pair samples when the original IQR is zero.
    """
    finite=pd.to_numeric(data,errors="coerce").replace([np.inf,-np.inf],np.nan).dropna()
    q=float(finite.quantile(.75)-finite.quantile(.25))
    sd=float(finite.std()) if len(finite)>1 else 0.
    return max(fraction*max(q,.10*sd),1e-5)


def stricter_pairs(df, water_iqr_fraction):
    """Exact same original Hungarian algorithm, plus PRE-match water calipers."""
    if water_iqr_fraction not in PRESPECIFIED_FRACS:
        raise ValueError("NOT_PRESPECIFIED_CALIPER")
    cols=list(BASE)
    f=df.dropna(subset=cols+["unit_id","block","light"]).copy()
    f=f[np.isfinite(f[cols+["light"]]).all(axis=1)]
    f=f.sort_values("unit_id").reset_index(drop=True)
    scales={k:robust_scale(f[k]) for k in cols}
    limits={k:maxgap(f[k],.5 if k=="built" else 1.) for k in cols}
    limits.update({k:water_caliper(f[k],water_iqr_fraction) for k in BASELINE_WATER})
    limits["start"]=5.
    limits["span"]=10.
    f["treated"]=f.light.gt(0)
    pairs=[]
    for block,z in f.groupby("block",sort=True):
        high=z[z.treated].sort_values("unit_id")
        low=z[~z.treated].sort_values("unit_id")
        if not len(high) or not len(low):
            continue
        d2=np.zeros((len(high),len(low)),dtype=float)
        accepted=np.ones_like(d2,dtype=bool)
        for col in cols:
            diff=high[col].to_numpy(float)[:,None]-low[col].to_numpy(float)[None,:]
            d2+=(diff/scales[col])**2
            if col in ("built","built0","light0","start","span",
                       "area0","freq0"):
                accepted &= np.abs(diff)<=limits[col]
        light=high.light.to_numpy(float)[:,None]-low.light.to_numpy(float)[None,:]
        accepted &= light>=1.
        penalty=1.e6
        costs=np.where(accepted,d2,penalty*2)
        costs=np.column_stack([costs,np.full((len(high),len(high)),penalty)])
        rr,cc=linear_sum_assignment(costs)
        for ii,jj in zip(rr,cc):
            if jj>=len(low) or not accepted[ii,jj]:
                continue
            a,b=high.iloc[ii],low.iloc[jj]
            item={"block":block,
                  "bright_site_id":int(a.unit_id),
                  "comparison_site_id":int(b.unit_id),
                  "light_difference":float(a.light-b.light)}
            for k in cols+WATER_EXTRA:
                item["d_"+k]=float(a[k]-b[k])
            # Outcome only accessed AFTER valid matches are frozen.
            for k,col in (("outcome_difference","rate_difference"),
                ("candidate_difference","candidate_rate"),
                ("reference_difference","reference_rate")):
                item[k]=float(a[col]-b[col])
            pairs.append(item)
    return pd.DataFrame(pairs),f,scales,limits

def audit(df,frac):
    p,eligible,scales,limits=stricter_pairs(df,frac)
    if p.empty:
        return {"status":"NO_COMMON_SUPPORT","pairs":0,"blocks":0,
                "caliper_fraction_IQR":frac}
    match_summary=summarize_pairs(p,eligible,scales)
    sensitivities=same_cohort_water_adjustment(p)
    critical={k:match_summary["all_matched_covariate_SMD"].get(k)
              for k in CORE_BALANCE}
    safe=(match_summary["status"]=="EXPLORATORY_COMMON_SUPPORT")
    unique_ids=pd.concat([p.bright_site_id,p.comparison_site_id])
    if unique_ids.duplicated().any():
        raise RuntimeError("HUNGARIAN_REPLACEMENT_VIOLATION")
    pairs_digest=hashlib.sha256(
        "\n".join(f"{int(z.bright_site_id)},{int(z.comparison_site_id)}"
            for z in p.sort_values(["block","bright_site_id"]).itertuples()
        ).encode()).hexdigest()
    return {
        "caliper_fraction_IQR":frac,
        "caliper_water_area0":float(limits["area0"]),
        "caliper_water_freq0":float(limits["freq0"]),
        "matched_id_set_sha256":pairs_digest,
        "match":match_summary,
        "critical_baseline_SMD":critical,
        "all_critical_SMD_within_0_25":all(
            np.isfinite(v) and abs(v)<=.25 for v in critical.values()),
        "same_subset_nested_water_change_sensitivity":sensitivities,
        "scientific_claim_gate":"EXPLORATORY_SOURCE_SUPPORTED_ONLY" if safe
                else "SOURCE_SUPPORT_STOP",
        "subset_not_same_187_pairs":True,
        "source_postselected_58_vs_305_plants":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="")
    args=p.parse_args()
    df,source_sha=read_source()
    orig,_,_=pair_sites(df,False)
    if len(orig)!=187 or orig.block.nunique()!=30:
        raise RuntimeError("FROZEN_ORIGINAL_187_PAIR_MISMATCH")
    original_digest=hashlib.sha256(
        "\n".join(f"{int(z.bright_site_id)},{int(z.comparison_site_id)}"
          for z in orig.sort_values(["block","bright_site_id"]).itertuples()
        ).encode()).hexdigest()
    results=[audit(df,frac) for frac in PRESPECIFIED_FRACS]
    report={
        "dataset":f"{SOURCE_COMMIT}:{SOURCE_PATH}",
        "source_sha256":source_sha,
        "original_reference_pairs":187,
        "original_reference_blocks":30,
        "original_reference_ids_sha256":original_digest,
        "original_reference_area0_SMD":.7224709411617956,
        "source_baseline_water_IQR":{k:float(df[k].quantile(.75)-df[k].quantile(.25)) for k in BASELINE_WATER},
        "source_baseline_water_sd":{k:float(df[k].std()) for k in BASELINE_WATER},
        "original_reference_freq0_SMD":.7445843410774747,
        "source_is_postselected_exploratory_not_causal":True,
        "matching_outcome_blind":True,
        "prespecified_order":"primary 0.5 IQR baseline water, sensitivity 0.25 IQR",
        "water_sensitivity_not_causal_mediation":True,
        "comparisons":results
    }
    text=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)
    print("ALAN_STRICT_BASELINE_WATER_MATCH_RESULTS\n"+text,flush=True)
    if args.out:
        o=Path(args.out);o.parent.mkdir(parents=True,exist_ok=True)
        o.write_text(text+"\n")
if __name__=="__main__":
    main()
