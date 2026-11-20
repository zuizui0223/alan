#!/usr/bin/env python3
"""Explain sensitivity of verified ALAN wet-plant light contrasts to spatial blocks.

The frozen original 1999-2023 matched wet-plant panel contains preselected 58
candidate and 305 reference wetland-affinity taxa. These are botanical survey
recording-rate associations, not measured population changes or causal ALAN
treatment effects. Do not reverse the outcome-selection gate.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from alan_light_development_overlap_20261119 import (
    read_source, pair_sites, WATER_EXTRA
)
from alan_strict_baseline_water_20261120 import stricter_pairs

SOURCE_SHA256="4d151698e0da5ee754326220d3bb817b141da00770299e7f8910366d8771e74f"
TERMS=["light_difference","d_built","d_built0","d_light0",
       "d_temp_trend","d_rain_trend","d_freq0","d_area0"]+["d_freq_trend","d_area_trend"]
OUTCOMES=("outcome_difference","candidate_difference","reference_difference")
SEED=20261120

def fit_slope(data, target, terms=TERMS):
    """Original no-intercept least squares, difference in pp/10DN/decade."""
    if not len(data) or target not in data:
        return float("nan")
    x=data[terms].to_numpy(float)
    y=data[target].to_numpy(float)
    if (not np.isfinite(x).all() or not np.isfinite(y).all() or
        len(data)<=len(terms)+5 or np.linalg.matrix_rank(x)<len(terms)):
        return float("nan")
    return 10.*float(np.linalg.lstsq(x,y,rcond=None)[0][0])

def light_residual_support(data, terms=TERMS):
    """Diagnose identification of light after nuisance residualization."""
    x=data[terms].to_numpy(float)
    if not np.isfinite(x).all():
        raise RuntimeError("NONFINITE_PREDICTORS")
    light=x[:,0]
    other=x[:,1:]
    if np.linalg.matrix_rank(other)<other.shape[1]:
        return dict(status="CONTROLS_COLLINEAR_STOP")
    resid=light-other@np.linalg.lstsq(other,light,rcond=None)[0]
    raw=float(light@light)
    remaining=float(resid@resid)
    return dict(status="DESCRIPTIVE_LINEAR_SUPPORT",
      light_after_controls_fraction_of_sumsquares=remaining/raw if raw>0 else None,
      light_variance_inflation_factor=raw/remaining if remaining>1e-10 else None,
      light_residual_rms_DN_per_decade=float(np.sqrt(np.mean(resid**2))),
      unadjusted_rms_DN_per_decade=float(np.sqrt(np.mean(light**2))))

def block_bootstrap_components(pairs,n_draws=1500,seed=SEED):
    """All outcome components re-use the same block draws, preserving covariance."""
    groups=[g for _,g in pairs.groupby("block",sort=True)]
    if len(groups)<10:
        return dict(status="BLOCK_SUPPORT_STOP",blocks=len(groups))
    rng=np.random.default_rng(seed)
    point={k:fit_slope(pairs,k) for k in OUTCOMES}
    if not all(np.isfinite(v) for v in point.values()):
        return dict(status="DESIGN_RANK_STOP")
    draw=[]
    for _ in range(n_draws):
        p=pd.concat([groups[j] for j in rng.integers(0,len(groups),size=len(groups))],
                    ignore_index=True)
        draw.append([fit_slope(p,k) for k in OUTCOMES])
    a=np.asarray(draw,dtype=float)
    good=a[np.isfinite(a).all(axis=1)]
    if len(good)<500:
        return dict(status="BOOTSTRAP_RANK_STOP",valid_draws=len(good))
    ci=np.quantile(good,[.025,.975],axis=0)
    # Algebraic identity must hold in point estimates and in each resample.
    if not np.allclose(good[:,1]-good[:,2],good[:,0],atol=1e-8) or not np.isclose(
        point["candidate_difference"]-point["reference_difference"],point["outcome_difference"],atol=1e-8):
        raise RuntimeError("CANDIDATE_REFERENCE_DECOMPOSITION_STOP")
    return dict(status="DESCRIPTIVE_COMPONENTS",
      source_matched_pairs=int(len(pairs)),spatial_blocks=int(len(groups)),
      coefficient_point_pp_per_10DN=point,
      block_bootstrap_95pct={k:ci[:,i].tolist() for i,k in enumerate(OUTCOMES)},
      bootstrap_complete_rank_draws=int(len(good)),
      outcome_identity_verified=True)

def leave_one_spatial_block_out(pairs):
    """Same pairs except the omitted 100km source block. Not a cross-validation."""
    all_slope=fit_slope(pairs,"outcome_difference")
    if not np.isfinite(all_slope):
        return dict(status="DESIGN_RANK_STOP")
    rows=[]
    for bid in sorted(pairs.block.unique()):
        p=pairs[pairs.block.ne(bid)]
        b=fit_slope(p,"outcome_difference")
        rows.append({"omitted_block":int(bid),
            "omitted_pairs":int((pairs.block==bid).sum()),
            "remaining_blocks":int(p.block.nunique()),
            "slope_without_block":float(b) if np.isfinite(b) else None,
            "change_from_full":float(b-all_slope) if np.isfinite(b) else None})
    valid=[z["slope_without_block"] for z in rows if z["slope_without_block"] is not None]
    if not valid:
        return dict(status="ALL_BLOCK_OMISSIONS_RANK_STOP")
    influential=max((z for z in rows if z["change_from_full"] is not None),
                     key=lambda z:abs(z["change_from_full"]))
    return dict(status="DESCRIPTIVE_LEAVE_ONE_BLOCK_OUT",
      n_total_blocks=len(rows),n_fitted=len(valid),
      full_sample_slope_pp_per_10DN=float(all_slope),
      n_negative=int(sum(v<0 for v in valid)),
      negative_fraction=float(np.mean(np.asarray(valid)<0)),
      min_slope=float(min(valid)),max_slope=float(max(valid)),
      median_slope=float(np.median(valid)),
      most_influential_block=influential,
      every_block_omission=rows)

def pair_keys(p):
    return set(zip(p.bright_site_id.astype(int),p.comparison_site_id.astype(int)))

def cohort_relation(original, current):
    a,b=pair_keys(original),pair_keys(current)
    olds=set(original.bright_site_id.astype(int))|set(original.comparison_site_id.astype(int))
    news=set(current.bright_site_id.astype(int))|set(current.comparison_site_id.astype(int))
    return dict(original_pairs=int(len(original)),current_pairs=int(len(current)),
      identical_original_pair_id_combinations=int(len(a&b)),
      original_pair_combinations_lost=int(len(a-b)),
      new_pair_combinations_not_in_original=int(len(b-a)),
      current_sites_shared_with_original=int(len(news&olds)),
      current_sites_not_in_original=int(len(news-olds)),
      unique_original_sites=int(len(olds)))

def analyze(data,label,frac,original):
    if label=="original":
        p=data
    else:
        p,_,_,_=stricter_pairs(data,frac)
    if p.empty: return dict(status="NO_MATCHED_PAIRS",label=label)
    if p.block.nunique()<10: return dict(status="TOO_FEW_BLOCKS",label=label)
    if len(pd.concat([p.bright_site_id,p.comparison_site_id]).unique())!=2*len(p):
        raise RuntimeError("SOURCE_MATCHED_SITE_REUSE")
    return dict(label=label,caliper_fraction=frac,
        pair_relation=cohort_relation(original,p),
        change_covariate_pair_mean={"d_freq_trend":float(p.d_freq_trend.mean()),
                                    "d_area_trend":float(p.d_area_trend.mean())},
        component_estimates=block_bootstrap_components(p),
        block_influence=leave_one_spatial_block_out(p),
        light_support=light_residual_support(p))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",default="")
    a=parser.parse_args()
    data,source_sha=read_source()
    if source_sha!=SOURCE_SHA256:
        raise RuntimeError("ORIGINAL_PINNED_SOURCE_HASH_CHANGED")
    orig,_,_=pair_sites(data,False)
    if len(orig)!=187 or orig.block.nunique()!=30:
        raise RuntimeError("FROZEN_ORIGINAL_PAIR_COUNT_CHANGED")
    results=[
      analyze(orig,"original",None,orig),
      analyze(data,"baseline_water_0p5",.5,orig),
      analyze(data,"baseline_water_0p25",.25,orig)]
    for result,expected,blocks in zip(results,[187,133,128],[30,26,25]):
        c=result["component_estimates"]
        if c.get("source_matched_pairs")!=expected or c.get("spatial_blocks")!=blocks:
            raise RuntimeError("FROZEN_MATCHED_SUBSAMPLE_CHANGED")
    report={"study":"ALAN_wet_habitat",
        "original_source_commit_and_sha256":f"6a9a484:{SOURCE_SHA256}",
        "outcomes":"paired differences in botanical survey recording-rate trend (pp/decade)",
        "slope_unit":"pp of recording-rate trend per 10 DN/decade brightening contrast",
        "source_selected_plants":True,"causal_ALAN_effect":False,
        "independent_botanical_validation":False,
        "cross_cohort_coefficient_difference_is_not_mediation":True,
        "critical_note":"No original 374 paired site coordinates or wetland vegetation covariate; this is spatial stability and outcome component diagnostic only",
        "cohorts":results}
    s=json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)
    print("ALAN_WET_BLOCK_COMPONENT_AUDIT\n"+s,flush=True)
    if a.out:
        p=Path(a.out);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(s+"\n")

if __name__=="__main__":
    main()
