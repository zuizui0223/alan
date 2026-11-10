#!/usr/bin/env python3
"""Empirical support, not biology: can dark fragmentation vary at similar brightness?

Consumes only the SHA-signed actual PCNL2024 5-km land-masked exposure audit
obtained from ORIGINAL 2025 photographed island site-years. Never creates
plant/insect abundance, species non-detection or a causal ALAN model.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA="32327abb37ee85f3d1d9e7446d9d742fb5adf1457dab60fae1aae791c29f7837"
SOURCE_SCHEMA="alan_dark_refugia_original_2025_sites_prior_2024_PCNL_GSHHG_v1"
THRESHOLDS=(0.,1.,3.)
SCHEMA="ALAN_dark_refugia_within_island_conditional_contrast_feasibility_v1"


def require(ok, msg):
    if not bool(ok):
        raise ValueError("ALAN_DARK_REFUGIA_SUPPORT: "+str(msg))


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for c in iter(lambda:f.read(4<<20),b""):
            h.update(c)
    return h.hexdigest()


def readable_source(raw, receipt, check_signature=True):
    info=json.loads(Path(receipt).read_text())
    require(info["schema"]==SOURCE_SCHEMA,"Different original 2025 site/raster cohort")
    if check_signature:
        require(info["source_output_csv_sha256"]==SOURCE_SHA and
                sha256(raw)==SOURCE_SHA,"Original raster-pixel source SHA changed")
    d=pd.read_csv(raw,dtype={"island":str})
    req={"island","latbin","lonbin","photo_observation_year",
         "threshold_PCNL_native","status","dark_land_fraction",
         "largest_dark_component_fraction_of_dark",
         "largest_dark_component_fraction_of_land","local_mean_PCNL_native",
         "n_valid_original_PCNL_pixels","n_dark_raster_pixels"}
    require(req.issubset(d.columns),"Not genuine original raster metrics")
    require(not d.duplicated(["island","latbin","lonbin","photo_observation_year",
                              "threshold_PCNL_native"]).any(),"Pseudo-replicated source")
    require(d.photo_observation_year.eq(2025).all(),"Wrong original observation year")
    require(set(d.threshold_PCNL_native)==set(THRESHOLDS),"Thresholds not prespecified")
    require(d.groupby(["island","latbin","lonbin"],observed=True).size().eq(3).all(),
            "Each original source site needs exactly three thresholds")
    good=d.loc[d.status.eq("observed_satellite_geometry_only")].copy()
    require(good.dark_land_fraction.between(0,1).all() and
            good.local_mean_PCNL_native.ge(0).all(),"Impossible original light measure")
    return d,good,info


def summarize(real, good):
    results={}
    for t in THRESHOLDS:
        x=good.loc[good.threshold_PCNL_native.eq(t)].copy()
        q=lambda z: {str(v):float(z.quantile(v)) for v in (.05,.25,.50,.75,.95)} if len(z) else {}
        x["site"]=(x.island+"|"+x.latbin.astype(str)+"|"+x.lonbin.astype(str))
        dark=x.loc[x.n_dark_raster_pixels.gt(0)].copy()
        # Require each compared site to have recorded dark land, as conditional
        # connectivity is undefined when all original pixels are bright.
        contrast=[]
        for island,g in dark.groupby("island",sort=False):
            if len(g)<2: continue
            for a,b in itertools.combinations(g.itertuples(index=False),2):
                if a.site==b.site:continue
                ds=abs(float(a.dark_land_fraction)-float(b.dark_land_fraction))
                dc=abs(float(a.largest_dark_component_fraction_of_dark)-
                       float(b.largest_dark_component_fraction_of_dark))
                dl=abs(np.log1p(float(a.local_mean_PCNL_native))-
                       np.log1p(float(b.local_mean_PCNL_native)))
                contrast.append((island,ds,dc,dl))
        pairs=pd.DataFrame(contrast,columns=["island","dark_fraction_difference",
                                            "connectivity_difference","logmean_light_difference"])
        # Strict within-island exposure-common-support tests, independent of
        # all species/traits/outcomes, not threshold-tuned p values.
        def support(min_dc,max_df,max_dl=None):
            if pairs.empty:return {"n_islands":0,"n_pairs":0}
            ok=(pairs.connectivity_difference>=min_dc)&(pairs.dark_fraction_difference<=max_df)
            if max_dl is not None:
                ok &= pairs.logmean_light_difference<=max_dl
            e=pairs.loc[ok]
            return {"n_islands":int(e.island.nunique()),"n_pairs":len(e)}
        results[str(t)]={
            "n_original_valid_sites":int(len(x)),"n_independent_islands":int(x.island.nunique()),
            "dark_fraction_quantiles":q(x.dark_land_fraction),
            "conditional_connectivity_quantiles":q(dark.largest_dark_component_fraction_of_dark),
            "n_sites_with_no_detected_dark_land":int(x.n_dark_raster_pixels.eq(0).sum()),
            "n_islands_with_two_or_more_valid_photo_sites":
                int(x.groupby("island").size().ge(2).sum()),
            "n_within_island_dark_site_pairs":len(pairs),
            "potential_matched_pairs_connectivity_delta0p10_darkfraction_delta0p10":
                support(.10,.10),
            "stronger_contrasts_connectivity_delta0p20_darkfraction_delta0p10":
                support(.20,.10),
            "joint_similar_total_light_connectivity_delta0p10_darkfraction_delta0p10_loglight_delta0p2":
                support(.10,.10,.20),
            "threshold_PCNL_native":t
        }
    # Check thresholds 0 and 1 coarsely, not assume they are independent replicates.
    d0=good.loc[good.threshold_PCNL_native.eq(0)].set_index(["island","latbin","lonbin"])
    d1=good.loc[good.threshold_PCNL_native.eq(1)].set_index(["island","latbin","lonbin"])
    common=d0.index.intersection(d1.index)
    same=(d0.loc[common,"dark_land_fraction"].to_numpy() ==
          d1.loc[common,"dark_land_fraction"].to_numpy()) if len(common) else []
    return {"schema":SCHEMA,"source_original_site_raster_SHA256":SOURCE_SHA,
            "not_species_moth_or_pollination":True,
            "thresholds_exploratory_sensitivities_not_independent_samples":True,
            "same_original_sites_threshold0_and1":len(common),
            "fraction_sites_exactly_same_dark_land_PCNL0_and1":
                float(np.mean(same)) if len(common) else None,
            "threshold_results":results}


def run(raw, receipt, out):
    a,b,_=readable_source(raw,receipt)
    result=summarize(a,b)
    out=Path(out)
    require(not out.exists(),"Refuse replace previously reported support")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_DARK_REFUGIA_TRUE_CONDITIONAL_SUPPORT",json.dumps(result,sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--receipt",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.receipt,a.out)
