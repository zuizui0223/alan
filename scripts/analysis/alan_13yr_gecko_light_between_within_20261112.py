#!/usr/bin/env python3
"""Decompose 13-year actual GBIF gecko/moth PHOTO × prior PCNL into space vs time.

Run on EXACT same 454 repeatedly photographed GSHHG island 0.05-degree grids.
The source is observational and positive-record only. Do not interpret a
GloBI food-web mechanism or causation from this decomposition.
"""
import argparse,json,math
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t
from alan_13year_moth_gecko_repeat_site_FE_20261112 import read_pair,repeated_sites
from alan_13year_moth_gecko_same_island_source_ALAN_20261112 import cohort


def fit(d,response="gecko",weights="island_equal"):
    a=d.copy()
    a["meanlight"]=a.groupby("site").loglight.transform("mean")
    a["withinlight"]=a.loglight-a.meanlight
    a["mean_insect_effort"]=a.groupby("site").log_insect.transform("mean")
    a["within_insect_effort"]=a.log_insect-a.mean_insect_effort
    a["mean_reptile_effort"]=a.groupby("site").log_reptile.transform("mean")
    a["within_reptile_effort"]=a.log_reptile-a.mean_reptile_effort
    cols=["meanlight","withinlight","mean_insect_effort",
          "within_insect_effort","mean_reptile_effort",
          "within_reptile_effort","latitude"]
    if response=="gecko":
        y=100*a.gecko_share
    elif response=="gecko_smoothed":
        y=100*(a.gecko_source_photo_count+.5)/(
            a.original_total_reptile_photo_source_events+1)
    elif response=="moth":
        y=a.moth_share
    else:
        raise ValueError(response)
    n=len(a)
    if weights=="island_equal":
        w=1/a.groupby("island").island.transform("size").to_numpy(float)
    elif weights=="reptile_count":
        w=(1/a.groupby("island").island.transform("size").to_numpy(float) *
           np.minimum(a.original_total_reptile_photo_source_events.to_numpy(float),10))
    else:
        raise ValueError(weights)
    group,levels=np.unique(a.islandyear,return_inverse=True)
    totals=np.bincount(levels,weights=w,minlength=len(group))
    V=np.column_stack((a[cols].to_numpy(float),y.to_numpy(float)))
    for j in range(V.shape[1]):
        means=np.bincount(levels,weights=w*V[:,j],minlength=len(group))/totals
        V[:,j]=(V[:,j]-means[levels])*np.sqrt(w)
    X,Y=V[:,:-1],V[:,-1]
    B,_,rank,_=np.linalg.lstsq(X,Y,rcond=1e-11)
    if rank!=len(cols):
        raise ValueError("Source within-between FE design rank deficient")
    residual=Y-X@B
    cluster,idx=np.unique(a.island,return_inverse=True);G=len(cluster)
    scores=np.column_stack([
        np.bincount(idx,weights=X[:,j]*residual,minlength=G)
        for j in range(len(cols))
    ])
    bread=np.linalg.inv(X.T@X)
    cov=bread@scores.T@scores@bread*(G/(G-1))*((len(a)-1)/(len(a)-len(cols)))
    tcrit=float(student_t.ppf(.975,G-1))
    estimates={}
    for key,r in (
       ("between_site_meanlight",np.r_[1.,np.zeros(len(cols)-1)]),
       ("within_site_light_change",np.r_[0.,1.,np.zeros(len(cols)-2)]),
       ("within_minus_between",np.r_[-1.,1.,np.zeros(len(cols)-2)])
    ):
        beta=float(r@B)
        se=math.sqrt(float(r@cov@r))
        estimates[key]={
            "beta":beta,
            "ci95":[beta-tcrit*se,beta+tcrit*se],
            "p_nominal":float(2*student_t.sf(abs(beta/se),G-1))
        }
    return {
        "source_cell_years":n,"sites":int(a.site.nunique()),
        "original_GSHHG_islands":G,"islandyear_FE":len(group),
        "response":response,"source_weights":weights,"effects":estimates
    }


def run(insects,reptiles):
    a=read_pair(insects,reptiles)
    sp=cohort(a)
    rep=repeated_sites(sp,.02)
    if (len(rep),rep.site.nunique(),rep.island.nunique())!=(1979,454,30):
        raise ValueError("The 13-year prior-source cohort changed")
    result={
        f"{resp}__{weight}":fit(rep,resp,weight)
        for resp,weight in (
            ("gecko","island_equal"),("gecko","reptile_count"),
            ("gecko_smoothed","island_equal"),("moth","island_equal")
        )
    }
    return {
        "schema":"alan_GBIF_13yr_same_source_site_between_within_priorPCNL_v1",
        "original_source_insects_SHA256":
           "64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852",
        "original_source_reptiles_SHA256":
           "03aa855c9e7691855bfee4baa099f0d3a31db8d53066ec36a95afda48c3ba339",
        "years":list(range(2013,2026)),
        "observational_same_site_support":{"sites":454,"cell_years":1979,"islands":30},
        "design":"year×island FE, mean source-site light versus light deviation from same site's own multi-year mean; photo effort decomposed into site means and within-site changes; island-cluster CR1. No random temporal lighting assignment.",
        "models":result,
        "limitations":{
            "posthoc_but_source_identity_FROZEN":True,
            "both_components_not_experimental_ALAN_effect":True,
            "observational_photo_ratios_not_true_abundance":True,
            "GloBI_links_not_field_moth_predation":True,
            "an_unrecorded_animal_is_not_biologically_absent":True,
            "no_local_effective_pollen_transfer_observed":True,
            "sampled_years_not_equally_represented":True,
            "no_private_species_observer_or_geographic_site_rows_exported":True
        }
    }


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--insects",required=True,type=Path)
    p.add_argument("--reptiles",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    x=p.parse_args()
    if x.out.exists():
        raise ValueError("Refuse overwrite")
    o=run(x.insects,x.reptiles)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(o,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("ALAN_GBIF_WITHIN_BETWEEN_SOURCE",json.dumps(o["models"]["gecko__island_equal"],sort_keys=True))
