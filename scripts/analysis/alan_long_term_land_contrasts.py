"""Audit the global screen and compare species-specific multi-decade slopes by landmass."""
import json
import sys
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

from alan_global_gbif_screen import model


def run(cache, results):
    out = results / "long_term_land_contrasts"
    out.mkdir(exist_ok=True)
    db = duckdb.connect(str(cache / "annual.duckdb"), read_only=True)
    data = db.sql("from trajectories").df()
    taxa = db.sql("from species_names").df()
    regional = pd.read_csv(results / "regional_light_traits/all_regions.csv.gz", dtype={"sp":str})
    supported = regional.query("blocks >= 3")
    shared = set.intersection(*(set(g.sp) for _,g in supported.groupby("land_class")))
    prior = pd.read_csv(results / "global_gbif_annual/long_term_screen.csv.gz",dtype={"sp":str})
    rows = []
    for sp,g in data.groupby("sp"):
        row = {"sp":sp}
        fits = {}
        for land,part in g.groupby("land_class"):
            row[land+"_coordinates"] = len(part)
            row[land+"_raw_blocks"] = part.block.nunique()
            fit = model(part)
            if fit is None:
                continue
            b,se,k,scores = fit
            row.update({land:b,land+"_se":se,land+"_blocks":k})
            fits[land] = (b,k,pd.Series(scores))
        if len(fits)==2:
            bi,ki,si = fits["island"]
            bm,km,sm = fits["mainland"]
            delta = bi-bm
            se = np.linalg.norm(si.subtract(sm,fill_value=0))
            df = min(ki,km)-1
            half = t.ppf(.975,df)*se
            row.update(delta=delta,low=delta-half,high=delta+half,p=2*t.sf(abs(delta)/se,df))
        rows.append(row)
    d = taxa.merge(prior[["sp","slope","q","brightening_mean_change_pp_per_decade"]],validate="1:1").merge(pd.DataFrame(rows),how="left",validate="1:1")
    d["shared_regional_1370"] = d.sp.isin(shared)
    pairs = d.p.notna()
    d["difference_q"] = np.nan
    if pairs.any():
        d.loc[pairs,"difference_q"] = multipletests(d.loc[pairs,"p"],method="fdr_bh")[1]
    d.to_csv(out/"all_species.csv.gz",index=False)
    both = d[pairs].sort_values("species")
    both.to_csv(out/"paired.csv",index=False)
    traits = pd.read_csv(results/"regional_light_traits/species_traits.csv",dtype={"sp":str})
    evidence = both.merge(traits[["sp","accepted_name","transport_review_status","night_exclusion_loss"]],on="sp",how="left",validate="1:1")
    updates = pd.read_csv(results/"mechanism_exemplars/transport_updates.csv")
    evidence.loc[evidence.accepted_name.isin(updates.plant),"transport_review_status"]="documented"
    evidence.to_csv(out/"paired_moth_evidence.csv",index=False)
    fig,ax=plt.subplots(figsize=(5,5))
    ax.scatter(both.mainland,both.island,c=np.where(both.slope.lt(0),"#19696a","#9f9a96"),alpha=.7)
    bounds = np.array([both[["mainland","island"]].min().min(),both[["mainland","island"]].max().max()])
    ax.plot(bounds,bounds,color="0.5",linestyle="--"); ax.axvline(0,color="0.8"); ax.axhline(0,color="0.8")
    ax.set(xlim=bounds,ylim=bounds,xlabel="Mainland long-term response",ylabel="Island long-term response",aspect="equal")
    fig.tight_layout(); fig.savefig(out/"paired.png",dpi=160); plt.close(fig)
    negative = both[both.slope.lt(0)]
    mainneg = both[both.mainland.lt(0)]
    audit=dict(input_species=len(taxa),regional_any=regional.sp.nunique(),regional_supported=supported.sp.nunique(),
               regional_shared=len(shared),trajectory_species=data.sp.nunique(),pooled_fitted=int(d.slope.count()),
               pooled_negative=int(d.slope.lt(0).sum()),pooled_negative_q05=int((d.slope.lt(0)&d.q.lt(.05)).sum()),
               island_fitted=int(d.island.count()),mainland_fitted=int(d.mainland.count()),paired=len(both),
               paired_from_1370=int(both.shared_regional_1370.sum()),difference_q05=int(both.difference_q.lt(.05).sum()),
               pooled_negative_paired=len(negative),pooled_negative_island_more_negative=int(negative.delta.lt(0).sum()),
               both_negative_island_more_negative=int(((negative.mainland<0)&(negative.island<negative.mainland)).sum()),
               pooled_negative_mainland_only=int(((negative.mainland<0)&(negative.island>=0)).sum()),
               pooled_negative_island_only=int(((negative.island<0)&(negative.mainland>=0)).sum()),
               mainland_negative_paired=len(mainneg),mainland_negative_amplified=int(mainneg.delta.lt(0).sum()),
               mainland_negative_weaker=int(((mainneg.delta>=0)&(mainneg.island<0)).sum()),
               mainland_negative_reversed=int(mainneg.island.ge(0).sum()))
    assert audit["regional_shared"]==1370 and len(d)==len(taxa) and d.sp.is_unique
    (out/"counts.json").write_text(json.dumps(audit,indent=2),encoding="utf8")
    print(audit)


if __name__=="__main__":
    run(*map(Path,sys.argv[1:3]))
