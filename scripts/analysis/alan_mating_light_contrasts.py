"""Equal-species contrasts in contributions to the frozen island-minus-mainland slope."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import primary
from alan_globi_crossland_light_robust_20261111 import fit

ROOT = Path(__file__).resolve().parents[2]
GAP = "beta_island_minus_mainland_slope_pp_per_log1p"


def analyze(s, hosts, out):
    traits = pd.read_csv(ROOT/"results/moth_pollen_evidence/shared_544_evidence.csv").set_index("plant")
    assert traits.index.is_unique and len(traits)==544
    hosts = hosts.loc[hosts.wcvp_accepted_name.isin(traits.index), ["key","wcvp_accepted_name"]].drop_duplicates()
    assert set(hosts.wcvp_accepted_name)==set(traits.index)
    frame = s.copy()
    rows = []

    def estimate(label, weights):
        # Linearity of the fixed design lets the same fit estimate a contrast
        # of species coefficients, with their shared spatial covariance retained.
        values = hosts.wcvp_accepted_name.map(weights).groupby(hosts.key).sum()
        frame["contrast_moth"] = frame.key.map(values).fillna(0)
        r = fit(frame,"contrast","moth",balance=True)
        assert r["estimable"] and r["source_5degree_clusters"]==78 and r["n_site_years"]==10698
        lo,hi = r["95CI_cluster_5deg"]
        rows.append(dict(estimate=label, species=len(weights),island=r["beta_island_pp_per_log1p"],
                         mainland=r["beta_mainland_pp_per_log1p"], gap=r[GAP],ci_low=lo,ci_high=hi,p=r["p_nominal"]))
        return r[GAP]

    pl = traits.glopl_median_log_response_ratio.dropna()
    assert len(pl)==64
    centered = pl-pl.mean()
    slope_weights = centered/centered.pow(2).sum()
    constant = pd.Series(1/len(pl),index=pl.index)
    slope = estimate("PL_per_unit_log_response_ratio",slope_weights)
    q = pl.quantile([.25,.75])
    low,high = [estimate("PL_prediction_"+label,constant+(v-pl.mean())*slope_weights)
                for label,v in zip(["Q25","Q75"],q)]
    change = estimate("PL_Q75_minus_Q25",slope_weights*(q.iloc[1]-q.iloc[0]))
    np.testing.assert_allclose(change,[high-low],atol=1e-12)
    known = traits.loc[traits.sc_source_class.ne("unknown"),"sc_source_class"]
    assert len(known)==246
    means = {}
    classes = {}
    for category,species in known.groupby(known):
        weights = pd.Series(1/len(species),index=species.index)
        classes[category] = weights
        means[category] = estimate(category,weights)
    contrast = estimate("SI_minus_SC",pd.concat([classes["SI"],-classes["SC"]]))
    np.testing.assert_allclose(contrast,means["SI"]-means["SC"],atol=1e-12)
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/"contrasts.csv",index=False)
    matched = hosts.merge(s[["key","region5","land_class"]],on="key",validate="many_to_one")
    support = {}
    for label, names in {"PL":pl.index,"all_classified":known.index,"SC":classes["SC"].index,"SI":classes["SI"].index}.items():
        data = matched[matched.wcvp_accepted_name.isin(names)]
        support[label] = dict(positive_site_years=data.groupby("land_class").key.nunique().to_dict(),
                              positive_regions=data.region5.nunique(),
                              positive_regions_both_sides=int(data.groupby("region5").land_class.nunique().eq(2).sum()))
    (out/"support.json").write_text(json.dumps(dict(pl_q25=q.iloc[0],pl_q75=q.iloc[1],support=support),indent=2)+"\n")


if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs","signed","out"):
        ap.add_argument(arg,type=Path)
    args = ap.parse_args()
    sites, _, plants, _ = primary(args.inputs,args.signed)
    analyze(sites,plants.loc[plants.all_moth],args.out)
