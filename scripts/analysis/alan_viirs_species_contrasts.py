"""Translate the frozen species slopes into contrasts on the observed VIIRS scale."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import primary

OUT = Path(__file__).resolve().parents[2] / "results/species_interpretation"


def run(inputs, signed):
    sites, _, members, _ = primary(inputs, signed)
    assert sites.key.is_unique and len(sites) == 10698
    sites = sites.assign(radiance=sites.viirs_median.clip(lower=0))
    profile = sites.groupby("land_class").radiance.quantile([.1,.25,.5,.75,.9]).unstack()
    profile.columns = ["p10","p25","p50","p75","p90"]
    profile["site_years"] = sites.land_class.value_counts()
    profile.to_csv(OUT / "viirs_observed_quantiles.csv")
    print(profile)
    s = pd.read_csv(OUT / "island_sign_reproductive_evidence.csv").set_index("plant")
    island = sites[sites.land_class.eq("island")]
    obs = members.merge(island[["key","radiance"]],on="key",validate="many_to_one")
    radiances = obs.groupby("wcvp_accepted_name").radiance
    s["island_recorded_min"] = radiances.min()
    s["island_recorded_max"] = radiances.max()
    anchor = profile.loc["island"]
    for start, end in [("p25","p50"),("p50","p75"),("p25","p75")]:
        delta = np.log1p(anchor[end]) - np.log1p(anchor[start])
        for side in ["island","mainland"]:
            s[side+"_"+start+"_"+end+"_pp"] = s[side+"_pp"] * delta
    s["recorded_span_covers_island_iqr"] = s.island_recorded_min.le(anchor.p25) & s.island_recorded_max.ge(anchor.p75)
    checks = radiances.size().reindex(s.index)
    np.testing.assert_array_equal(checks,s.island_site_years)
    s.to_csv(OUT / "viirs_species_contrasts.csv")
    neg = s[s.island_pp.lt(0)]
    print(neg.filter(regex="island_p[25]").median(), neg.recorded_span_covers_island_iqr.sum())
    print(neg.island_p25_p75_pp.quantile([.25,.5,.75]))
    for group, cohort in neg.groupby(neg.transport_review_status):
        print(group,len(cohort),cohort.island_p25_p75_pp.median())
    print(s.loc[["Lythrum salicaria","Silene latifolia","Asclepias syriaca"],["island_p25_p75_pp","mainland_p25_p75_pp","island_recorded_min","island_recorded_max"]].to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("signed", type=Path)
    args = parser.parse_args()
    run(args.inputs,args.signed)
