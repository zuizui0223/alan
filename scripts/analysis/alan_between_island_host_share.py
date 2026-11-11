"""Observed host shares and between-island comparisons within frozen region-years."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from alan_gecko_primary_fiveyear import primary


def summarize(s):
    s = s.copy()
    s["pct"] = 100 * s.all_moth / s.n_species
    raw = s.groupby("land_class").pct.agg(["size", "mean", "median"])
    islands = s.loc[s.land_class.eq("island")].copy()
    islands["land_id"] = islands.target_id.str.rsplit("|", n=1).str[-1]
    assert islands.land_id.str.fullmatch(r"[0-9]+").all()
    predictors = ["loglight", "logpop", "logphotos", "logspecies", "bio5", "lograin", "elevation_km"]
    strata = islands.groupby(["region5", "region_year", "land_id"])[["pct", *predictors]].mean().reset_index()
    pairs = strata[strata.groupby("region_year").land_id.transform("size").ge(2)]
    model = smf.ols("pct ~ C(region_year) + " + " + ".join(predictors), data=pairs).fit(
        cov_type="cluster", cov_kwds={"groups": pairs.region5}, use_t=True)
    assert np.linalg.matrix_rank(model.model.exog) == model.model.exog.shape[1]
    return dict(observed_site_year_shares=raw.to_dict(orient="index"), total_islands=islands.land_id.nunique(),
                matched_islands=pairs.land_id.nunique(), matched_regions=pairs.region5.nunique(),
                matched_region_years=pairs.region_year.nunique(), island_region_year_rows=len(pairs),
                between_island_pp_per_log1p=model.params.loglight,
                ci95_cluster_region=model.conf_int().loc["loglight"].tolist())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["inputs", "signed", "output"]:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    frame, _, _, _ = primary(args.inputs, args.signed)
    args.output.write_text(json.dumps(summarize(frame), indent=2) + "\n", encoding="utf8")
