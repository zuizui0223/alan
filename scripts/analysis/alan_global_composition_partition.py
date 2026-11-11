"""Additively partition the frozen community response by cross-land plant membership."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import ID, primary
from alan_globi_crossland_light_robust_20261111 import fit


def run(inputs, signed, output):
    sites, _, plants, _ = primary(inputs, signed)
    hosts = plants.loc[plants.all_moth].merge(sites[["key", "land_class"]], on="key", validate="many_to_one")
    occurrence = hosts.groupby("wcvp_accepted_name").land_class.nunique()
    shared = set(occurrence[occurrence.eq(2)].index)
    previous = pd.read_csv(output.parent / "shared_544.csv")
    assert shared == set(previous.plant) and len(shared) == 544
    counts = hosts[hosts.wcvp_accepted_name.isin(shared)].groupby("key")[ID].nunique()
    sites["shared_moth"] = sites.key.map(counts).fillna(0)
    sites["remaining_moth"] = sites.all_moth - sites.shared_moth
    assert sites.remaining_moth.ge(0).all()
    result = {scope: fit(sites, scope, "moth", balance=True) for scope in ("all", "shared", "remaining")}
    for term in ("beta_island_minus_mainland_slope_pp_per_log1p", "beta_mainland_pp_per_log1p", "beta_island_pp_per_log1p"):
        np.testing.assert_allclose(result["all"][term], result["shared"][term] + result["remaining"][term], atol=1e-10, rtol=0)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs", "signed", "output"):
        parser.add_argument(arg, type=Path)
    run(**vars(parser.parse_args()))
