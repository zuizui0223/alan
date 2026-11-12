"""Localize the largest signed species contributions in the unchanged five-year fit."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from alan_shared_species_contributions import RESULTS, influence, primary


def run(inputs, signed, out):
    ledger = pd.read_csv(RESULTS / "species_interpretation/shared_544.csv")
    assert len(ledger) == ledger.plant.nunique() == 544
    ranked = ledger.sort_values(["gap_pp", "plant"])
    focal = pd.concat([ranked.head(3), ranked.tail(3).iloc[::-1]]).set_index("plant")
    focal["rank_in_tail"] = [1, 2, 3, 1, 2, 3]
    sites, _, members, _ = primary(inputs, signed)
    design = sites[["key", "land_class", "year", "region5", "region_year", "viirs_median"]].copy()
    design["contribution_pp"] = 100 * influence(sites)[1] / sites.n_species
    rows = members.loc[members.wcvp_accepted_name.isin(focal.index),
                       ["key", "wcvp_accepted_name"]].rename(columns={"wcvp_accepted_name": "plant"})
    assert not rows.duplicated().any()
    rows = rows.merge(design, on="key", validate="m:1")
    rows["paired"] = rows.groupby(["plant", "region_year"]).land_class.transform("nunique").eq(2)
    by_year = rows.groupby(["plant", "year"]).agg(
        gap_contribution_pp=("contribution_pp", "sum"),
        island_site_years=("land_class", lambda x: (x == "island").sum()),
        mainland_site_years=("land_class", lambda x: (x == "mainland").sum()))
    levels = pd.MultiIndex.from_product([focal.index, sorted(sites.year.unique())], names=["plant", "year"])
    by_year = by_year.reindex(levels, fill_value=0)
    annual = by_year.gap_contribution_pp.unstack()
    np.testing.assert_allclose(annual.sum(axis=1).loc[focal.index], focal.gap_pp, atol=1e-10)
    by_region = rows.groupby(["plant", "region5"]).contribution_pp.sum()
    largest = by_region.loc[by_region.abs().groupby(level=0).idxmax()].reset_index().set_index("plant")
    focal["largest_region5"] = largest.region5
    focal["largest_region_pp"] = largest.contribution_pp
    focal["regions_recorded"] = by_region.groupby(level=0).size()
    focal["largest_region_share_abs"] = by_region.abs().groupby(level=0).max() / by_region.abs().groupby(level=0).sum()
    focal["years_contributing_same_sign"] = np.sign(annual).eq(np.sign(focal.gap_pp), axis=0).sum(axis=1)
    focal["paired_region_year_contribution_pp"] = rows[rows.paired].groupby("plant").contribution_pp.sum()
    focal["paired_region_year_contribution_pp"] = focal.paired_region_year_contribution_pp.fillna(0)
    focal["remaining_contribution_pp"] = focal.gap_pp - focal.paired_region_year_contribution_pp
    columns = ["rank_in_tail", "source_class", "island_pp", "mainland_pp", "gap_pp", "island_site_years",
               "mainland_site_years", "region_years_with_species_on_both_sides", "regions_recorded",
               "years_contributing_same_sign", "largest_region5", "largest_region_pp", "largest_region_share_abs",
               "paired_region_year_contribution_pp", "remaining_contribution_pp"]
    out.mkdir(parents=True, exist_ok=True)
    focal[columns].to_csv(out / "ranked_contributions.csv")
    by_year.to_csv(out / "annual_contributions.csv")
    # Identify the region most often ranked first by absolute regional contribution.
    common_region = focal.largest_region5.value_counts().idxmax()
    plants_here = focal.index[focal.largest_region5.eq(common_region)]
    local = rows[rows.region5.eq(common_region) & rows.plant.isin(plants_here)]
    radiance = local.groupby(["region5", "plant", "land_class"]).viirs_median.agg(["size", "median", "min", "max"])
    radiance.to_csv(out / "common_region_radiance.csv")
    jointly_recorded = local.merge(local, on="key", suffixes=("_a", "_b"))
    pairs = jointly_recorded[jointly_recorded.plant_a.lt(jointly_recorded.plant_b)]
    print("Co-photographed in the most frequent region:", pairs.groupby(["plant_a", "plant_b"]).size().to_dict())
    print(focal[columns].round(3).to_string())
    print("Absolute contribution fraction:", focal.gap_pp.abs().sum() / ledger.gap_pp.abs().sum())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "signed", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    run(args.inputs, args.signed, args.output)
