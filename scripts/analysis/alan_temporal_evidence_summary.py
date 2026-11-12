"""Summarize dated source support and measured VIIRS without fitting new models."""
import argparse
from itertools import combinations
from pathlib import Path

import pandas as pd


def summarize(cache):
    sites = pd.read_parquet(cache / "universe.parquet").drop_duplicates("key")
    land = pd.read_parquet(cache / "strict_site_light.parquet", columns=["key", "land_id"])
    sites = sites.merge(land, on="key", validate="1:1")
    sites["cell"] = sites.key.str.split("|", regex=False).str[-1]
    assert len(sites) == 10698 and sites.region5.nunique() == 78
    assert set(sites.year) == {2019, 2021, 2023, 2024, 2025}
    years = sites.year.nunique()
    full = sites[sites.groupby(["land_id", "cell"]).year.transform("nunique").eq(years)]
    radiance = full.pivot(index=["land_id", "cell", "land_class"], columns="year", values="viirs_median")
    assert radiance.shape == (113, years)
    print("Median VIIRS by year in the same 113 cells:")
    print(radiance.groupby(level="land_class").median())
    change = radiance[2025] - radiance[2019]
    print("Within-cell 2025 minus 2019, by land:")
    print(change.groupby(level="land_class").agg(["size", "median", lambda y: (y > 0).sum()]))

    target = pd.read_parquet(cache / "target.parquet", columns=["key", "plant"])
    assert len(target) == 8466 and target.plant.nunique() == 544
    target = target.drop_duplicates().merge(
        sites[["key","year","region5","land_class","land_id","cell"]], on="key", validate="m:1"
    )
    paired = target.groupby(["year","region5","plant"]).land_class.nunique().eq(2)
    print("Species seen on both sides within the same region-year:")
    print(paired[paired].reset_index().groupby("year").plant.nunique())
    groups = ["region5", "plant", "land_class", "cell", "land_id"]
    detected = target.groupby(groups).year.agg(lambda x: sorted(set(x)))
    repeated = detected[detected.str.len().ge(2)].reset_index()
    by_land = repeated.groupby(["region5","plant"]).land_class.nunique()
    candidates = repeated.merge(by_land[by_land.eq(2)].rename("classes"), on=["region5","plant"])
    print("Both sides have a cell with the species recorded at least twice:")
    print(candidates[["region5","plant","land_class","year"]].to_string(index=False))
    candidates["year_pair"] = candidates.year.map(lambda ys: list(combinations(ys,2)))
    pairs = candidates.explode("year_pair").groupby(["region5","plant","year_pair"]).land_class.nunique()
    print("Identical year pair on both sides:", int(pairs.eq(2).sum()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    summarize(parser.parse_args().cache)
