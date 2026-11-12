"""Describe conditional plant composition slopes within 21_26, separately from the global attribution."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from alan_gecko_primary_fiveyear import primary
from alan_shared_species_contributions import influence

FOCAL = ["Asclepias syriaca", "Cirsium vulgare"]
BASE = ["island", "loglight", "light_island", "logpop", "logphotos", "logspecies",
        "bio5", "lograin", "elevation_km"]
COVER = ["built_frac", "open_frac", "water_frac"]
SEASON = ["season_sin", "season_cos"]


def summaries(inputs, signed, photos, out):
    s, _, members, _ = primary(inputs, signed)
    s = s[s.region5.eq("21_26")].copy()
    assert len(s) == 835
    cover = pd.read_csv(out / "background_landcover.csv")
    s = s.merge(cover, on=["key", "year", "land_class"], validate="1:1")
    assert len(s) == 835 and s.valid_frac.min() > .99
    s["open_frac"] = s.bare_frac + s.crops_frac + s.range_frac
    dates = pd.read_parquet(photos, columns=["region5", "key", "observed_on"]).query("region5=='21_26'")
    dates = dates[["key", "observed_on"]].drop_duplicates().copy()
    dt = pd.to_datetime(dates.observed_on)
    angle = 2 * np.pi * (dt.dt.dayofyear - 1) / (365 + dt.dt.is_leap_year.astype(int))
    dates["season_sin"], dates["season_cos"] = np.sin(angle), np.cos(angle)
    season = dates.groupby("key")[SEASON].mean()
    s = s.join(season, on="key", validate="1:1")
    assert s[BASE + COVER + SEASON].notna().all().all()
    groups = s.groupby(["region_year", "land_class"]).key.transform("size")
    weights = groups.groupby(s.region_year).transform("min") / groups
    vals = members[members.wcvp_accepted_name.isin(FOCAL)]
    y = np.column_stack([100 * s.key.isin(vals.loc[vals.wcvp_accepted_name.eq(t), "key"]) / s.n_species for t in FOCAL])
    assert (y > 0).sum(axis=0).tolist() == [25, 13]
    years = pd.get_dummies(s.year, dtype=float)
    rows = []
    for name, extra in [("local_base", []), ("local_cover", COVER), ("local_cover_and_dates", COVER + SEASON)]:
        cols = BASE + extra
        x = pd.concat([s[cols], years], axis=1).to_numpy()
        w = np.sqrt(weights.to_numpy())[:, None]
        b, _, rank, _ = np.linalg.lstsq(x * w, y * w, rcond=1e-12)
        assert rank == x.shape[1]
        if not extra:
            # The fixed-effect implementation used to attribute the published model.
            np.testing.assert_allclose(b[[1, 2], :], influence(s) @ y, rtol=1e-9, atol=1e-10)
        for j, plant in enumerate(FOCAL):
            main, gap = b[cols.index("loglight"), j], b[cols.index("light_island"), j]
            rows.append(dict(specification=name, plant=plant, island=main + gap, mainland=main,
                             island_minus_mainland=gap, n=len(s), region_clusters=1, positive_island=int(
                                 np.sum((y[:, j] > 0) & s.island.astype(bool))), positive_mainland=int(
                                 np.sum((y[:, j] > 0) & ~s.island.astype(bool)))))
    pd.DataFrame(rows).to_csv(out / "background_local_slopes.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print("Within-side rank correlations between VIIRS and built fraction:")
    print(s.groupby("land_class")[["viirs_median", "built_frac"]].corr(method="spearman"))
    print("Distinct observation dates per cell:", dates.groupby("key").size().describe())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("inputs", "signed", "photos", "output"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    summaries(a.inputs, a.signed, a.photos, a.output)
