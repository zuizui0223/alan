"""Compare photographed conspecifics at brighter and dimmer sites on the same island in the same year."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


def load_waves(directory):
    receipt = json.loads(next(directory.rglob("*signed_receipt.json")).read_text())
    frames = []
    for wave, digest in receipt["wave_output_SHA256"].items():
        path = next(directory.rglob(f"*wave{wave}.csv.gz"))
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
        frame = pd.read_csv(path)
        assert len(frame) == 8192 and frame.wave.eq(int(wave)).all()
        frames.append(frame)
    photos = pd.concat(frames, ignore_index=True)
    assert len(photos) == photos.photo_id.nunique() == photos.observation_id.nunique() == 24576
    assert photos.strict9.equals(photos.strict9_white.notna())
    assert photos.strict9.sum() == 8803
    assert photos.strict9_white.dropna().between(0, 1).all()
    np.testing.assert_allclose(photos.log_light, np.log1p(photos.viirs_median), atol=1e-7)
    return photos


def compare(photos, outcome, ratio):
    keys = ["land_id", "year", "wcvp_accepted_plant_name_id"]
    assert photos[keys + ["target_id", "log_light"]].notna().all().all()
    assert photos.groupby("target_id").log_light.nunique().max() == 1
    # Select extremes on exposure before checking colour, retaining all tied sites.
    sites = photos.groupby(keys + ["target_id"], as_index=False).agg(
        light=("log_light", "first"), response=(outcome, "mean"))
    strata = sites.groupby(keys).light
    low = strata.transform("min")
    high = strata.transform("max")
    supported = (high > low) & (high - low >= np.log(ratio))
    endpoints = sites.loc[supported & ((sites.light == low) | (sites.light == high))].copy()
    endpoints["bright"] = endpoints.light.eq(high.loc[endpoints.index])
    pairs = endpoints.groupby(keys + ["bright"]).response.mean().unstack()
    complete = pairs.dropna()
    contrasts = (complete[True] - complete[False]) * 100
    islands = contrasts.groupby(level="land_id").mean()
    n = len(islands)
    assert n > 1
    estimate = islands.mean()
    margin = t.ppf(.975, n - 1) * islands.sem()
    return dict(outcome=outcome, min_ratio_one_plus_radiance=ratio,
                eligible_species_island_years=len(pairs), complete_pairs=len(complete),
                species=complete.index.get_level_values(keys[-1]).nunique(), islands=n,
                brighter_minus_dimmer_pp=estimate, ci_low=estimate-margin, ci_high=estimate+margin)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = load_waves(args.source)
    results = [compare(data, outcome, ratio) for ratio in (2, 1, 4)
               for outcome in ("strict9_white", "strict9")]
    pd.DataFrame(results).to_csv(args.output, index=False)
