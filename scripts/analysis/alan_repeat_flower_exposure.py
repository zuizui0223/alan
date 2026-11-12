"""Assign monitoring coordinates to the primary region-years and sample same-year VIIRS."""
import argparse
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

from mois_global_2020_GSHHG_scope_20261107 import assign_points, prepare_land
from mois_global_multipoint_numeric_20261107 import MONTHS, VIIRS, api_params, daterange, extract_samples


def expose(locations, frame, coast, cache):
    cache.mkdir(exist_ok=True, parents=True)
    cols = ["latitude", "longitude"]
    points = locations[cols + ["year"]].dropna().drop_duplicates().copy()
    points["region5"] = (np.floor((points.longitude + 180) / 5).astype(int).astype(str) + "_" +
                         np.floor((points.latitude + 90) / 5).astype(int).astype(str))
    sites = frame.drop_duplicates("key")
    assert len(sites) == 10698
    thresholds = sites.groupby(["region5", "year"]).viirs_median.median().rename("threshold")
    assert len(thresholds) == 195
    points = points.merge(thresholds, on=["region5", "year"], validate="m:1")
    unique = points[cols].drop_duplicates().reset_index(drop=True)
    unique["unit_id"] = unique.index
    land = assign_points(unique, prepare_land(gpd.read_file(coast / "GSHHS_h_L1.shp")),
                         gpd.read_file(coast / "GSHHS_h_L2.shp"))
    points = points.merge(unique.merge(land, on="unit_id"), on=cols, validate="m:1")
    points = points[points.land_class.isin(["mainland","island"])].sort_values(["year", *cols]).reset_index(drop=True)
    points["key"] = (points.year.astype(str) + "|" +
                     np.floor((points.longitude + 180) * 20).astype(int).astype(str) + "_" +
                     np.floor((points.latitude + 90) * 20).astype(int).astype(str))
    plans = [(int(y), m, g.iloc[k:k+96]) for y,g in points.groupby("year") for m in MONTHS
             for k in range(0,len(g),96)]

    def get(plan):
        year, month, group = plan
        rows = group.rename(columns={"longitude":"query_longitude","latitude":"query_latitude"}).to_dict("records")
        params = api_params("viirs", month, rows)
        params["time"] = daterange(year, month)
        tag = hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()
        path = cache / (tag + ".json")
        for attempt in range(4):
            try:
                if path.exists():
                    raw = path.read_bytes()
                else:
                    response = requests.get(VIIRS, params=params, timeout=60)
                    response.raise_for_status()
                    raw = response.content
                parsed = extract_samples("viirs", month, rows, json.loads(raw), year=year)
                if not path.exists():
                    path.write_bytes(raw)
                break
            except (requests.RequestException,ValueError):
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)
        return pd.Series([r["value"] for r in parsed], index=group.index, name=month)
    if len(points):
        with ThreadPoolExecutor(max_workers=4) as pool:
            monthly = pd.concat([v.to_frame().rename_axis("row").reset_index().melt("row")
                                 for v in pool.map(get, plans)])
        measured = monthly.pivot(index="row",columns="variable",values="value").apply(pd.to_numeric)
        points["valid_months"] = measured.notna().sum(axis=1)
        points["radiance"] = measured.median(axis=1).where(points.valid_months.ge(3))
        points["light"] = np.where(points.radiance.gt(points.threshold),"bright","dark")
        points.loc[points.radiance.isna(),"light"] = pd.NA
    return points


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["locations", "photo_frame", "coast", "cache", "output"]:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    expose(pd.read_parquet(args.locations), pd.read_parquet(args.photo_frame), args.coast, args.cache).to_parquet(args.output,index=False)
