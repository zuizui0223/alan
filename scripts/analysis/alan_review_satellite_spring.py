"""Read first-cycle vegetation mid-greenup at the six 2021 Holodiscus stations."""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

YEARS = (2019, 2021, 2023, 2024, 2025)
API = "https://modis.ornl.gov/rst/api/v1/MCD12Q2/subset"


def greenup_date(bands):
    day = bands.get("MidGreenup.Num_Modes_01", 32767)
    qa = bands.get("QA_Detailed.Num_Modes_01", 65535)
    cycles = bands.get("NumCycles", 0)
    if not (1 <= cycles <= 7 and 11138 <= day <= 32766 and 0 <= qa < 65535):
        return pd.NaT
    if ((int(qa) >> 2) & 3) > 1:
        return pd.NaT
    return pd.Timestamp("1970-01-01") + pd.Timedelta(days=day)


def summarize(visits, exposure, cache, output):
    cache.mkdir(parents=True, exist_ok=True)
    plant = pd.read_parquet(visits).query("plant == 'Holodiscus discolor' and year == 2021")
    plant = plant.merge(pd.read_parquet(exposure).query("region5 == '11_27'"),
                        on=["year", "latitude", "longitude"], validate="m:1")
    stations = plant[["station", "latitude", "longitude", "land_class", "light"]].drop_duplicates()
    assert len(stations) == stations.station.nunique() == 6

    def extract(station):
        path = cache / f"{station.station}.json"
        if not path.exists():
            params = dict(latitude=station.latitude, longitude=station.longitude,
                          startDate="A2019001", endDate="A2025001", kmAboveBelow=0, kmLeftRight=0)
            response = requests.get(API, params=params, timeout=120)
            response.raise_for_status()
            path.write_bytes(response.content)
        raw = path.read_bytes()
        source = json.loads(raw)
        assert source["nrows"] == source["ncols"] == 1
        for key in ("latitude", "longitude"):
            assert abs(source[key] - getattr(station, key)) < 1e-5
        rows = []
        for year in YEARS:
            sub = [r for r in source["subset"] if r["modis_date"] == f"A{year}001"]
            bands = {r["band"]: pd.to_numeric(r["data"][0], errors="coerce") for r in sub}
            rows.append(dict(year=year, land_class=station.land_class, light_in_2021=station.light,
                             year_returned=bool(sub), date=greenup_date(bands)))
        return rows, hashlib.sha256(raw).hexdigest()

    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(extract, stations.itertuples(index=False)))
    frame = pd.DataFrame([row for rows, _ in results for row in rows])
    table = frame.groupby(["year", "land_class", "light_in_2021"]).agg(
        stations=("year_returned", "size"), with_source=("year_returned", "sum"),
        good_or_best=("date", "count"), earliest_midgreenup=("date", "min"),
        latest_midgreenup=("date", "max"))
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "holodiscus_satellite_spring.csv")
    provenance = dict(product="MCD12Q2", service=API,
                      private_response_sha256=sorted(d for _, d in results))
    (output / "holodiscus_satellite_source.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(table.to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["visits", "exposure", "cache", "output"]:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    summarize(args.visits, args.exposure, args.cache, args.output)
