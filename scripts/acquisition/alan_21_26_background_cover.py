"""Recover annual land cover at all original 21_26 cells using aligned 0.001-degree exports."""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import rasterio
from rasterio.io import MemoryFile

SERVICE = "https://ic.imagery1.arcgis.com/arcgis/rest/services/Sentinel2_10m_LandCover/ImageServer"
CLASSES = {1: "water", 2: "trees", 4: "flooded", 5: "crops", 7: "built", 8: "bare",
           9: "snow", 10: "clouds", 11: "range"}


def download(tile, directory, ids):
    year, west, south = tile
    path = directory / f"{year}_{west}_{south}.tif"
    if not path.exists():
        params = dict(f="image", bbox=f"{west},{south},{west+1},{south+1}", bboxSR=4326,
                      imageSR=4326, size="1000,1000", format="tiff", pixelType="U8",
                      interpolation="RSP_NearestNeighbor", renderingRule=json.dumps(dict(rasterFunction="None")),
                      mosaicRule=json.dumps(dict(mosaicMethod="esriMosaicLockRaster", lockRasterIds=[ids[year]])))
        r = requests.get(SERVICE + "/exportImage", params=params, timeout=180)
        r.raise_for_status()
        with MemoryFile(r.content) as mem, mem.open() as ds:
            assert ds.count == 1 and ds.width == ds.height == 1000
        tmp = path.with_suffix(".part")
        tmp.write_bytes(r.content)
        tmp.replace(path)
    return path


def run(frame, directory, output):
    directory.mkdir(exist_ok=True, parents=True)
    x = pd.read_parquet(frame, columns=["key", "region5", "year", "land_class"]).query("region5=='21_26'").drop_duplicates().copy()
    assert x.key.nunique() == len(x) == 835
    x[["gx", "gy"]] = x.key.str.split("|", regex=False).str[1].str.split("_", expand=True).astype(int)
    x["west"], x["south"] = x.gx // 20 - 180, x.gy // 20 - 90
    groups = x.groupby(["year", "west", "south"], sort=True)
    receipt = directory / "attributes.json"
    if not receipt.exists():
        resp = requests.get(SERVICE + "/query", params=dict(where="1=1", outFields="OBJECTID,Name,Year",
                            returnGeometry="false", f="json"), timeout=45)
        resp.raise_for_status()
        receipt.write_text(json.dumps(resp.json(), indent=2), encoding="utf-8")
    attributes = [item["attributes"] for item in json.loads(receipt.read_text())["features"]]
    ids = {a["Year"]: a["OBJECTID"] for a in attributes if str(a["Year"]) in a["Name"]}
    assert set(x.year).issubset(ids)
    tiles = list(groups.groups)
    with ThreadPoolExecutor(max_workers=6) as pool:
        paths = list(pool.map(lambda t: download(t, directory, ids), tiles))
    records, hashes = [], []
    for tile, path in zip(tiles, paths):
        year, west, south = tile
        with rasterio.open(path) as ds:
            assert ds.crs.to_epsg() == 4326
            np.testing.assert_allclose([*ds.bounds], [west,south,west+1,south+1], atol=1e-8, rtol=0)
            img = ds.read(1)
        assert set(np.unique(img)) <= (set(CLASSES) | {0, 255})
        weights = np.cos(np.deg2rad(south + 1 - (np.arange(1000) + .5) / 1000))
        for cell in groups.get_group(tile).itertuples():
            row, col = 50*(19-cell.gy % 20), 50*(cell.gx % 20)
            block = img[row:row+50, col:col+50]
            area = np.bincount(block.ravel(), weights=np.repeat(weights[row:row+50],50), minlength=256)
            valid = sum(area[c] for c in CLASSES if c != 10)
            r = dict(key=cell.key, year=year, land_class=cell.land_class, valid_frac=valid/area.sum())
            for code, label in CLASSES.items():
                if code != 10:
                    r[label+"_frac"] = area[code]/valid if valid else np.nan
            records.append(r)
        hashes.append(dict(year=year, west=west, south=south, objectid=ids[year], sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    output.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(records).sort_values("key").to_csv(output/"background_landcover.csv", index=False)
    pd.DataFrame(hashes).to_csv(output/"background_tile_receipts.csv", index=False)
    print("Annual assets:", ids, "tiles:", len(tiles), "summarized rows:", len(records))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("frame", "cache", "output"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    run(a.frame, a.cache, a.output)
