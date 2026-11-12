"""Join original focal plant observations to same-year Esri 10 m land cover.

Usage: python scripts/analysis/alan_21_26_landcover.py WORK_DIRECTORY OUTPUT_DIRECTORY
Public outputs omit observation IDs and exact coordinates. The private cache retains
the source TIFFs and per-observation neighbourhoods for reproducing seasonal joins.
"""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import rasterio
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

SERVICE = "https://ic.imagery1.arcgis.com/arcgis/rest/services/Sentinel2_10m_LandCover/ImageServer"
REGION = "21_26"
SPECIES = ("Asclepias syriaca", "Cirsium vulgare")
YEARS = (2019, 2021, 2023, 2024, 2025)
# Cloud, source nodata and unmapped pixels excluded from valid-area fractions.
CLASSES = {1: "water", 2: "trees", 4: "flooded_vegetation", 5: "crops",
           7: "built", 8: "bare", 9: "snow_ice", 11: "rangeland"}
PROJECT = Transformer.from_crs(4326, 5070, always_xy=True)


def session():
    http = requests.Session()
    retry = Retry(total=4, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
    http.mount("https://", HTTPAdapter(max_retries=retry))
    return http


def load_observations(work):
    photos = pd.read_parquet(work / "phenology_run/target.parquet")
    photos = photos[photos.region5.eq(REGION) & photos.plant.isin(SPECIES)]
    coordinates = pd.read_csv(work / "mois/signed/canonical/canonical_observations_climate.csv.gz",
                              usecols=["photo_id", "observation_id", "longitude", "latitude",
                                       "source_subarm", "observed_on"])
    coordinates = coordinates[coordinates.source_subarm.eq("historical_main")]
    photos = photos.merge(coordinates, on=["photo_id", "observation_id", "observed_on"],
                          how="left", validate="1:1")
    assert not photos[["longitude", "latitude"]].isna().any().any()
    assert set(photos.year) == set(YEARS)
    keys = (photos.year.astype(str) + "|" + np.floor((photos.longitude+180)*20).astype(int).astype(str)
            + "_" + np.floor((photos.latitude+90)*20).astype(int).astype(str))
    assert photos.key.eq(keys).all()
    assert photos.groupby("observation_id")[["longitude", "latitude", "plant", "key"]].nunique().max().max() == 1
    return photos.drop_duplicates("observation_id").copy()


def projected_cell(key, margin):
    gx, gy = map(int, key.split("|")[1].split("_"))
    west, south = gx/20-180, gy/20-90
    east, north = west+.05, south+.05
    corners = [(west,south), (east,south), (east,north), (west,north), (west,south)]
    ring = []
    for left, right in zip(corners, corners[1:]):
        xs, ys = PROJECT.transform(np.linspace(left[0],right[0],21), np.linspace(left[1],right[1],21))
        ring.extend(zip(xs,ys))
    xs, ys = np.array(ring).T
    # Export enough margin for the full neighbourhood of observations on edges.
    bbox = (10*int(np.floor((xs.min()-margin)/10)), 10*int(np.floor((ys.min()-margin)/10)),
            10*int(np.ceil((xs.max()+margin)/10)), 10*int(np.ceil((ys.max()+margin)/10)))
    return {"type": "Polygon", "coordinates": [ring]}, bbox


def class_fractions(values):
    valid = np.isin(values, list(CLASSES))
    nv = int(valid.sum())
    out = {"valid_pct": 100*nv/len(values)}
    out.update({name+"_pct": 100*np.count_nonzero(values == code)/nv if nv else np.nan
                for code, name in CLASSES.items()})
    return out


def extract(group, raster_ids, cache):
    key, rows = group
    year = int(rows.year.iat[0])
    obj = raster_ids[year]
    max_radius = int(rows.radius_m.max()) if "radius_m" in rows else 100
    polygon, bbox = projected_cell(key, max_radius+10)
    west, south, east, north = bbox
    width, height = (east-west)//10, (north-south)//10
    suffix = "" if max_radius == 100 else f"_{max_radius}m"
    tif = cache / (key.replace("|","_") + suffix + ".tif")
    params = dict(f="image", format="tiff", bbox=",".join(map(str,bbox)), bboxSR=5070, imageSR=5070,
                  size=f"{width},{height}", pixelType="U8", interpolation="RSP_NearestNeighbor",
                  mosaicRule=json.dumps({"mosaicMethod":"esriMosaicLockRaster","lockRasterIds":[obj]}),
                  renderingRule='{"rasterFunction":"None"}', adjustAspectRatio="false")
    if not tif.exists():
        with session() as http:
            response = http.get(SERVICE+"/exportImage", params=params, timeout=120)
            response.raise_for_status()
            if "tiff" not in response.headers.get("Content-Type","").lower():
                raise ValueError(response.text[:400])
            temp = tif.with_suffix(".tmp")
            temp.write_bytes(response.content)
            temp.replace(tif)
    with rasterio.open(tif) as source:
        assert source.crs.to_epsg() == 5070 and source.count == 1
        assert (source.width,source.height) == (width,height)
        transform = from_origin(west,north,10,10)
        assert source.transform.almost_equals(transform)
        image = source.read(1)
    if not set(np.unique(image)).issubset({0,10,255,*CLASSES}):
        raise ValueError(f"Unexpected land-cover category in {key}")
    mask = geometry_mask([polygon], image.shape, transform, invert=True)
    cell = class_fractions(image[mask])
    xs = west+10*(np.arange(width)+0.5)
    ys = north-10*(np.arange(height)+0.5)
    observations = []
    for row in rows.itertuples(index=False):
        x,y = PROJECT.transform(row.longitude,row.latitude)
        distance = getattr(row, "radius_m", 100)
        radius = (xs-x)**2 + (ys[:,None]-y)**2 <= distance**2
        assert 0.95 <= radius.sum() / (np.pi*distance**2/100) <= 1.05
        observations.append({"observation_id":row.observation_id, "key":key, "plant":row.plant,
                             "year":year, "land_class":row.land_class, **class_fractions(image[radius])})
    manifest = {"key":key, "year":year, "raster_id":obj, "sha256":hashlib.sha256(tif.read_bytes()).hexdigest(),
                "cell_pixels":int(mask.sum()), "cell_valid_pct":cell["valid_pct"]}
    return {"key":key, **{"cell_"+k:v for k,v in cell.items()}}, observations, manifest



def island_accuracy(work, observations, ids, cache, output):
    context = json.loads((work / "local_21_26_private/photo_context.json").read_bytes())
    islands = observations[observations.land_class.eq("island")].copy()
    assert len(islands) == 10
    api_refs = {p["id"]: row for row in context["results"] for p in row["photos"]}
    uuids = sorted({api_refs[pid]["uuid"] for pid in islands.photo_id})
    api_file = cache.parent/"island_public_positions.json"
    if not api_file.exists():
        url = "https://api.inaturalist.org/v2/observations/" + ",".join(uuids)
        response = session().get(url, params={"fields": "uuid,geojson,positional_accuracy,geoprivacy,obscured,photos.id,observed_on"}, timeout=60)
        response.raise_for_status()
        api_file.write_bytes(response.content)
    data = json.loads(api_file.read_bytes())["results"]
    assert sorted(row["uuid"] for row in data) == uuids
    positions = {p["id"]: row for row in data for p in row["photos"]}
    evidence = []
    for row in islands.itertuples(index=False):
        pos = positions[row.photo_id]
        assert pos["observed_on"] == row.observed_on and pos["uuid"] == api_refs[row.photo_id]["uuid"]
        assert pos["positional_accuracy"] == api_refs[row.photo_id]["positional_accuracy"]
        assert pos["geoprivacy"] in (None, "open") and pos["obscured"] is False
        x,y = PROJECT.transform(*pos["geojson"]["coordinates"])
        ox,oy = PROJECT.transform(row.longitude,row.latitude)
        offset = float(np.hypot(x-ox, y-oy))
        evidence.append({"observation_id": row.observation_id, "public_position_offset_m":offset,
                         "reported_accuracy_m":pos["positional_accuracy"]})
    islands = islands.merge(pd.DataFrame(evidence), on="observation_id",validate="1:1")
    islands["radius_m"] = islands.reported_accuracy_m.clip(lower=100)
    with ThreadPoolExecutor(max_workers=3) as threads:
        records = list(threads.map(lambda g: extract(g,ids,cache), islands.groupby("key")))
    cover = pd.DataFrame([r for _,rows,_ in records for r in rows])
    index = ["observation_id", "key", "plant"]
    cover = cover.merge(islands[index+["radius_m","reported_accuracy_m","public_position_offset_m"]],on=index)
    assert len(cover) == 10 and not cover.duplicated(["key","plant"]).any()
    cover.drop(columns="observation_id").to_csv(output/"island_accuracy_cover.csv",index=False)
    pd.DataFrame([m for _,_,m in records]).to_csv(output/"island_accuracy_rasters.csv",index=False)
    print("Maximum difference from current public position (m):",cover.public_position_offset_m.max())
    print(cover.groupby("plant")[["built_pct","rangeland_pct","water_pct"]].median())


def run(work, output):
    private = work / "landcover_21_26"
    cache = private / "raster"
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    docs = session().get(SERVICE+"/query", params={"where":"1=1","f":"json","outFields":"OBJECTID,Year,Name",
                         "returnGeometry":"false"},timeout=45).json()
    catalog = pd.DataFrame([r["attributes"] for r in docs["features"]]).set_index("Year")
    assert catalog.index.is_unique
    raster_ids = {year:int(catalog.loc[year,"OBJECTID"]) for year in YEARS}
    for year in YEARS:
        assert str(year) in catalog.loc[year,"Name"]
    observations = load_observations(work)
    with ThreadPoolExecutor(max_workers=6) as threads:
        records = list(threads.map(lambda g: extract(g,raster_ids,cache), observations.groupby("key")))
    cells = pd.DataFrame([c for c,_,_ in records])
    photos = pd.DataFrame([r for _,rows,_ in records for r in rows])
    group = ["key","plant","year","land_class"]
    means = photos.groupby(group)[[c for c in photos if c.endswith("_pct")]].mean().add_prefix("photo_")
    means["observations"] = photos.groupby(group).size()
    table = means.reset_index().merge(cells, on="key", validate="m:1")
    assert len(table)==38 and table.key.nunique()==37 and photos.observation_id.is_unique
    table.to_csv(output/"landcover.csv", index=False)
    pd.DataFrame([m for _,_,m in records]).to_csv(output/"landcover_raster_sources.csv",index=False)
    photos.to_csv(private/"observation_cover.csv",index=False)
    island_accuracy(work, observations, raster_ids, cache, output)
    print(table.groupby(["plant","land_class"])[["observations","cell_built_pct","photo_built_pct",
                                                  "cell_rangeland_pct","photo_rangeland_pct"]].median())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.work,args.output)
