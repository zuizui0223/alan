#!/usr/bin/env python3
"""Source-verified local dark-refuge structure around real ALAN island photo sites.

Pilot: 2025 botanical photo site-years, frozen 2024 PCNL v3 (PRE-observation),
original GSHHG 2.3.7 island polygons (inland lakes removed), actual raster
pixels. This only measures the *geography of satellite-detected darkness*.
It does NOT measure habitat, insect movement, ground lux, moth abundance,
effective pollination, flowering change, or a causal ALAN response.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from rasterio.windows import Window, from_bounds
from scipy import ndimage
from shapely import box, contains_xy, union_all

from scripts.analysis.alan_island_area_light_20261108 import source_land_area
from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land

RECORD = "https://zenodo.org/records/17013414"
ANALYSIS_YEAR = 2025
RASTER_YEAR = 2024
THRESHOLDS = (0.0, 1.0, 3.0)
GEOD = Geod(ellps="WGS84")
SOURCE_N = 85318
PREFIX = "ALAN_DARK_REFUGIA: "
SCHEMA = "alan_dark_refugia_original_2025_sites_prior_2024_PCNL_GSHHG_v1"


def require(ok, message):
    if not bool(ok):
        raise ValueError(PREFIX + str(message))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_2025_sites(joined, receipt, max_per_island):
    """Use no biological outcome when selecting the deterministic pilot."""
    report = json.loads(Path(receipt).read_text())
    require(report["schema"] == "ALAN_full93m_real_PCNL_33year_1992_2024_source_join_v1",
            "Wrong previously verified historical source")
    require(report["source_site_year_rows"] == SOURCE_N and
            report["original_source_sha256"] ==
            "e7800476ae6d35acb9ef71525f32521596f884733e2e8ee9f87fda7acb8e9eab",
            "Original 93m-photo site frame provenance differs")
    require(sha256(joined) == report["output_sha256"],
            "Previously signed real PCNL join bytes changed")
    d = pd.read_csv(joined, compression="infer", dtype={"island":str}, low_memory=False)
    need = {"island", "latbin", "lonbin", "year", "lat", "lon",
            "historic_ntl_complete", "all_accepted_observation_uuids",
            "all_accepted_recorded_species"}
    require(len(d) == SOURCE_N and need.issubset(d.columns),
            "Expected true 85,318 distinct original island site-years")
    require(not d.duplicated(["island", "latbin", "lonbin", "year"]).any(),
            "Duplicate island/site/year")
    require(max_per_island > 0, "Per-island pilot cap must be positive")
    d = d.loc[(d.year == ANALYSIS_YEAR) &
              d.historic_ntl_complete.astype(str).str.lower().eq("true")].copy()
    require(len(d) and d.island.notna().all() and
            d.lat.between(-90, 90).all() and
            d.lon.between(-180, 180).all(), "Invalid pilot coordinates")
    d["pilot_sort"] = [
        hashlib.sha256(f"{i}|{a}|{b}|2025".encode()).hexdigest()
        for i, a, b in zip(d.island, d.latbin, d.lonbin)
    ]
    # Deterministic, response-blind, capped per island; NOT random population sampling.
    d = (d.sort_values(["island", "pilot_sort"], kind="stable")
           .groupby("island", sort=False, group_keys=False)
           .head(max_per_island).drop(columns="pilot_sort").reset_index(drop=True))
    require(not d.duplicated(["island", "latbin", "lonbin", "year"]).any(),
            "Pilot photo cell duplicated")
    return d


def dark_refuge_metrics(src, island_polygon, lon, lat, radius_km=5.,
                        thresholds=THRESHOLDS, min_valid_land_pixels=16,
                        min_coverage=0.80):
    """Raster-native 4-neighbour dark-component metrics in a geodesic buffer.

    Ocean and inland lakes are excluded by caller's source-island polygon.
    Dark = calibrated PCNL <= fixed threshold; *not true ecological darkness*.
    Areas use a latitude-corrected pixel-weight approximation; small radius
    and the same product are essential for spatial comparability.
    """
    require(src.crs.to_epsg() == 4326 and src.count == 1,
            "Require one original global EPSG:4326 annual PCNL band")
    require(np.isfinite([lon, lat, radius_km]).all() and
            -180.0 <= lon <= 180.0 and -80.0 <= lat <= 80.0 and
            0 < radius_km <= 25, "Unsupported geographic coordinates/radius")
    require(abs(lon) <= 179.5, "Antimeridian requires separate wrapped sampling")
    require(0 < min_coverage <= 1 and min_valid_land_pixels >= 1,
            "Invalid land coverage gate")
    require(bool(island_polygon.is_valid) and not island_polygon.is_empty,
            "Unverified original island polygon")
    bounds_lat = radius_km / 110.57
    bounds_lon = radius_km / (111.32 * np.cos(np.deg2rad(lat)))
    win = from_bounds(lon - bounds_lon, lat - bounds_lat,
                      lon + bounds_lon, lat + bounds_lat, src.transform)
    # Add one native pixel outside the approximate geodesic bounding box.
    win = Window(int(np.floor(win.col_off)) - 1, int(np.floor(win.row_off)) - 1,
                 int(np.ceil(win.width)) + 3, int(np.ceil(win.height)) + 3)
    require(win.width * win.height <= 100000, "Polar or raster-scale window too large")
    values = src.read(1, window=win, boundless=True, masked=True)
    # Convert integer (e.g. original VIIRS Int16) masked values to float
    # BEFORE filling with NaN. Otherwise NumPy cannot fill NaN into Int16.
    raw = values.astype(float).filled(np.nan)
    transform = src.window_transform(win)
    rr, cc = np.indices(raw.shape)
    xx = transform.c + (cc + .5) * transform.a + (rr + .5) * transform.b
    yy = transform.f + (cc + .5) * transform.d + (rr + .5) * transform.e
    _, _, dm = GEOD.inv(np.full(xx.shape, lon), np.full(xx.shape, lat), xx, yy)
    within = dm <= 1000.0 * radius_km
    land = within & contains_xy(island_polygon, xx, yy)
    # PCNL nodata/invalid pixels are unknown, not 'dark'.
    valid = land & np.isfinite(raw) & (raw >= 0)
    n_land, n_valid = int(land.sum()), int(valid.sum())
    common = {
        "radius_km": radius_km,
        "n_original_island_land_pixels": n_land,
        "n_valid_original_PCNL_pixels": n_valid,
        "valid_land_fraction": n_valid / n_land if n_land else None,
    }
    if n_valid < min_valid_land_pixels or not n_land or \
            n_valid / n_land < min_coverage:
        return [{**common, "threshold_PCNL_native": float(t),
                 "status": "insufficient_original_raster_land_support"}
                for t in thresholds]
    # Approximate per-pixel geographic surface area, not an ecological weighting.
    weights = np.maximum(np.cos(np.deg2rad(yy)), 0)
    denom = float(weights[valid].sum())
    require(denom > 0, "Non-positive geodetic land area")
    mean_light = float(np.sum(raw[valid] * weights[valid]) / denom)
    result = []
    for t in thresholds:
        require(np.isfinite(t) and t >= 0, "Darkness threshold must be fixed numeric")
        dark = valid & (raw <= t)
        dark_area = float(weights[dark].sum())
        labels, count = ndimage.label(dark, structure=np.array(
            [[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8))
        if count:
            component_area = np.bincount(
                labels[dark], weights=weights[dark], minlength=count + 1)
            biggest = float(component_area[1:].max())
            nearest_km = float(dm[dark].min()) / 1000
        else:
            biggest = 0.0
            nearest_km = None  # right-censored above radius, not a zero-distance.
        dark_fraction = dark_area / denom
        result.append({
            **common, "threshold_PCNL_native": float(t),
            "status": "observed_satellite_geometry_only",
            "dark_land_fraction": dark_fraction,
            "largest_dark_component_fraction_of_land": biggest / denom,
            "largest_dark_component_fraction_of_dark":
                biggest / dark_area if dark_area else None,
            "n_four_connected_dark_components": int(count),
            "nearest_detected_dark_pixel_km": nearest_km,
            "distance_to_dark_censored_at_radius": not bool(count),
            "local_mean_PCNL_native": mean_light,
            "n_dark_raster_pixels": int(dark.sum()),
        })
    return result


def run(joined, joined_receipt, coast, raster, out, radius, max_per_island):
    out = Path(out)
    require(not out.exists(), "Do not overwrite a real source-derived result")
    sites = source_2025_sites(joined, joined_receipt, max_per_island)
    coast = Path(coast)
    # These routines enforce original GSHHG archive SHA256 and every extracted
    # shapefile's source digest before any geometry can enter the analysis.
    land = load_signed_land(coast)
    land = land.set_index("land_id", drop=False)
    require(sites.island.astype(str).isin(land.index).all() and
            land.loc[sites.island.astype(str), "land_class"].eq("island").all(),
            "Original island identities do not match genuine GSHHG land")
    lakes = gpd.read_file(coast / "GSHHS_h_L2.shp")
    if lakes.crs is None:
        raise ValueError(PREFIX + "GSHHG lake CRS missing")
    lakes = lakes.to_crs(4326)
    index = lakes.sindex
    features = []
    skipped = []
    with rasterio.open(raster) as src:
        require(src.crs.to_epsg() == 4326, "Unexpected original PCNL CRS")
        for island, group in sites.groupby("island", sort=True):
            shape = land.loc[str(island), "geometry"]
            nearby = index.query(shape, predicate="intersects")
            if len(nearby):
                pieces = [v for v in lakes.geometry.iloc[nearby] if v.intersects(shape)]
                if pieces:
                    shape = shape.difference(union_all(pieces))
            if shape.is_empty or not shape.is_valid:
                skipped.append({"island":island, "reason":"invalid_after_inland_water_subtraction"})
                continue
            for row in group.itertuples(index=False):
                try:
                    values = dark_refuge_metrics(
                        src, shape, float(row.lon), float(row.lat), radius_km=radius)
                except ValueError as err:
                    skipped.append({"island":island, "latbin":row.latbin,
                                    "lonbin":row.lonbin, "reason":str(err)})
                    continue
                for v in values:
                    features.append({
                        "island":str(island), "latbin":row.latbin,
                        "lonbin":row.lonbin, "lat":row.lat, "lon":row.lon,
                        "photo_observation_year": ANALYSIS_YEAR,
                        "source_raster_year": RASTER_YEAR,
                        "source_all_accepted_observation_uuids":
                            row.all_accepted_observation_uuids,
                        **v,
                    })
    result = pd.DataFrame(features)
    require(not result.empty, "No island radar neighbourhoods obtained")
    require(not result.duplicated(["island", "latbin", "lonbin",
                                   "photo_observation_year",
                                   "threshold_PCNL_native"]).any(),
            "Repeated island-neighborhood / threshold is pseudoreplication")
    out.mkdir(parents=True)
    filename = out / "real_original_island_dark_refugia_2025_photo_sites_PCNL2024.csv.gz"
    result.to_csv(filename, index=False, compression="gzip")
    success = result.loc[result.status.eq("observed_satellite_geometry_only")]
    audit = {}
    for t in THRESHOLDS:
        d = success.loc[success.threshold_PCNL_native.eq(t)]
        audit[str(t)] = {
            "source_valid_photo_sites": len(d),
            "source_valid_islands": int(d.island.nunique()),
            "median_dark_land_fraction": float(d.dark_land_fraction.median()) if len(d) else None,
            "median_dark_connectivity_among_sites_with_dark":
                float(d.largest_dark_component_fraction_of_dark.median())
                if d.largest_dark_component_fraction_of_dark.notna().any() else None,
            "sites_with_zero_dark_pixels": int(d.n_dark_raster_pixels.eq(0).sum()),
            "sites_with_dark_continuity_measured":
                int(d.largest_dark_component_fraction_of_dark.notna().sum()),
        }
    report = {
        "schema":SCHEMA, "zenodo_source": RECORD,
        "product":"Li-Cao PCNL v3, original 2024 raster, native corrected value",
        "GSHHG_source":"SHA-verified 2.3.7 high-res L1 islands minus L2 inland waters",
        "source_original_observed_2025_photo_sites_selected": len(sites),
        "source_islands_selected": int(sites.island.nunique()),
        "sites_with_any_returned_metric": int(result[["island","latbin","lonbin"]].drop_duplicates().shape[0]),
        "actual_satellite_metrics_by_threshold":audit,
        "n_failed_or_unsupported_sites":len(skipped),
        "failed_reason_examples":skipped[:15],
        "radius_km":radius, "photo_observation_year": ANALYSIS_YEAR,
        "satellite_year_strictly_before_photo":RASTER_YEAR < ANALYSIS_YEAR,
        "source_output_csv_sha256":sha256(filename),
        "source_join_sha256":sha256(joined),
        "original_PCNl_raster_not_VIIRS_500m":True,
        "sampling":"hash-sorted pre-outcome cap per island; exploratory site pilot",
        "no_insect_or_pollination_or_plant_population_effect_estimated":True,
        "no_ocean_or_inland_lake_dark_pixels_admitted":True,
        "dark_zero_PCNL_is_not_verified_no_artificial_ground_light":True,
    }
    (out / "dark_refugia_original_real_source_receipt.json").write_text(
        json.dumps(report, sort_keys=True, indent=2) + "\n")
    print("ALAN_DARK_REFUGIA_REAL_SOURCE", json.dumps(report, sort_keys=True), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--joined", required=True, type=Path)
    p.add_argument("--joined-receipt", required=True, type=Path)
    p.add_argument("--coast", required=True, type=Path)
    p.add_argument("--raster", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--radius", type=float, default=5.)
    p.add_argument("--max-per-island", type=int, default=8)
    a = p.parse_args()
    run(a.joined, a.joined_receipt, a.coast, a.raster, a.out, a.radius,
        a.max_per_island)


if __name__ == "__main__":
    main()
