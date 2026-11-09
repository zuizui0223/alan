#!/usr/bin/env python3
"""Source-only island census for published moth-visited plant taxa + historic ALAN.

Inventory consumes the 14 signed original-photo-ID rosters PER source group
from ALAN PR #22. No RGB, flower ROI, species absence or local effective
moth pollination is inferred. The labels are GloBI-positive visitation evidence,
NOT moth dependence. One biological observation UUID is counted once despite
multiple photos. Spatial assignment uses genuine source geopoints and the
SHA-verified GSHHG landmass polygons; ambiguous/sea points are omitted.

Long-term exposure uses a SINGLE cross-sensor harmonized PCNL annual product
(1992--2024, Zenodo v3), with years strictly PRIOR to the photographic year.
Native DMSP digital numbers and VIIRS radiance are never added together.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sqlite3
from collections import Counter
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

YEARS = {2019, 2021, 2023, 2024, 2025}
FIELDS = {"photo_id", "observation_uuid", "inat_taxon_id", "observed_on",
          "latitude", "longitude", "valid_observed_date", "valid_geopoint",
          "moth_database_pollinates", "walton2020_nocturnal_pollen_carriage",
          "butterfly_record_positive"}
TABLE = """CREATE TABLE IF NOT EXISTS taxa (
    island TEXT NOT NULL, latbin INTEGER NOT NULL, lonbin INTEGER NOT NULL,
    year INTEGER NOT NULL, species TEXT NOT NULL,
    n_photos INTEGER NOT NULL, n_observations INTEGER NOT NULL,
    n_pollinates_photos INTEGER NOT NULL, n_walton_photos INTEGER NOT NULL,
    n_butterfly_copositive_photos INTEGER NOT NULL,
    representative_lon REAL NOT NULL, representative_lat REAL NOT NULL,
    PRIMARY KEY(island,latbin,lonbin,year,species)
)"""
UPSERT = """INSERT INTO taxa VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
ON CONFLICT(island,latbin,lonbin,year,species) DO UPDATE SET
n_photos=n_photos+excluded.n_photos,
n_observations=n_observations+excluded.n_observations,
n_pollinates_photos=n_pollinates_photos+excluded.n_pollinates_photos,
n_walton_photos=n_walton_photos+excluded.n_walton_photos,
n_butterfly_copositive_photos=n_butterfly_copositive_photos+excluded.n_butterfly_copositive_photos"""


def require(condition, message):
    if not condition:
        raise ValueError("ALAN_MOTH_LONGTERM: " + message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def accepted_taxa(wcvp_csv, accepted_csv):
    """Only uniquely WCVP-accepted original iNat names are admitted."""
    from scripts.analysis import alan_full93m_globi_moth_original_photo_census_20261109 as src
    require(sha256(wcvp_csv) == src.NAMES_SHA, "original exact WCVP mapping changed")
    require(sha256(accepted_csv) == src.ACCEPTED_SHA, "original accepted-name catalog changed")
    m = pd.read_csv(wcvp_csv, dtype=str, compression="gzip", keep_default_na=False,
                    usecols=["inat_taxon_id", "wcvp_match_state", "wcvp_accepted_plant_name_id"])
    a = pd.read_csv(accepted_csv, dtype=str, compression="gzip", keep_default_na=False,
                    usecols=["wcvp_accepted_plant_name_id", "wcvp_accepted_name"])
    require(m.inat_taxon_id.is_unique, "iNat taxon IDs not unique")
    a = a.drop_duplicates()
    require(a.wcvp_accepted_plant_name_id.is_unique, "WCVP accepted IDs ambiguous")
    m = m.loc[m.wcvp_match_state.eq("unique_exact_name_to_accepted_id")]
    m = m.merge(a, on="wcvp_accepted_plant_name_id", how="left", validate="many_to_one")
    m = m.loc[m.wcvp_accepted_name.notna() & m.wcvp_accepted_name.ne("")]
    return m.set_index("inat_taxon_id").wcvp_accepted_name.to_dict()


def load_islands(coast_dir):
    from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land
    land = load_signed_land(coast_dir)
    land = land.loc[land.land_class.eq("island"), ["land_id", "geometry"]].copy()
    require(land.land_id.is_unique, "GSHHG island land_id not unique")
    require(land.crs is not None, "GSHHG geometry missing CRS")
    return land.to_crs("EPSG:4326").rename(columns={"land_id": "island"})


def join_geography(chunk, islands):
    """Assign exact observed geopoints, reject multipart overlaps and sea points."""
    points = gpd.GeoDataFrame(chunk.copy(),
        geometry=gpd.points_from_xy(chunk.longitude, chunk.latitude), crs="EPSG:4326")
    matched = gpd.sjoin(points, islands, how="inner", predicate="within")
    # A photo on multiple polygons is ambiguous, not multiple independent records.
    ambiguous = matched.photo_id.duplicated(keep=False)
    ambiguous_count = int(matched.loc[ambiguous, "photo_id"].nunique())
    matched = matched.loc[~ambiguous].drop(columns=["geometry", "index_right"])
    require(matched.photo_id.is_unique, "Ambiguous photographic island membership")
    return pd.DataFrame(matched), ambiguous_count


def normalized_frame(frame, taxon_lookup):
    require(FIELDS.issubset(frame.columns), "Original moth-positive roster columns missing")
    require(frame.photo_id.is_unique, "Duplicate original photo ID within file chunk")
    # These files ONLY contain positively documented moth-associated plant taxa.
    x = frame.loc[frame.valid_observed_date.eq("true") & frame.valid_geopoint.eq("true")].copy()
    x["year"] = pd.to_numeric(x.observed_on.str.slice(0, 4), errors="coerce")
    x = x.loc[x.year.isin(YEARS)].copy()
    x["species"] = x.inat_taxon_id.map(taxon_lookup)
    x = x.loc[x.species.notna()].copy()
    x["latitude"] = pd.to_numeric(x.latitude, errors="coerce")
    x["longitude"] = pd.to_numeric(x.longitude, errors="coerce")
    x = x.loc[x.latitude.between(-90, 90) & x.longitude.between(-180, 180)].copy()
    x["year"] = x.year.astype(int)
    return x


def aggregate_source_chunk(frame, islands, seen_uuids):
    """Return grouped SQL tuples; only first UUID in the disjoint group counts."""
    matched, ambiguous = join_geography(frame, islands)
    if matched.empty:
        return [], 0, ambiguous
    matched["latbin"] = np.floor((matched.latitude + 90) * 20).astype(int)
    matched["lonbin"] = np.floor((matched.longitude + 180) * 20).astype(int)
    first = ~matched.observation_uuid.duplicated(keep="first")
    first &= ~matched.observation_uuid.isin(seen_uuids)
    seen_uuids.update(matched.observation_uuid.tolist())
    matched["one_observation"] = first.astype(int)
    for flag in ["moth_database_pollinates", "walton2020_nocturnal_pollen_carriage",
                 "butterfly_record_positive"]:
        require(matched[flag].isin(["true", "false"]).all(), "GloBI evidence flag invalid: " + flag)
        matched[flag] = matched[flag].eq("true").astype(int)
    group = matched.groupby(["island", "latbin", "lonbin", "year", "species"],
                            sort=False, observed=True).agg(
        n_photos=("photo_id", "size"),
        n_observations=("one_observation", "sum"),
        n_pollinates_photos=("moth_database_pollinates", "sum"),
        n_walton_photos=("walton2020_nocturnal_pollen_carriage", "sum"),
        n_butterfly_copositive_photos=("butterfly_record_positive", "sum"),
        representative_lon=("longitude", "first"),
        representative_lat=("latitude", "first")).reset_index()
    return list(group.itertuples(index=False, name=None)), len(matched), ambiguous


def run_inventory(args):
    files = sorted(Path(args.rosters).glob(f"moth-positive-original-tier-*-group-{args.group:02d}.csv.gz"))
    require(len(files) == 14, "Must have 14 signed source-tier rosters for one group")
    require(len({int(p.name.split("-tier-")[1][:2]) for p in files}) == 14,
            "Duplicate or missing original tier in group")
    require(not args.db.exists(), "Refusing to overwrite source group database")
    lookup = accepted_taxa(args.wcvp, args.accepted)
    islands = load_islands(args.coast)
    db = sqlite3.connect(args.db)
    db.execute(TABLE)
    seen_uuids = set()
    stats = Counter()
    input_hashes = {}
    for path in files:
        input_hashes[path.name] = sha256(path)
        for raw in pd.read_csv(path, compression="gzip", dtype=str, keep_default_na=False,
                               chunksize=args.chunksize):
            stats["all_original_moth_positive_roster_photo_ids"] += len(raw)
            work = normalized_frame(raw, lookup)
            stats["fiveyear_geopoint_WCVP_source_photos"] += len(work)
            if work.empty:
                continue
            rows, n_assigned, n_ambiguous = aggregate_source_chunk(work, islands, seen_uuids)
            stats["photos_contained_in_exact_single_GSHHG_island"] += n_assigned
            stats["ambiguous_landmass_photo_ids"] += n_ambiguous
            db.executemany(UPSERT, rows)
        db.commit()
    stats["genuine_unique_observation_uuids_included"] = len(seen_uuids)
    db.close()
    receipt = {"schema": "ALAN_GloBI_moth_positive_island_original_photo_metadata_inventory_v1",
               "group": args.group, "counts": dict(stats), "input_sha256": input_hashes,
               "original_RGB_downloads": 0, "photo_ROI_colour_inferred": 0,
               "moth_positive_is_NOT_local_moth_dependency": True,
               "missing_GloBI_records_are_NOT_moth_negative": True,
               "all_source_93m_denominator_NOT_included": True}
    args.db.with_suffix(".json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


def export_csv(db, sql, path):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        out = csv.writer(f)
        cur = db.execute(sql)
        out.writerow([column[0] for column in cur.description])
        for row in cur:
            out.writerow(row)


def run_collect(args):
    """Merge exactly 16 disjoint signed observation-hash groups."""
    databases = [Path(args.db_dir) / f"moth-island-group-{g:02d}.sqlite" for g in range(16)]
    require(all(p.is_file() and p.with_suffix(".json").is_file() for p in databases),
            "Missing one or more of 16 original group databases/receipts")
    require(not args.out_dir.exists(), "Refusing to overwrite global source aggregation")
    args.out_dir.mkdir(parents=True)
    merged = sqlite3.connect(args.out_dir / "moth_island_source.sqlite")
    merged.execute(TABLE)
    total_roster = 0
    for group, path in enumerate(databases):
        receipt = json.loads(path.with_suffix(".json").read_text())
        require(receipt["group"] == group and len(receipt["input_sha256"]) == 14,
                "Source group receipt is incomplete")
        total_roster += receipt["counts"]["all_original_moth_positive_roster_photo_ids"]
        with sqlite3.connect(path) as inp:
            cur = inp.execute("SELECT * FROM taxa")
            while True:
                part = cur.fetchmany(4000)
                if not part:
                    break
                merged.executemany(UPSERT, part)
            merged.commit()
    require(total_roster == 34_346_428, "Full signed 34,346,428 moth-positive IDs not conserved")
    export_csv(merged, """SELECT island,year,species,SUM(n_photos) AS photo_ids,
        SUM(n_observations) AS observation_uuids,
        SUM(n_pollinates_photos) AS pollinates_photo_ids,
        SUM(n_walton_photos) AS nocturnal_pollen_carriage_photo_ids
        FROM taxa GROUP BY island,year,species ORDER BY island,year,species""",
        args.out_dir / "moth_island_species_year.csv")
    export_csv(merged, """SELECT island,latbin,lonbin,year,
        SUM(n_photos) AS photo_ids,SUM(n_observations) AS observation_uuids,
        COUNT(DISTINCT species) AS moth_positive_species,
        MIN(representative_lon) AS lon,MIN(representative_lat) AS lat
        FROM taxa GROUP BY island,latbin,lonbin,year
        ORDER BY island,latbin,lonbin,year""",
        args.out_dir / "moth_island_site_year.csv")
    export_csv(merged, """SELECT island,year,COUNT(DISTINCT species) AS moth_positive_species,
        SUM(n_photos) AS photo_ids,SUM(n_observations) AS observation_uuids,
        COUNT(DISTINCT printf('%d:%d',latbin,lonbin)) AS source_sites
        FROM taxa GROUP BY island,year ORDER BY island,year""",
        args.out_dir / "moth_island_year.csv")
    merged.close()
    (args.out_dir / "source_receipt.json").write_text(json.dumps({
        "schema": "ALAN_GloBI_moth_positive_16group_geography_conserved_v1",
        "all_moth_positive_original_source_IDs": total_roster, "group_count": 16,
        "not_a_plant_flora_fraction": True, "missing_GloBI_not_nonmoth": True,
        "source_island_year_records_only": True}, indent=2) + "\n")


def annual_exposure(site_year, raster_directory, years=range(1992, 2025)):
    """PCNL annual raster sample; right-censored strictly before photo year."""
    import rasterio
    from rasterio.warp import transform
    x = site_year.copy()
    require({"year", "lat", "lon"}.issubset(x.columns), "Missing actual-site coordinates")
    require(x.year.isin(YEARS).all(), "Unexpected observation year")
    n = len(x)
    exposure = np.zeros(n, dtype=float)
    observed = np.zeros(n, dtype=np.int16)
    start, stop = min(years), max(years)
    require(start == 1992 and stop == 2024 and len(set(years)) == 33,
            "Annual harmonized PCNL coverage must be 1992 through 2024")
    files = [Path(raster_directory) / f"PCNL{yr}.tif" for yr in years]
    require(all(p.is_file() for p in files), "Missing annual PCNL harmonized raster")
    coords = list(zip(x.lon.astype(float), x.lat.astype(float)))
    for yr, file in zip(years, files):
        active = np.flatnonzero(x.year.to_numpy() > yr)
        if not len(active):
            continue
        with rasterio.open(file) as src:
            require(src.count == 1 and src.crs is not None, "PCNL annual raster invalid")
            c = [coords[i] for i in active]
            if src.crs.to_epsg() != 4326:
                xx, yy = transform("EPSG:4326", src.crs, [p[0] for p in c], [p[1] for p in c])
                c = list(zip(xx, yy))
            # Iterating in bounded chunks avoids materializing the whole raster.
            for start_idx in range(0, len(active), 16000):
                targets = active[start_idx:start_idx + 16000]
                sub = c[start_idx:start_idx + 16000]
                values = list(src.sample(sub, indexes=1, masked=True))
                numeric = np.asarray([float(v[0]) if not np.ma.is_masked(v[0])
                                       else np.nan for v in values], dtype=float)
                good = np.isfinite(numeric) & (numeric >= 0)
                exposure[targets[good]] += np.log1p(numeric[good])
                observed[targets[good]] += 1
    expected = x.year.to_numpy(dtype=int) - 1992
    complete = (observed == expected)
    x["pcnl_years_expected"] = expected
    x["pcnl_years_observed"] = observed
    x["pcnl_1992_to_previous_year_sum_log1p"] = np.where(complete, exposure, np.nan)
    x["pcnl_1992_to_previous_year_mean_log1p"] = np.where(
        complete, exposure / np.maximum(expected, 1), np.nan)
    x["historic_ntl_complete"] = complete
    return x


def run_exposure(args):
    require(not args.output.exists(), "Refusing to overwrite PCNL exposure output")
    x = pd.read_csv(args.site_year)
    require(x[["island", "latbin", "lonbin", "year"]].duplicated().sum() == 0,
            "Same site/year repeated; collect sources before joining")
    out = annual_exposure(x, args.pcnl_dir)
    out.to_csv(args.output, index=False)
    summary = out.groupby(["island", "year"], observed=True).agg(
        n_site_years=("year", "size"),
        observed_moth_positive_species_site_sum=("moth_positive_species", "sum"),
        original_photo_ids=("photo_ids", "sum"),
        observation_uuids=("observation_uuids", "sum"),
        complete_exposure_sites=("historic_ntl_complete", "sum"),
        median_pcnl_prior_mean_log1p=("pcnl_1992_to_previous_year_mean_log1p", "median")
    ).reset_index()
    summary.to_csv(args.output.with_name("moth_island_year_with_longterm_ntl.csv"), index=False)
    args.output.with_suffix(".json").write_text(json.dumps({
        "schema": "ALAN_PCNL_preobservation_historical_light_site_year_v1",
        "raster_product": "PCNL v3 (Li and Cao, Zenodo 17013414), 1992-2024",
        "years": [1992, 2024], "strictly_pre_photography_exposure": True,
        "complete_annual_support_required": True,
        "PCNL_not_raw_DMSP_DN_or_VIIRS_radiance": True,
        "not_illuminance_at_flower": True,
        "site_summed_species_counts_NOT_independent_unique_island_species": True,
        "moth_visitor_not_moth_dependency": True,
        "no_local_species_declines_or_plant_absences_inferred": True
    }, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("inventory")
    a.add_argument("--rosters", type=Path, required=True)
    a.add_argument("--group", type=int, required=True, choices=range(16))
    a.add_argument("--wcvp", type=Path, required=True)
    a.add_argument("--accepted", type=Path, required=True)
    a.add_argument("--coast", type=Path, required=True)
    a.add_argument("--db", type=Path, required=True)
    a.add_argument("--chunksize", type=int, default=50000)
    b = sub.add_parser("collect")
    b.add_argument("--db-dir", type=Path, required=True)
    b.add_argument("--out-dir", type=Path, required=True)
    c = sub.add_parser("exposure")
    c.add_argument("--site-year", type=Path, required=True)
    c.add_argument("--pcnl-dir", type=Path, required=True)
    c.add_argument("--output", type=Path, required=True)
    x = p.parse_args()
    if x.command == "inventory":
        require(x.chunksize > 0, "chunksize must be positive")
        run_inventory(x)
    elif x.command == "collect":
        run_collect(x)
    else:
        run_exposure(x)


if __name__ == "__main__":
    main()
