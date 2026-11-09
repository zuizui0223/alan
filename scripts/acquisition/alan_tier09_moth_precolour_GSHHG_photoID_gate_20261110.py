#!/usr/bin/env python3
"""SHA-gated tier09 GloBI moth-positive source-only island/photo overlap audit.

This is NOT RGB inference or a same-year VIIRS comparison. It produces a
photo-ID-deduplicated, accepted-taxon GSHHG island site/year candidate frame.
Source cohort: signed 16 independent tier09 group CSVs from the completed
93,209,351 original photo-ID extraction. No GloBI negative/obligate moth claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import alan_moth_species_historical_pcnl_inventory_20261109 as original
from scripts.analysis import alan_isolation_area_light_20261108 as area
from scripts.acquisition import alan_full93m_extract_GloBI_moth_positive_original_photoIDs_20261109 as extraction
from scripts.acquisition import alan_tier08_24576_balanced_area_light_photo_source_20261109 as tier08

YEARS = ("2019", "2021", "2023", "2024", "2025")
GROUPS = 16
TIER = 9
EXPECTED_SOURCE_TIER09_MOTH = 1760604
EXPECTED_SOURCE_TIER09_FIVEYEAR = 1069431
EXPECTED_PRIOR_PHOTOS = 112428
KNOWN_INPUTS = {
    "canonical": ("f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c", "old"),
    "source06": ("e64f77d91ab9c19db7e38861baf96753395de5f1189f2599a1ac5cf1f55d8acd", "tier06"),
    "source07": ("7198455819f1ac260a6cbeacd9f4536c30d766595daa64dc5205c75e9d5a1c62", "tier07"),
    "previous910": (tier08.PREVIOUS910_SHA, "tier08_910"),
    "prior24576": (tier08.PRECOLOUR_MANIFEST_SHA if hasattr(tier08, "PRECOLOUR_MANIFEST_SHA") else "38d089278cfafd51d88ef9936e4882b9d115b2f8d0eefd0a9f4e4da33c", "tier08_24576"),
}
SCHEMA = "ALAN_tier09_moth_original_photoIDs_GSHHG_precolour_gate_v1"


def check(yes, why):
    if not bool(yes):
        raise ValueError("ALAN_TIER09_SOURCE_GATE: " + why)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_gz(d, path):
    d.to_csv(path, index=False, lineterminator="\n",
             compression={"method": "gzip", "mtime": 0})


def prior_ids(paths, out):
    out = Path(out)
    check(not out.exists(), "Refuse overwrite of prior exact-photo source")
    known = set()
    input_sha = {}
    for name, p in paths.items():
        if name == "prior24576_manifest":
            continue
        p = Path(p)
        expected = KNOWN_INPUTS[name][0]
        if name == "prior24576":
            # The source queue has its own signed manifest. Its compressed
            # CSV bytes are hashed and bound by this new receipt.
            check(sha(paths["prior24576_manifest"]) == expected,
                  "Frozen 24576 queue manifest source changed")
        else:
            check(sha(p) == expected, "Prior original source SHA drift: " + name)
        ids = pd.read_csv(p, compression="gzip", usecols=["photo_id"],
                          dtype=str, keep_default_na=False).photo_id
        check(ids.str.fullmatch(r"\d+").all(), "Non-numeric photo IDs: " + name)
        input_sha[name] = sha(p)
        known.update(ids)
    check(len(known) == EXPECTED_PRIOR_PHOTOS,
          "Actual photo-ID union must reconcile to 112428 original images")
    out.mkdir(parents=True)
    dest = out / "already_attempted_exact_photo_ids.csv.gz"
    write_gz(pd.DataFrame({"photo_id": sorted(known, key=int)}), dest)
    receipt = {
        "schema": SCHEMA, "purpose": "prior_attempted_photo_IDs_only",
        "known_exact_photo_ids": len(known),
        "known_photo_csv_sha256": sha(dest),
        "input_sha256": input_sha,
        "old_numeric_observation_ID_to_new_UUID_crosswalk_not_claimed": True,
        "RGB_colour_not_inferred": True,
    }
    (out / "known_photo_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    return receipt


def unique_precolour_candidates(raw, known, accepted, islands):
    """Retain strict-source candidate IDs; never infer light or native status."""
    check(raw.photo_id.is_unique and raw.photo_id.str.fullmatch(r"\d+").all(),
          "Repeated or malformed tier09 photo IDs")
    check(raw.tier.eq(str(TIER)).all(), "Non-tier09 original photo-ID row")
    check(raw.observation_uuid.ne("").all(), "Missing true original UUID")
    x = raw.loc[raw.observed_on.str.slice(0, 4).isin(YEARS) &
                raw.valid_observed_date.eq("true") &
                raw.valid_geopoint.eq("true") &
                raw.frozen_loader_supported.eq("true")].copy()
    x = x.loc[~x.photo_id.isin(known)].copy()
    x["year"] = x.observed_on.str.slice(0, 4)
    x["species"] = x.inat_taxon_id.map(accepted)
    x = x.loc[x.species.notna()].copy()
    x["latitude"] = pd.to_numeric(x.latitude, errors="coerce")
    x["longitude"] = pd.to_numeric(x.longitude, errors="coerce")
    x = x.loc[x.latitude.between(-90, 90) &
              x.longitude.between(-180, 180)].copy()
    if x.empty:
        return x, {"ambiguous": 0, "valid_island_photos": 0}
    # Uses the identical SHA-pinned GSHHG 2.3.7 L1 geometry as existing study.
    import geopandas as gpd
    points = gpd.GeoDataFrame(
        x, geometry=gpd.points_from_xy(x.longitude, x.latitude),
        crs="EPSG:4326")
    m = gpd.sjoin(points, islands, how="inner", predicate="within")
    repeated = m.photo_id.duplicated(keep=False)
    ambiguous = int(m.loc[repeated, "photo_id"].nunique())
    m = pd.DataFrame(m.loc[~repeated].drop(
        columns=["geometry", "index_right"]))
    check(m.photo_id.is_unique, "Island overlap ambiguity")
    if m.empty:
        return m, {"ambiguous": ambiguous, "valid_island_photos": 0}
    # One photograph per true UUID; deterministic, not colour-chosen.
    m = m.sort_values("photo_id", key=lambda s: s.map(int),
                      kind="stable").drop_duplicates("observation_uuid")
    m["latbin"] = np.floor((m.latitude + 90) * 20).astype(int)
    m["lonbin"] = np.floor((m.longitude + 180) * 20).astype(int)
    check(m.area_km2.gt(0).all(), "GSHHG source island area invalid")
    for flag in ("moth_database_pollinates", "walton2020_nocturnal_pollen_carriage",
                 "butterfly_record_positive"):
        check(m[flag].isin(["true", "false"]).all(), "Bad source flag " + flag)
    keep = ["photo_id", "observation_uuid", "inat_taxon_id", "species",
            "image_extension", "photo_license", "year", "land_id", "area_km2",
            "latbin", "lonbin", "moth_database_pollinates",
            "walton2020_nocturnal_pollen_carriage", "butterfly_record_positive"]
    return m[keep], {"ambiguous": ambiguous, "valid_island_photos": len(m)}


def run_group(group, source_csv, source_receipt, known_csv, known_receipt,
              coast, wcvp, accepted_names, out):
    check(0 <= group < GROUPS and not Path(out).exists(), "Group range or output exists")
    sr = json.loads(Path(source_receipt).read_text())
    check(sr.get("schema") == extraction.SCHEMA and sr.get("source_group") == group,
          "Wrong signed 93m group receipt")
    row = {int(v["tier"]): v for v in sr["source_tier_rows"]}[TIER]
    check(sha(source_csv) == row["moth_original_photo_file_SHA256"],
          "Extracted original tier09 photo roster SHA drift")
    prior = json.loads(Path(known_receipt).read_text())
    check(prior["schema"] == SCHEMA and
          prior["known_exact_photo_ids"] == EXPECTED_PRIOR_PHOTOS and
          sha(known_csv) == prior["known_photo_csv_sha256"],
          "Prior attempted photo roster not source verified")
    raw = pd.read_csv(source_csv, dtype=str, compression="gzip",
                      keep_default_na=False)
    check(len(raw) == row["moth_positive_original_photo_IDs"] and
          raw.photo_id.is_unique, "Tier09 source photo count or uniqueness drift")
    known = set(pd.read_csv(known_csv, compression="gzip",
                            dtype=str).photo_id)
    accepted = original.accepted_taxa(wcvp, accepted_names)
    land = area.load_signed_land(coast)
    islands = land.loc[land.land_class.eq("island"),
                       ["land_id", "area_km2", "geometry"]].to_crs("EPSG:4326")
    result, audit = unique_precolour_candidates(raw, known, accepted, islands)
    out = Path(out)
    out.mkdir(parents=True)
    roster = out / f"tier09_moth_island_group_{group:02d}.csv.gz"
    write_gz(result, roster)
    count_fiveyear = int(raw.observed_on.str.slice(0, 4).isin(YEARS).sum())
    receipt = {
        "schema": SCHEMA, "group": group,
        "source_tier": TIER,
        "source_total_moth_photos": len(raw),
        "source_fiveyear_moth_photos": count_fiveyear,
        "exact_prior_photo_ID_overlap": int(raw.photo_id.isin(known).sum()),
        "one_photo_per_UUID_island_candidates": len(result),
        "GSHHG_island_count": int(result.land_id.nunique()) if len(result) else 0,
        "candidate_csv_sha256": sha(roster),
        "original_extraction_SHA256": sha(source_csv),
        "prior_exact_photo_ids_SHA256": sha(known_csv),
        "input_93m_group_receipt_SHA256": sha(source_receipt),
        "GSHHG_SHA_pinned_via_original_loader": True,
        "year_matched_VIIRS_NOT_YET_JOINED": True,
        "RGB_measured": 0,
        "moth_source_positive_not_effective_nocturnal_dependence": True,
        "no_missing_GloBI_edge_interpreted_as_moth_absence": True,
        **audit,
    }
    (out / f"tier09_moth_island_group_{group:02d}.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print("ALAN_TIER09_SIGNED_MOTH_PRECOLOUR_GROUP", json.dumps(receipt, sort_keys=True),
          flush=True)


def collect(folder, out):
    out, folder = Path(out), Path(folder)
    check(not out.exists(), "Refuse overwrite of finished tier09 audit")
    all_ = []
    receipts = []
    for i in range(GROUPS):
        p = folder / f"tier09_moth_island_group_{i:02d}.csv.gz"
        r = json.loads((folder / f"tier09_moth_island_group_{i:02d}.json").read_text())
        check(r["schema"] == SCHEMA and r["group"] == i and
              sha(p) == r["candidate_csv_sha256"] and
              r["year_matched_VIIRS_NOT_YET_JOINED"], "Missing signed tier09 group")
        all_.append(pd.read_csv(p, dtype=str, compression="gzip",
                                keep_default_na=False))
        receipts.append(r)
    check(sum(r["source_total_moth_photos"] for r in receipts)
          == EXPECTED_SOURCE_TIER09_MOTH, "16 group tier09 source total drift")
    check(sum(r["source_fiveyear_moth_photos"] for r in receipts)
          == EXPECTED_SOURCE_TIER09_FIVEYEAR, "16 group fiveyear total drift")
    d = pd.concat(all_, ignore_index=True)
    check(d.photo_id.is_unique and d.observation_uuid.is_unique,
          "Original photo/UUID duplicated across source groups")
    if not d.empty:
        site = d.groupby(["land_id", "year", "latbin", "lonbin"],
                         observed=True).agg(
            photos=("photo_id", "size"),
            taxa=("species", "nunique"),
            area_km2=("area_km2", "first")).reset_index()
        valid = site.loc[(site.photos >= 5) & (site.taxa >= 3)].copy()
        pairs = valid.groupby(["land_id", "year"]).size()
        pair_islands = set(pairs[pairs >= 2].index.get_level_values(0))
        area = pd.to_numeric(d.area_km2)
        bins = pd.cut(area, [0, 20, 100, 1000, 10000, float("inf")],
                      labels=["le20", "20_100", "100_1000",
                              "1000_10000", "gt10000"])
        counts = {str(k): int(v) for k, v in bins.value_counts(sort=False).items()}
    else:
        valid = pd.DataFrame()
        pair_islands = set()
        counts = {}
    out.mkdir(parents=True)
    roster = out / "tier09_moth_GSHHG_exact_photoID_novel_precolour.csv.gz"
    write_gz(d, roster)
    result = {
        "schema": SCHEMA, "source_only": True,
        "all_source_tier09_moth_positive": EXPECTED_SOURCE_TIER09_MOTH,
        "five_target_year_tier09_moth_positive": EXPECTED_SOURCE_TIER09_FIVEYEAR,
        "already_attempted_exact_photo_ids_union": EXPECTED_PRIOR_PHOTOS,
        "tier09_photo_ids_overlapping_prior":
            sum(r["exact_prior_photo_ID_overlap"] for r in receipts),
        "unique_original_UUID_island_candidates_after_photo_dedup": len(d),
        "source_islands": int(d.land_id.nunique()) if len(d) else 0,
        "source_site_years": int(d[["land_id", "year", "latbin", "lonbin"]].drop_duplicates().shape[0]) if len(d) else 0,
        "source_5photo_3taxa_site_years": len(valid),
        "source_islands_with_two_sameyear_qualified_sites_NO_VIIRS_YET": len(pair_islands),
        "island_area_photo_candidates": counts,
        "GloBI_positive_not_local_moth_dependence": True,
        "sameyear_VIIRS_light_gradient_NOT_VERIFIED": True,
        "RGB_downloads": 0, "numeric_flower_ROI9": 0,
        "no_new_RGB_wave_automatically_dispatched": True,
        "signed_source_candidate_roster_sha256": sha(roster),
        "source_group_receipts": {str(x["group"]): x["candidate_csv_sha256"]
                                   for x in receipts},
    }
    (out / "tier09_moth_precolour_geography_gate.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n")
    print("ALAN_TIER09_GLOBAL_GSHHG_MOTH_SOURCE_ONLY_GATE",
          json.dumps(result, sort_keys=True), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="mode", required=True)
    q = sub.add_parser("prior")
    for n in KNOWN_INPUTS:
        q.add_argument("--" + n, required=True, type=Path)
    q.add_argument("--prior24576_manifest", required=True, type=Path)
    q.add_argument("--out", required=True, type=Path)
    q = sub.add_parser("group")
    for n in ("source_csv", "source_receipt", "known_csv", "known_receipt",
              "coast", "wcvp", "accepted_names", "out"):
        q.add_argument("--" + n, required=True, type=Path)
    q.add_argument("--group", required=True, type=int)
    q = sub.add_parser("collect")
    q.add_argument("--folder", required=True, type=Path)
    q.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    if a.mode == "prior":
        prior_ids({k: getattr(a, k) for k in [*KNOWN_INPUTS, "prior24576_manifest"]}, a.out)
    elif a.mode == "group":
        run_group(a.group, a.source_csv, a.source_receipt, a.known_csv,
                  a.known_receipt, a.coast, a.wcvp, a.accepted_names, a.out)
    else:
        collect(a.folder, a.out)


if __name__ == "__main__":
    main()
