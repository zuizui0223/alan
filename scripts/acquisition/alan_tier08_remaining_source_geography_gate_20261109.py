#!/usr/bin/env python3
"""Audit exact prior-photo-ID deduplicated tier08 remainder BEFORE another RGB wave.

This is source-metadata geography support only. It never downloads RGB, calculates
floral ROI9, or treats missing GloBI records as absence. Do not use the 30k
remaining opportunities as a power guarantee for the island-cluster CI.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.acquisition import alan_tier08_24576_balanced_area_light_photo_source_20261109 as source

PRECOLOUR_MANIFEST_SHA = "38d089278cfafd51d88ef9936e4882b9d115b2f8d0eefd0a9f4e4f2d0e4da33c"
EXPECTED_EXTRA_SOURCE = 55038
EXPECTED_NEW_MEASURED = 24576
EXPECTED_REMAINING = 30462


def require(value, message):
    if not bool(value):
        raise ValueError("ALAN_TIER08_REMAINDER_AUDIT: " + message)


def qualified_same_year_spatial_support(df):
    """Potential source-only pairs, not genuine measurement-qualified ROI sites."""
    if df.empty:
        return {"source_sites_with_5_photos_3_taxa": 0,
                "source_site_eligible_islands": 0,
                "island_years_with_two_sites_sameyear_light_gap_0p05": 0,
                "islands_with_two_sites_sameyear_light_gap_0p05": 0,
                "paired_eligible_source_sites": 0}
    require(df.photo_id.is_unique and df.observation_id.is_unique,
            "duplicate original photo or independent observation ID")
    sites = df.groupby(["land_id", "year", "target_id"], sort=False).agg(
        photos=("photo_id", "size"),
        uuids=("observation_id", "nunique"),
        accepted_taxa=("wcvp_accepted_plant_name_id", "nunique"),
        light=("log_light", "first"),
    ).reset_index()
    valid = sites.loc[(sites.photos >= 5) & (sites.uuids >= 5) &
                      (sites.accepted_taxa >= 3)].copy()
    if valid.empty:
        return {"source_sites_with_5_photos_3_taxa": 0,
                "source_site_eligible_islands": 0,
                "island_years_with_two_sites_sameyear_light_gap_0p05": 0,
                "islands_with_two_sites_sameyear_light_gap_0p05": 0,
                "paired_eligible_source_sites": 0}
    valid["light"] = pd.to_numeric(valid.light, errors="raise")
    pair = valid.groupby(["land_id", "year"]).agg(
        nsites=("target_id", "nunique"),
        min_light=("light", "min"),
        max_light=("light", "max"))
    pair = pair.loc[(pair.nsites >= 2) &
                    ((pair.max_light - pair.min_light) >= 0.05)]
    return {
        "source_sites_with_5_photos_3_taxa": int(len(valid)),
        "source_site_eligible_islands": int(valid.land_id.nunique()),
        "island_years_with_two_sites_sameyear_light_gap_0p05": int(len(pair)),
        "islands_with_two_sites_sameyear_light_gap_0p05":
            int(pair.reset_index().land_id.nunique()),
        "paired_eligible_source_sites": int(pair.nsites.sum()),
    }


def summarize(novel, measured):
    require(novel.photo_id.is_unique and novel.observation_id.is_unique,
            "source-only distinct-photo-ID/UUID requirement")
    require(measured.photo_id.is_unique and measured.observation_id.is_unique,
            "measured queue is not distinct")
    require(measured.photo_id.isin(novel.photo_id).all(),
            "frozen original 24,576 photos are not subset of verified novel cohort")
    remain = novel.loc[~novel.photo_id.isin(measured.photo_id)].copy()
    require(remain.photo_id.is_unique and remain.observation_id.is_unique and
            not remain.photo_id.isin(measured.photo_id).any(),
            "remaining original photo IDs overlap previously measured")
    area = pd.to_numeric(remain.area_km2, errors="raise")
    group = pd.cut(area, bins=[0,20,100,1000,10000,float("inf")],
                   labels=["area_le20", "area20_100", "area100_1000",
                           "area1000_10000", "area_gt10000"])
    support = qualified_same_year_spatial_support(remain)
    prior_support = qualified_same_year_spatial_support(measured)
    stats = {
        "schema": "ALAN_tier08_after24576_original_photo_remaining_geography_support_v1",
        "source_only_not_real_RGB_or_floral_ROI9": True,
        "all_prior_24576_are_signed_RGB_attempts_not_all_numeric": True,
        "tier08_eligible_after_exact_photo_ID_source_dedup": int(len(novel)),
        "verified_prior_24576_new_original_photo_IDs": int(len(measured)),
        "tier08_unattempted_exact_photo_ID_candidates": int(len(remain)),
        "remaining_distinct_original_observation_UUIDs": int(remain.observation_id.nunique()),
        "remaining_original_GSHHG_islands": int(remain.land_id.nunique()),
        "remaining_original_island_site_years": int(remain.target_id.nunique()),
        "remaining_photo_attempts_by_absolute_island_area_km2":
            {str(k): int(v) for k, v in group.value_counts(sort=False).items()},
        "remaining_source_only_geography": support,
        "previous_24576_source_only_geography": prior_support,
        "remaining_islands_not_already_in_24576_source_queue":
            int(len(set(remain.land_id) - set(measured.land_id))),
        "remaining_photo_sites_not_already_in_24576_source_queue":
            int(len(set(remain.target_id) - set(measured.target_id))),
        "historical_numeric_ID_to_new_UUID_observation_independence_unverified": True,
        "single_photographed_species_not_local_moth_pollination_proof": True,
        "source_only_pair_counts_not_upper_bound_when_combined_with_existing_successes": True,
        "do_not_disproportionately_spend_next_RGB_wave_on_large_old_islands": True,
        "tier09_requires_a_separate_SHApinned_source_and_photo_ID_dedup": True
    }
    return remain, stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for n in ("master", "previous910", "canonical", "source06", "source07", "prior24576",
              "prior24576_manifest", "out"):
        p.add_argument("--"+n, required=True, type=Path)
    args = p.parse_args()
    require(source.sha(args.prior24576_manifest) == PRECOLOUR_MANIFEST_SHA,
            "changed prior 24576 queue source manifest SHA")
    novel, audit = source.load_verified(args.master, args.previous910,
                                        args.canonical, args.source06, args.source07)
    measured = pd.read_csv(args.prior24576, compression="gzip", dtype=str,
                           keep_default_na=False)
    require(len(novel) == EXPECTED_EXTRA_SOURCE and len(measured) == EXPECTED_NEW_MEASURED,
            "prior source frame/queue counts changed")
    remain, stats = summarize(novel, measured)
    require(len(remain) == EXPECTED_REMAINING,
            "remainder after prior photo ID dedup does not reconcile")
    stats["prior_source_exact_photo_ID_overlap_excluded"] = int(
        audit["already_measured_exact_original_photo_ids_across_source_arms"])
    stats["frozen_prior_source_manifest_sha256"] = PRECOLOUR_MANIFEST_SHA
    stats["tier08_original_master_SHA256"] = source.MASTER_SHA
    out = args.out
    require(not out.exists(), "Never overwrite a finished source-only audit")
    out.mkdir(parents=True)
    path = out/"tier08_remaining30462_photo_IDs_source_only.csv.gz"
    remain.sort_values("photo_id", kind="stable").to_csv(
        path, compression={"method":"gzip","mtime":0},
        index=False, lineterminator="\n")
    stats["remaining_roster_sha256"] = source.sha(path)
    (out/"tier08_remaining_geographic_support_audit.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print("ALAN_TIER08_REMAINDER_PHOTO_ID_DISJOINT_GEOGRAPHIC_SUPPORT",
          json.dumps(stats, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
