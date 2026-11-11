#!/usr/bin/env python3
"""Source-only audit of tier09 qualified moth-photo cells absent from signed VIIRS.

This does not sample nightlight, download images, or fit ALAN models. A missing
existing-grid VIIRS value means UNKNOWN, never dark. A prospective two-cell
pair is not an observed light gradient. Selection is before observing ROI9.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

SCHEMA = "ALAN_tier09_missing_existing_VIIRS_precolour_geographic_priority_v1"


def require(condition, why):
    if not bool(condition):
        raise ValueError("ALAN_TIER09_MISSING_VIIRS: " + why)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def unique_file(root, filename):
    found = list(Path(root).rglob(filename))
    require(len(found) == 1, f"expected exactly one {filename}, found {len(found)}")
    return found[0]


def truth(s, column):
    x = s.astype(str).str.lower().map({"true": True, "false": False})
    require(x.notna().all(), f"missing or malformed {column}")
    return x.astype(bool)


def audit(source, original, source_receipt, original_receipt):
    sr = json.loads(Path(source_receipt).read_text())
    fr = json.loads(Path(original_receipt).read_text())
    require(sr.get("schema") == "ALAN_tier09_moth_sameyear_VIIRS_existingcell_source_gate_v1",
            "not the frozen tier09 same-year source receipt")
    require(fr.get("schema") == "ALAN_true112428_original_photo_RGB_GloBI_moth_guild_stratified_v1",
            "not the verified full 112428 ROI9 model")
    require(sha(source) == sr["candidate_status_SHA256"],
            "tier09 original photo roster SHA mismatch")
    geo = fr["deduplicated_geography_files"]
    require(sha(original) == geo["photo_id_SHA256"],
            "full model original photo-ID roster SHA mismatch")
    d = pd.read_csv(source, compression="gzip",
                    dtype={"photo_id": str, "year": str, "land_id": str,
                           "cell_id": str, "species": str}, low_memory=False)
    p = pd.read_csv(original, compression="gzip",
                    dtype={"photo_id": str, "year": str, "land_id": str,
                           "target_id": str}, low_memory=False)
    require(len(d) == 78214 and d.photo_id.is_unique,
            "tier09 photo counts or photo ID uniqueness drift")
    require(len(p) == 112428 and p.photo_id.is_unique,
            "measured union counts or photo ID uniqueness drift")
    require(not set(d.photo_id).intersection(p.photo_id),
            "new tier09 photos include already attempted exact photo ID")
    numeric = truth(p.real_numeric_ROI9, "full model numeric ROI9")
    require(int(numeric.sum()) == 40914,
            "previous genuine numeric original RGB count drift")
    d["viirs_eligible"] = truth(d.viirs_eligible, "existing cell VIIRS")
    require(pd.to_numeric(d.area_km2, errors="coerce").notna().all() and
            pd.to_numeric(d.area_km2, errors="coerce").gt(0).all(),
            "unverified positive island area")
    key = ["land_id", "year", "cell_id"]
    site = d.groupby(key, observed=True, sort=True).agg(
        photo_ids=("photo_id", "nunique"),
        taxa=("species", "nunique"),
        area_km2=("area_km2", "first"),
        viirs_valid=("viirs_eligible", "all"),
    ).reset_index()
    site = site.loc[site.photo_ids.ge(5) & site.taxa.ge(3)].copy()
    require(len(site) == 3174 and int(site.viirs_valid.sum()) == 2482,
            "signed 3174 qualified / 2482 VIIRS matched site-years drift")
    site["target_id"] = (site.year.astype(str) + "|" +
                         site.cell_id.astype(str) + "|" +
                         site.land_id.astype(str))
    require(site.target_id.is_unique, "duplicate target GSHHG same-year cell")
    numeric_targets = set(p.loc[numeric, "target_id"])
    measured_islands = set(p.loc[numeric, "land_id"])
    measured_iy = set(zip(p.loc[numeric, "land_id"], p.loc[numeric, "year"]))
    site["existing_true_numeric_ROI9"] = site.target_id.isin(numeric_targets)
    site["new_numeric_cell"] = ~site.existing_true_numeric_ROI9
    site["missing_existing_VIIRS"] = ~site.viirs_valid
    iy = site.groupby(["land_id", "year"], observed=True).agg(
        sites=("cell_id", "nunique"),
        valid_viirs_sites=("viirs_valid", "sum"),
        missing_viirs_sites=("missing_existing_VIIRS", "sum"),
    ).reset_index()
    # This flag only says two precolour source cells could be resampled:
    # it does NOT assert distinct VIIRS values, a gradient, or effective moth pollination.
    iy["possible_two_cell_recovery"] = (
        iy.sites.ge(2) &
        ((iy.missing_viirs_sites.ge(2)) |
         (iy.missing_viirs_sites.ge(1) & iy.valid_viirs_sites.ge(1)))
    )
    site = site.merge(iy, on=["land_id", "year"], how="left",
                      validate="many_to_one")
    site["novel_numeric_island"] = ~site.land_id.isin(measured_islands)
    site["novel_numeric_island_year"] = [
        (land, year) not in measured_iy for land, year in
        zip(site.land_id, site.year)]
    missing = site.loc[site.missing_existing_VIIRS].copy()
    missing["area_band"] = pd.cut(
        pd.to_numeric(missing.area_km2), bins=[0, 20, 100, float("inf")],
        labels=["le20km2", "20to100km2", "gt100km2"], include_lowest=True)
    missing = missing.sort_values(
        ["area_km2", "land_id", "year", "cell_id"], kind="mergesort")
    band = {}
    for name in ["le20km2", "20to100km2", "gt100km2"]:
        sub = missing.loc[missing.area_band.eq(name)]
        band[name] = {
            "missing_qualified_cells": int(len(sub)),
            "distinct_islands": int(sub.land_id.nunique()),
            "possible_two_cell_recovery_cells": int(sub.possible_two_cell_recovery.sum()),
            "not_yet_numeric_ROI9_cells": int(sub.new_numeric_cell.sum()),
            "unmeasured_numeric_islands": int(sub.loc[sub.novel_numeric_island, "land_id"].nunique()),
            "unmeasured_numeric_island_years": int(sub.loc[sub.novel_numeric_island_year, ["land_id", "year"]].drop_duplicates().shape[0]),
        }
    summary = {
        "schema": SCHEMA,
        "signed_tier09_candidate_status_SHA256": sha(source),
        "signed_complete112428_photo_status_SHA256": sha(original),
        "source_only_no_new_VIIRS_or_original_RGB_or_ROI9": True,
        "photo_id_overlap": 0,
        "qualified_source_cells": int(len(site)),
        "qualified_existingVIIRS_cells": int(site.viirs_valid.sum()),
        "qualified_missing_or_ineligibleVIIRS_cells": int(len(missing)),
        "missing_cells_by_area": band,
        "missing_distinct_islands": int(missing.land_id.nunique()),
        "missing_possible_two_cell_island_years": int(
            iy.loc[iy.possible_two_cell_recovery, ["land_id", "year"]]
            .merge(missing[["land_id", "year"]].drop_duplicates(),
                   on=["land_id", "year"]).drop_duplicates().shape[0]),
        "unmeasured_numeric_islands_in_missing_frame": int(
            missing.loc[missing.novel_numeric_island, "land_id"].nunique()),
        "all_missing_values_remain_unknown_not_zero": True,
        "not_a_verified_light_gradient_or_causal_ALAN_effect": True,
        "not_verified_moth_obligate_or_nighttime_pollinated": True,
        "no_historical_numeric_to_UUID_crosswalk": True,
        "bulk_RGB_authorized": False,
    }
    return summary, missing


def run(tier09_root, full_root, out):
    tier09 = unique_file(tier09_root,
        "tier09_moth_original_photo_existingVIIRS_cell_status.csv.gz")
    tier09_receipt = unique_file(tier09_root,
        "tier09_moth_VIIRS_existingcell_gate.json")
    photo = unique_file(full_root,
        "ALAN_112428_original_photoID_sitecell_status.csv.gz")
    full_receipt = unique_file(full_root,
        "ALAN_112428_true_RGB_GloBI_moth_strata.json")
    out = Path(out)
    require(not out.exists(), "refuse overwrite")
    result, missing = audit(tier09, photo, tier09_receipt, full_receipt)
    out.mkdir(parents=True)
    csv = out / "tier09_missing_VIIRS_qualified_cells_source_priority.csv.gz"
    missing.to_csv(csv, index=False, compression="gzip")
    result["precolour_qualified_missing_cells_SHA256"] = sha(csv)
    (out / "tier09_missing_VIIRS_priority_receipt.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n")
    print("ALAN_TIER09_MISSING_VIIRS_TRUE_SOURCE_PRIORITY", json.dumps(result, sort_keys=True),
          flush=True)
    return result


if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--tier09-root", type=Path, required=True)
    a.add_argument("--full-root", type=Path, required=True)
    a.add_argument("--out", type=Path, required=True)
    args = a.parse_args()
    run(args.tier09_root, args.full_root, args.out)
