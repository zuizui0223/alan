#!/usr/bin/env python3
"""Conservative tier09 moth-photo support in *already sampled* same-year VIIRS cells.

Not a new global VIIRS raster sampling, not measured RGB, and not a causal fit.
Every exact source cell-year must have a signed original VIIRS observation;
other tier09 cells retain an explicit 'not yet sampled' status.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA = "ALAN_tier09_moth_sameyear_VIIRS_existingcell_source_gate_v1"
SOURCE_SCHEMA = "ALAN_tier09_moth_original_photoIDs_GSHHG_precolour_gate_v1"
SOURCE_N = 78214
YEARS = (2019, 2021, 2023, 2024, 2025)
SHA_2024 = "6e0a077d496e90848d05d043bffae57e4772519f4572a457221c5ef4da0fb813"
LIGHT_MIN_LOG_SPAN = 0.05

def require(test, reason):
    if not bool(test):
        raise ValueError("ALAN_TIER09_EXISTING_VIIRS: " + reason)

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(4 << 20), b""):
            h.update(part)
    return h.hexdigest()

def write_gz(d, path):
    d.to_csv(path, index=False, lineterminator="\n",
             compression={"method": "gzip", "mtime": 0})

def source_photos(path, receipt_path):
    receipt = json.loads(Path(receipt_path).read_text())
    require(receipt["schema"] == SOURCE_SCHEMA and
            receipt["unique_original_UUID_island_candidates_after_photo_dedup"] == SOURCE_N and
            receipt["sameyear_VIIRS_light_gradient_NOT_VERIFIED"] is True and
            receipt["RGB_downloads"] == 0, "Prior GSHHG source-only receipt drift")
    require(sha256(path) == receipt["signed_source_candidate_roster_sha256"],
            "Tier09 exact-photo roster SHA drift")
    use = ["photo_id", "observation_uuid", "species", "year", "land_id",
           "latbin", "lonbin", "area_km2", "moth_database_pollinates",
           "walton2020_nocturnal_pollen_carriage", "butterfly_record_positive"]
    d = pd.read_csv(path, compression="gzip", usecols=use, dtype=str, keep_default_na=False)
    require(len(d) == SOURCE_N and d.photo_id.is_unique and d.observation_uuid.is_unique,
            "Original source photo or UUID count/uniqueness drift")
    require(d.year.isin([str(y) for y in YEARS]).all(), "Unexpected source year")
    require(d.land_id.ne("").all() and d.species.ne("").all(), "Missing island or species")
    x = pd.to_numeric(d.lonbin, errors="raise").astype(int)
    y = pd.to_numeric(d.latbin, errors="raise").astype(int)
    require(x.between(0, 7199).all() and y.between(0, 3599).all(),
            "Source 0.05-degree cell out of bounds")
    d["cell_id"] = x.astype(str) + "_" + y.astype(str)
    d["target_id"] = d.year + "|" + d.cell_id + "|" + d.land_id
    require(pd.to_numeric(d.area_km2).gt(0).all(), "Invalid GSHHG island area")
    return d, receipt

def source_viirs(year, filename, receipt=None):
    year = int(year)
    require(year in YEARS, "Unapproved year")
    digest = sha256(filename)
    if year in (2019, 2021, 2023):
        r = json.loads(Path(receipt).read_text())
        require(r["source_year"] == year and r["values_SHA256"] == digest,
                "Signed original quarterly VIIRS values/year drift")
    elif year == 2024:
        require(digest == SHA_2024, "Signed 2024 VIIRS/Global2 table drift")
    else:
        r = json.loads(Path(receipt).read_text())
        require(r["VIIRS_year"] == r["source_photo_year"] == year and
                r["output_values_SHA256"] == digest and
                r["original_target_count"] == 73463,
                "Signed 2025 original VIIRS values/year drift")
    columns = ["target_id", "viirs_median", "viirs_valid_months"]
    v = pd.read_csv(filename, usecols=columns, dtype={"target_id": str}, low_memory=False)
    require(v.target_id.is_unique, "Original VIIRS target ID repeated")
    parts = v.target_id.str.split("|", regex=False, n=2, expand=True)
    require(parts.shape[1] == 3 and parts[0].eq(str(year)).all(),
            "VIIRS target year or GSHHG target syntax drift")
    v["year"], v["cell_id"], v["land_id"] = parts[0], parts[1], parts[2]
    v["viirs_median"] = pd.to_numeric(v.viirs_median, errors="coerce")
    v["viirs_valid_months"] = pd.to_numeric(v.viirs_valid_months, errors="coerce")
    require(v.viirs_valid_months.between(0, 12).fillna(False).all(),
            "Original valid-month counts malformed")
    # The 0.05-degree cell may contain more than one landmass.
    # Never use a mean light value or arbitrarily choose one of those landmasses.
    ambiguous = v.duplicated(["year", "cell_id"], keep=False)
    v = v.loc[~ambiguous].copy()
    return v[columns + ["year", "cell_id", "land_id"]], {
        "year": year, "viirs_original_values_SHA256": digest,
        "original_target_rows": int(len(ambiguous)),
        "ambiguous_landmass_cell_rows_excluded": int(ambiguous.sum()),
        "unambiguous_viirs_source_cells": int(len(v)),
    }

def score_support(photos, viirs):
    require(photos.target_id.is_unique is False or len(photos) > 0,
            "Empty original-photo source")
    require(viirs.target_id.is_unique, "VIIRS source duplicated")
    merged = photos.merge(viirs[["target_id", "viirs_median", "viirs_valid_months"]],
                          on="target_id", how="left", validate="many_to_one", indicator=True)
    merged["existingcell_source_hit"] = merged["_merge"].eq("both")
    merged["viirs_eligible"] = (merged.existingcell_source_hit &
        merged.viirs_valid_months.ge(3) & merged.viirs_median.notna())
    merged["log_light"] = np.where(merged.viirs_eligible,
        np.log1p(merged.viirs_median.clip(lower=0)), np.nan)
    return merged

def cohort_support(frame):
    site = frame.groupby(["land_id", "year", "cell_id"], observed=True).agg(
        photos=("photo_id", "size"), accepted_species=("species", "nunique"),
        area_km2=("area_km2", "first"), light=("log_light", "first"),
        viirs_eligible=("viirs_eligible", "all")).reset_index()
    qualified = site.loc[site.photos.ge(5) & site.accepted_species.ge(3)].copy()
    light_sites = qualified.loc[qualified.viirs_eligible & qualified.light.notna()].copy()
    if light_sites.empty:
        pairs = pd.DataFrame(columns=["land_id", "year", "sites", "span"])
        eligible = pairs.copy()
    else:
        pairs = light_sites.groupby(["land_id", "year"]).agg(
            sites=("cell_id", "nunique"),
            light_min=("light", "min"), light_max=("light", "max")).reset_index()
        pairs["span"] = pairs.light_max - pairs.light_min
        eligible = pairs.loc[pairs.sites.ge(2) & pairs.span.ge(LIGHT_MIN_LOG_SPAN)].copy()
    ids = set(eligible.land_id) if len(eligible) else set()
    areas = (light_sites.drop_duplicates("land_id").set_index("land_id").area_km2
             if len(light_sites) else pd.Series(dtype=float))
    light_islands = light_sites.land_id.nunique() if len(light_sites) else 0
    small = 0
    if ids:
        area = light_sites.drop_duplicates("land_id").set_index("land_id").area_km2.astype(float)
        small = int(area.loc[list(ids)].le(20).sum())
    return {
        "photos": int(len(frame)), "source_site_years": int(len(site)),
        "source_qualified_5photo_3taxa_sites": int(len(qualified)),
        "qualified_existing_valid_VIIRS_sites": int(len(light_sites)),
        "qualified_existing_valid_VIIRS_islands": int(light_islands),
        "supported_sameyear_island_year_pairs_with_loglight_span_ge_0p05": int(len(eligible)),
        "independent_supported_islands": int(len(ids)),
        "supported_small_islands_le20km2": small,
        "existing_VIIRS_absent_or_ineligible_qualified_sites": int(len(qualified) - len(light_sites)),
    }, eligible

def run(source, receipt, viirs_paths, viirs_receipts, output):
    output = Path(output)
    require(not output.exists(), "Refuse overwrite of signed audit")
    d, prior = source_photos(source, receipt)
    tables, signed = [], []
    for y in YEARS:
        v, meta = source_viirs(y, viirs_paths[y], viirs_receipts.get(y))
        tables.append(v)
        signed.append(meta)
    full = pd.concat(tables, ignore_index=True)
    require(full.target_id.is_unique, "VIIRS target duplicated across years")
    m = score_support(d, full)
    groups = {
        "GloBI_moth_record_positive": m,
        "GloBI_database_pollinates": m.loc[m.moth_database_pollinates.eq("true")],
        "Walton_nocturnal_pollen_carriage": m.loc[m.walton2020_nocturnal_pollen_carriage.eq("true")],
        "moth_record_without_butterfly_positive": m.loc[m.butterfly_record_positive.eq("false")],
    }
    results, pairs = {}, []
    for name, frame in groups.items():
        a, b = cohort_support(frame)
        results[name] = a
        if len(b):
            b = b.copy()
            b["cohort"] = name
            pairs.append(b)
    output.mkdir(parents=True)
    outcsv = output / "tier09_moth_original_photo_existingVIIRS_cell_status.csv.gz"
    write_gz(m[["photo_id", "observation_uuid", "year", "land_id", "cell_id",
                "area_km2", "species", "moth_database_pollinates",
                "walton2020_nocturnal_pollen_carriage", "butterfly_record_positive",
                "existingcell_source_hit", "viirs_eligible", "log_light"]], outcsv)
    if pairs:
        pd.concat(pairs, ignore_index=True).to_csv(
            output / "tier09_moth_supporting_island_years.csv", index=False)
    report = {
        "schema": SCHEMA, "source_only": True, "source_roster_photo_count": int(len(d)),
        "original_tier09_source_SHA256": sha256(source),
        "preverified_GSHHG_and_dedup_source_receipt_SHA256": sha256(receipt),
        "sameyear_original_VIIRS_verification": signed,
        "source_only_existing_cell_VIIRS_values_not_new_raster_sampling": True,
        "unmatched_cells_are_UNKNOWN_not_dark": True,
        "observational_light_span_threshold_log1p": LIGHT_MIN_LOG_SPAN,
        "cohorts": results, "new_RGB_downloads": 0, "new_numeric_ROI9": 0,
        "cohort_selection_pre_colour": True,
        "white_ROI_is_not_native_plant_white_prevalence": True,
        "GloBI_record_positive_is_not_moth_obligate": True,
        "ready_for_bulk_RGB": False,
        "candidate_status_SHA256": sha256(outcsv),
    }
    (output / "tier09_moth_VIIRS_existingcell_gate.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("ALAN_TIER09_EXISTING_SAMEYEAR_VIIRS_SOURCE_GATE",
          json.dumps(report, sort_keys=True), flush=True)
    return report

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    for year in YEARS:
        p.add_argument("--viirs" + str(year), type=Path, required=True)
        if year != 2024:
            p.add_argument("--receipt" + str(year), type=Path, required=True)
    a = p.parse_args()
    run(a.source, a.receipt,
        {y: getattr(a, "viirs" + str(y)) for y in YEARS},
        {y: getattr(a, "receipt" + str(y)) for y in YEARS if y != 2024},
        a.out)

if __name__ == "__main__":
    main()
