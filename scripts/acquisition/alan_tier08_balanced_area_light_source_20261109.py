#!/usr/bin/env python3
"""ALAN tier08: genuine new source UUIDs and pre-colour matched island/light pairs.

Uses only the immutable MOIS 93,209,351-original-photo source registry's disjoint
tier08 (observation-hash probability interval [0.0256, 0.0512)). No flower
colours, native status, unobserved plants or new satellite values are imputed.
This intake is distinct from historical-main and tier06+tier07 image data.
A 256-photo conditional hash reference arm and a separate area/light-paired
source-targeting arm are frozen BEFORE any tier08 RGB/ROI measurement. They must
NOT be pooled as a single unbiased ecological or native-flora sample.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.acquisition import alan_tier07_1193829_independent_island_source_20261108 as previous
from scripts.analysis import alan_island_area_light_20261108 as areas
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel

REGISTRY_SHA = "4ed5139e510ba8442d4cdb03f2b769d94fac942832b1d544d7a7925bc431255b"
PRIOR06_SHA = "d0b7f7d9f117036d05fbc7ff8b1e62847ef05e7bbc07efc659c9f6cf895d30f9"
PRIOR07_SHA = "e38161652d3ec85928fa0f4fad91d70c123a7cb2af32f65f199f42504f76e056"
SCHEMA = "ALAN_tier08_precolour_independent_uuid_balanced_light_pair_source_v1"
TIER = 8
SEED = "ALAN_TIER08_NEW_INDEPENDENT_SOURCE_20261109"
PHOTO_SOURCE = "mois_complete_full93m_original_photo_target_registry_v1"
YEARS = (2019, 2021, 2023, 2024, 2025)
REFERENCE_N = 256
MAX_PAIR_ISLANDS_PER_AREA_THIRD = 12
MAX_PHOTOS_PER_SITE = 12

def require(check, msg):
    if not bool(check):
        raise ValueError("ALAN_TIER08_SOURCE: " + msg)

def rank(label, item):
    return hashlib.sha256((SEED + "|" + label + "|" + str(item)).encode()).hexdigest()

def choose_per_observation(rows, known):
    """Single photo per genuinely new UUID, chosen colour-blind and order-independently."""
    best = {}
    for row in rows:
        obs = str(row["observation_uuid"])
        if obs in known:
            continue
        photo = str(row["photo_id"])
        h = rank("photo-selection", obs + "|" + photo)
        if obs not in best or h < best[obs][0]:
            best[obs] = (h, dict(row))
    return [v[1] for v in best.values()]

def validate_registry(registry):
    require(panel.sha(registry) == REGISTRY_SHA,
            "Signed complete full93m registry SHA mismatch")
    manifest = json.loads(Path(registry).read_text())
    require(manifest.get("schema") == PHOTO_SOURCE and
            manifest.get("exact_original_photo_IDs_in_all_disjoint_tiers") == 93209351 and
            manifest.get("source_only_NO_photo_pixels_or_hue_measured") is True,
            "Source is not the original photo-only full93m registry")
    hashes = manifest["published_tier_archive_SHA256"]
    names = [f"tier-{TIER:02d}-group-{g:02d}.tar" for g in range(16)]
    require(all(len(hashes.get(n, "")) == 64 for n in names),
            "Incomplete 16 original signed tier08 archive hashes")
    required_n = int(manifest["new_target_counts_by_tier"][str(TIER)])
    require(required_n > 1193829,
            "Tier08 original disjoint photo sampling-frame count unexpectedly small")
    return {n: hashes[n] for n in names}, required_n

def read_tier08(source, registry, known):
    """Stream all 16 original signed tar files; never download image pixels."""
    hashes, required_n = validate_registry(registry)
    original_photo_ids = set()
    chosen = {}
    group_counts = []
    for group in range(16):
        name = f"tier-{TIER:02d}-group-{group:02d}.tar"
        path = Path(source) / name
        require(path.is_file() and panel.sha(path) == hashes[name],
                "Original tier08 tar missing or checksum drift: " + name)
        parts = set()
        count = 0
        with tarfile.open(path, "r:") as archive:
            for member in archive:
                if not member.isfile() or not member.name.endswith(".csv.gz"):
                    continue
                match = re.search(r"part-(\d{3})\.csv\.gz$", member.name)
                require(match is not None, "Unexpected source member in tar")
                part = int(match.group(1))
                require(part // 8 == group and part not in parts,
                        "Source partition duplicated or assigned to wrong group")
                parts.add(part)
                with gzip.GzipFile(fileobj=archive.extractfile(member)) as zipped:
                    reader = csv.DictReader(io.TextIOWrapper(
                        zipped, encoding="utf-8", newline=""))
                    require(tuple(reader.fieldnames or ()) == previous.COLUMNS,
                            "Original photo-index columns drift")
                    for row in reader:
                        photo = row["photo_id"]
                        require(photo.isdecimal() and
                                photo not in original_photo_ids and
                                row["tier"] == str(TIER),
                                "Repeated original photo or mismatched tier")
                        original_photo_ids.add(photo)
                        count += 1
                        obs = row["observation_uuid"]
                        if obs in known:
                            continue
                        score = rank("photo-selection", obs + "|" + photo)
                        if obs not in chosen or score < chosen[obs][0]:
                            chosen[obs] = (score, dict(row))
        require(len(parts) == 8, "Eight authentic source shards required per group")
        group_counts.append(count)
    require(sum(group_counts) == len(original_photo_ids) == required_n,
            "Source tier08 photo count does not match frozen 93m manifest")
    rows = [v[1] for v in chosen.values()]
    require(len(rows) == len(chosen) and
            len(set(x["observation_uuid"] for x in rows)) == len(rows),
            "One original photograph per never-before-used observation UUID")
    audit = {
        "disjoint_tier08_original_photo_IDs": required_n,
        "new_original_observation_UUIDs_after_canonical_tier06_tier07_exclusion": len(rows),
        "signed_16_original_tar_group_photo_counts": group_counts,
        "tier08_observation_hash_probability_interval": [0.0256, 0.0512],
        "full_registry_sha256": REGISTRY_SHA,
        "one_colour_blind_original_photo_per_uuid": True,
        "RGB_original_image_bytes_downloaded": 0,
    }
    return rows, audit

def source_env_light(paths, canonical, wcvp):
    # Original observation-year photo-cell support only; no VIIRS extrapolation.
    parts = []
    for year in YEARS:
        if year in (2019, 2021, 2023):
            full, _ = panel.year_historic(year, paths, canonical, wcvp, return_all=True)
        elif year == 2024:
            full, _ = panel.year2024(paths, canonical, wcvp, return_all=True)
        else:
            full, _ = panel.year2025(paths, canonical, wcvp, return_all=True)
        require(full.target_id.is_unique and "viirs_median" in full,
                "Original same-year VIIRS source not uniquely indexed")
        parts.append(full[["target_id", "viirs_median"]].copy())
    env = pd.concat(parts, ignore_index=True)
    require(env.target_id.is_unique, "Same source photo cell/year repeated")
    return env

def build_paired_queue(master):
    """Pair two genuinely distinct sites within an island and SAME source year.

    Select islands across three observed area thirds. Select minimum and
    maximum observed year-matched VIIRS sites, never years of the same cell.
    This is adaptive geography sampling, not a colour response test.
    """
    m = master.copy()
    m["site_id"] = m.target_id.astype(str).str.split("|", regex=False).str[1]
    require(m.site_id.notna().all(), "Missing original 0.05 degree photo site")
    site = m.groupby(["land_id", "year", "target_id", "site_id"], observed=True).agg(
        attempts=("observation_id", "nunique"),
        taxa=("wcvp_accepted_plant_name_id", "nunique"),
        log_light=("log_light", "first"),
        area_km2=("area_km2", "first"),
    ).reset_index()
    # Actual photo success is still unknown; provide enough attempts and
    # accepted taxa for potential strict measured-photo qualification.
    eligible = site.loc[site.attempts.ge(10) & site.taxa.ge(4)].copy()
    candidates = []
    for (land, year), d in eligible.groupby(["land_id", "year"], observed=True):
        if d.site_id.nunique() < 2:
            continue
        low = d.sort_values(["log_light", "target_id"], kind="stable").iloc[0]
        high = d.sort_values(["log_light", "target_id"], ascending=[False, True],
                             kind="stable").iloc[0]
        span = float(high.log_light - low.log_light)
        if low.site_id == high.site_id or span < 0.05:
            continue
        candidates.append({
            "land_id": land, "year": int(year),
            "low_target_id": low.target_id, "high_target_id": high.target_id,
            "log_light_span": span,
            "area_km2": float(low.area_km2),
            "min_site_opportunities": int(min(low.attempts, high.attempts)),
            "pair_rank": rank("pair", str(land) + "|" + str(year)),
        })
    if not candidates:
        return m.iloc[0:0].copy(), {
            "candidate_area_light_paired_islands": 0,
            "selected_paired_islands": 0,
            "selected_paired_photos": 0,
            "selected_pairs_by_area_third": {},
        }
    pairs = pd.DataFrame(candidates).sort_values(
        ["min_site_opportunities", "log_light_span", "pair_rank"],
        ascending=[False, False, True], kind="stable"
    ).drop_duplicates("land_id", keep="first")
    pairs = pairs.sort_values(["area_km2", "land_id"], kind="stable").copy()
    pairs["area_third"] = pd.cut(
        pairs.area_km2.rank(method="first", pct=True),
        bins=[0.0, 1/3, 2/3, 1.0],
        labels=["small", "medium", "large"], include_lowest=True)
    selected_pairs = pd.concat([
        d.sort_values(["min_site_opportunities", "pair_rank"],
                      ascending=[False, True], kind="stable"
                      ).head(MAX_PAIR_ISLANDS_PER_AREA_THIRD)
        for _, d in pairs.groupby("area_third", observed=True)
    ], ignore_index=True)
    selected = []
    for p in selected_pairs.itertuples(index=False):
        for target in (p.low_target_id, p.high_target_id):
            sub = m.loc[m.target_id.eq(target)].copy()
            sub["_hash"] = sub.observation_id.map(lambda obs: rank("within-site", obs))
            sub = sub.sort_values("_hash", kind="stable")
            first = sub.drop_duplicates("wcvp_accepted_plant_name_id", keep="first")
            rest = sub.loc[~sub.observation_id.isin(set(first.observation_id))]
            out = pd.concat([first, rest], ignore_index=True).head(MAX_PHOTOS_PER_SITE)
            out["acquisition_arm"] = "area_stratified_sameyear_two_site_light_pair"
            out["area_third"] = str(p.area_third)
            out["source_light_span"] = float(p.log_light_span)
            selected.append(out.drop(columns=["_hash"]))
    queue = pd.concat(selected, ignore_index=True)
    require(queue.observation_id.is_unique and queue.photo_id.is_unique and
            queue.groupby("land_id").site_id.nunique().ge(2).all(),
            "Targeted paired-site queue repeats UUIDs or lacks distinct sites")
    return queue, {
        "candidate_area_light_paired_islands": int(len(pairs)),
        "selected_paired_islands": int(selected_pairs.land_id.nunique()),
        "selected_paired_site_years": int(queue.target_id.nunique()),
        "selected_paired_photos": int(len(queue)),
        "selected_pairs_by_area_third": {
            str(x): int(n) for x, n in
            selected_pairs.groupby("area_third", observed=True).size().items()
        },
        "min_required_observed_light_span_log1p_viirs": 0.05,
        "precolour_min_source_attempts_per_site": 10,
        "precolour_min_wcvp_accepted_taxa_per_site": 4,
    }

def execute(paths, coast, tier_source, tier_registry, prior06, prior07, out):
    output = Path(out)
    require(not output.exists(), "Do not overwrite frozen tier08 original source")
    require(panel.sha(prior06) == PRIOR06_SHA and
            panel.sha(prior07) == PRIOR07_SHA,
            "Prior original tier06/tier07 UUID manifests have changed")
    p6 = pd.read_csv(prior06, compression="gzip", dtype=str,
                     usecols=["observation_id", "photo_id"])
    p7 = pd.read_csv(prior07, compression="gzip", dtype=str,
                     usecols=["observation_id", "photo_id"])
    require(len(p6) == 8942 and len(p7) == 27507 and
            p6.observation_id.is_unique and p7.observation_id.is_unique and
            not set(p6.observation_id).intersection(set(p7.observation_id)),
            "Prior measured-photo cohort source identities drift")
    canonical, wcvp = panel.load_frozen_photos(paths)
    known = set(canonical.observation_id.astype(str)) | (
        set(p6.observation_id) | set(p7.observation_id))
    rows, audit = read_tier08(tier_source, tier_registry, known)
    geo, geo_audit = previous.source_gshhg(rows, coast)
    ready = previous.environmental_support(geo, canonical, wcvp, paths)
    island = ready.land_class.eq("island")
    valid = ready.has_original_year_ALAN_environment
    accepted = ready.accepted_taxon
    loader = ready.supported_image_extension & ready.source_loader_supported
    master = ready.loc[island & valid & accepted & loader].copy()
    require(master.observation_uuid.is_unique and master.unit_id.is_unique and
            not set(master.observation_uuid).intersection(known) and
            len(master) >= 256, "No genuine tier08 novel island photos to freeze")
    master["photo_id"] = master.unit_id.astype(str)
    master["observation_id"] = master.observation_uuid.astype(str)
    master["source_arm"] = "alan_93m_tier08_new_independent_island_original_photo"
    master["photo_url_large"] = master.apply(
        lambda r: "https://inaturalist-open-data.s3.amazonaws.com/photos/"
                  + str(r.photo_id) + "/large." + str(r.image_extension).lower(), axis=1)
    master["selection_probability"] = ""  # Enriched conditional frame, not global π.
    land_areas = areas.source_land_area(coast)
    island_areas = land_areas.loc[land_areas.land_class.eq("island"),
                                   ["land_id", "area_km2"]]
    master = master.merge(island_areas, on="land_id", validate="many_to_one")
    env = source_env_light(paths, canonical, wcvp)
    master = master.merge(env, on="target_id", validate="many_to_one")
    require(master.viirs_median.notna().all() and
            master.area_km2.gt(0).all(), "New observations lack real year light/area")
    master["log_light"] = np.log1p(
        pd.to_numeric(master.viirs_median).clip(lower=0))
    master = master.sort_values(["land_id", "target_id", "observation_id"],
                                kind="stable").reset_index(drop=True)
    paired, pair_audit = build_paired_queue(master)
    untouched = master.loc[
        ~master.observation_id.isin(set(paired.observation_id))].copy()
    untouched["_rank"] = untouched.observation_id.map(lambda obs: rank("reference", obs))
    reference = untouched.sort_values("_rank", kind="stable").head(REFERENCE_N).copy()
    require(len(reference) == REFERENCE_N,
            "Not enough new independent original reference observation UUIDs")
    reference["acquisition_arm"] = "hash_random_within_envmatched_tier08_reference"
    reference["area_third"] = ""
    reference["source_light_span"] = np.nan
    reference = reference.drop(columns=["_rank"])
    queue = pd.concat([reference, paired], ignore_index=True)
    require(queue.observation_id.is_unique and queue.photo_id.is_unique,
            "Acquisition arms must use disjoint original observations and photos")
    output.mkdir(parents=True)
    keep = ["photo_id", "observation_id", "inat_taxon_id",
            "wcvp_accepted_plant_name_id", "land_id", "target_id", "year",
            "observed_on", "latitude", "longitude", "photo_license",
            "image_extension", "photo_url_large", "source_arm",
            "selection_probability", "area_km2", "viirs_median", "log_light"]
    master[keep].to_csv(
        output / "tier08_all_new_envmatched_island_observations.csv.gz",
        compression={"method": "gzip", "mtime": 0}, index=False)
    queued_cols = keep + ["acquisition_arm", "area_third", "source_light_span"]
    queue[queued_cols].to_csv(
        output / "tier08_precolour_reference_and_area_light_paired_queue.csv.gz",
        compression={"method": "gzip", "mtime": 0}, index=False)
    counts = {
        "new_island_envmatched_original_observations": int(len(master)),
        "new_island_envmatched_GSHHG_islands": int(master.land_id.nunique()),
        "new_island_envmatched_source_cell_years": int(master.target_id.nunique()),
        "all_original_new_tier08_observation_UUIDs": int(audit[
            "new_original_observation_UUIDs_after_canonical_tier06_tier07_exclusion"]),
        "sameyear_environment_missing_for_new_accepted_island_observations": int(
            (island & accepted & ~valid & loader).sum()),
        "reference_actual_original_photos_unmeasured": int(len(reference)),
        "source_paired_actual_original_photos_unmeasured": int(len(paired)),
        "total_precolour_source_photo_queue": int(len(queue)),
        "new_source_image_pixels_downloaded": 0,
        "new_flower_ROI9_measurements": 0,
    }
    report = {
        "schema": SCHEMA,
        "observation_counts": counts,
        "source_registry_tier08_audit": audit,
        "source_GSHHG_audit": geo_audit,
        "paired_geographic_original_opportunities": pair_audit,
        "upstream_photo_source_date": "2026-09-27",
        "original_reference_source_probability_tier_interval": [0.0256, 0.0512],
        "reference_draw_within_remaining_eligible_UUIDs": "256 fixed smallest SHA256 ranks, independent of flower colour",
        "reference_further_conditional_source_random_fraction": float(
            REFERENCE_N / len(untouched)),
        "selection_policy": "reference and targeted arms separately labelled; conditional sampling not globally representative",
        "photo_attempts_are_not_numerical_flower_ROI_success": True,
        "no_unmeasured_colour_inferred": True,
        "prior_historical_main_and_tier06_tier07_original_photo_models_unchanged": True,
        "master_sha256": panel.sha(
            output / "tier08_all_new_envmatched_island_observations.csv.gz"),
        "queue_sha256": panel.sha(
            output / "tier08_precolour_reference_and_area_light_paired_queue.csv.gz"),
    }
    (output / "ALAN_tier08_immutable_precolour_source_receipt.json").write_text(
        json.dumps(report, sort_keys=True, indent=2) + "\n")
    print("ALAN_NEW_INDEPENDENT_TIER08_AREA_LIGHT_MATCHED_PRECOLOUR_SOURCE_VERIFIED",
          json.dumps(report, sort_keys=True), flush=True)
    return report

def main():
    p = argparse.ArgumentParser(description=__doc__)
    keys = ("canonical", "wcvp",
      "geo2019", "viirs2019", "pop2019",
      "geo2021", "viirs2021", "pop2021",
      "geo2023", "viirs2023", "pop2023",
      "geo2024", "cn2024", "uc2024", "uc2024receipt",
      "geo2025", "geo2025receipt", "viirs2025")
    for key in keys:
        p.add_argument("--" + key.replace("_", "-"), dest=key,
                       type=Path, required=True)
    for name in ("coast", "tier-source", "tier-registry",
                 "prior-tier06", "prior-tier07", "out"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    execute({k: getattr(a, k) for k in keys}, a.coast,
            a.tier_source, a.tier_registry,
            a.prior_tier06, a.prior_tier07, a.out)

if __name__ == "__main__":
    main()
