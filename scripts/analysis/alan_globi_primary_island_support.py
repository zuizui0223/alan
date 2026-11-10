#!/usr/bin/env python3
"""GloBI-first island plant support audit, without fd-based trait filtering.

Requires three frozen GitHub Actions artifact ZIP files. This source evaluates
recorded moth/butterfly flower links, NOT demonstrated effective pollination,
exclusive guilds, local island interactions or established naturalization.
No underlying restricted species-to-island rows are exported.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pandas as pd

SOURCES = {
    "GloBI": ("GloBI_stage2_signed_extension.csv.gz",
              "a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994"),
    "island_flora": ("fixed/canonical/input/chapter1_status_flora.csv.gz",
                     "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"),
    "island_photo": ("island_year_ROI_free_source_support.csv",
                     "acaf9287a29ef01e2fc79fd4b3450b4c09bb5eee16737e55b4326a80359cbdf8"),
}
MOTH = "moth_family_visit"
BUTTERFLY = "butterfly_family_visit"


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def locked_member(zip_path, source_name):
    member, expected_sha = SOURCES[source_name]
    with zipfile.ZipFile(zip_path) as archive:
        raw = archive.read(member)
    require(hashlib.sha256(raw).hexdigest() == expected_sha,
            "Source member SHA256 mismatch: " + source_name)
    return pd.read_csv(io.BytesIO(gzip.decompress(raw) if member.endswith(".gz") else raw),
                       dtype=str, keep_default_na=False)


def validate(globi, flora, photos):
    require(len(globi) == 61967, "Unexpected GloBI source frame")
    require({"plant_name", "guild", "relation", "study_title", "source_family"} <= set(globi),
            "GloBI schema drift")
    require(set(globi.guild) <= {MOTH, BUTTERFLY} and
            set(globi.relation) <= {"visitsFlowersOf", "pollinates"},
            "Unexpected GloBI interaction type or guild")
    require({"island_id", "accepted_species", "origin_status", "status_conflict"} <= set(flora),
            "Island status schema drift")
    require(not flora.duplicated(["island_id", "accepted_species"]).any(),
            "Duplicate island/species record")
    require(set(flora.origin_status) <= {"native", "introduced", "unresolved"},
            "Unexpected floristic origin status")
    require(set(flora.status_conflict) == {"False"},
            "Source island status conflicts")
    require({"island_id", "year", "photos_all_original_no_roi_gate",
             "photos_with_positive_moth_visit_record", "same_year_viirs_matched_photos",
             "light_source"} <= set(photos), "Photo source schema drift")
    require(not photos.duplicated(["island_id", "year"]).any(),
            "Duplicate photo island/year")
    require(set(photos.light_source) <= {"NO_same_year_light",
                                         "quarter_months_NOT_annual"},
            "Photo source includes unexpected light product; update contract")


def scenario(globi, flora, photos):
    flags = globi[["plant_name", "guild"]].drop_duplicates().groupby(
        "plant_name").guild.agg(frozenset)
    moth = set(flags.index[flags.map(lambda s: MOTH in s)])
    butterfly = set(flags.index[flags.map(lambda s: BUTTERFLY in s)])
    moth_record_only = set(flags.index[flags.map(lambda s: s == {MOTH})])
    butterfly_record_only = set(flags.index[flags.map(lambda s: s == {BUTTERFLY})])
    a = flora[flora.accepted_species.isin(moth)]
    b = flora[flora.accepted_species.isin(butterfly)]
    m_only = flora[flora.accepted_species.isin(moth_record_only)]
    b_only = flora[flora.accepted_species.isin(butterfly_record_only)]
    intro = a[a.origin_status.eq("introduced")]
    a_is = a.groupby("island_id").accepted_species.nunique()
    b_is = b.groupby("island_id").accepted_species.nunique()
    mo_is = m_only.groupby("island_id").accepted_species.nunique()
    bo_is = b_only.groupby("island_id").accepted_species.nunique()

    def matched_above(x, y, n):
        return int(sum(x.loc[i] >= n and y.loc[i] >= n
                       for i in (set(x.index) & set(y.index))))

    p_moth = photos.loc[
        pd.to_numeric(photos.photos_with_positive_moth_visit_record).gt(0), "island_id"]
    return {
        "source_interaction_rows": int(len(globi)),
        "distinct_moth_positive_plant_names": len(moth),
        "distinct_butterfly_positive_plant_names": len(butterfly),
        "moth_record_without_butterfly_record_NOT_true_exclusivity": len(moth_record_only),
        "butterfly_record_without_moth_record_NOT_true_exclusivity": len(butterfly_record_only),
        "plants_with_both_recorded_guilds": len(moth & butterfly),
        "exact_moth_plant_name_matches_in_island_flora": int(a.accepted_species.nunique()),
        "moth_plant_island_name_rows": int(len(a)),
        "islands_with_at_least_1_moth_linked_plant": int(a_is.ge(1).sum()),
        "islands_with_at_least_3_moth_linked_plants": int(a_is.ge(3).sum()),
        "islands_with_at_least_5_moth_linked_plants": int(a_is.ge(5).sum()),
        "islands_with_at_least_10_moth_linked_plants": int(a_is.ge(10).sum()),
        "islands_with_both_positive_guilds_at_least_3_plants_each_SHARED_ALLOWED":
            matched_above(a_is, b_is, 3),
        "islands_with_3_species_each_of_record_only_evidence_NOT_absence":
            matched_above(mo_is, bo_is, 3),
        "introduced_origin_moth_species_NOT_naturalization": int(intro.accepted_species.nunique()),
        "introduced_origin_moth_islands_NOT_naturalization": int(intro.island_id.nunique()),
        "introduced_origin_moth_island_name_rows_NOT_naturalization": int(len(intro)),
        "moth_flora_rows_by_origin": {s: int(a.origin_status.eq(s).sum())
                                      for s in ("native", "introduced", "unresolved")},
        "moth_flora_islands_overlapping_photo_island_ID_NOT_polygon_verified":
            len(set(a.island_id) & set(photos.island_id)),
        "moth_flora_islands_overlapping_positive_moth_photo_island_ID":
            len(set(a.island_id) & set(p_moth)),
    }


def photo_coverage(photos):
    p = photos.copy()
    for field in ("photos_all_original_no_roi_gate",
                  "photos_with_positive_moth_visit_record",
                  "same_year_viirs_matched_photos"):
        p[field] = pd.to_numeric(p[field], errors="raise")
        require(p[field].ge(0).all(), "Invalid negative photo count: " + field)
    group = p.groupby("island_id").agg(
        total=("photos_all_original_no_roi_gate", "sum"),
        moth=("photos_with_positive_moth_visit_record", "sum"))
    year_support = p[p.photos_with_positive_moth_visit_record > 0].groupby(
        "island_id").year.nunique()
    proxy = p[p.same_year_viirs_matched_photos > 0]
    require(proxy.light_source.eq("quarter_months_NOT_annual").all(),
            "Quarterly photo-cell proxy cannot be called annual VIIRS")
    return {
        "source_photo_islands": len(group),
        "source_photo_island_years": len(p),
        "all_original_plant_photo_observation_opportunities":
            int(p.photos_all_original_no_roi_gate.sum()),
        "positive_GloBI_moth_tagged_PLANT_photos_NOT_moth_visits":
            int(p.photos_with_positive_moth_visit_record.sum()),
        "islands_with_at_least_1_GloBI_moth_tagged_photo": int(group.moth.ge(1).sum()),
        "islands_with_at_least_5_GloBI_moth_tagged_photos": int(group.moth.ge(5).sum()),
        "islands_with_at_least_20_GloBI_moth_tagged_photos": int(group.moth.ge(20).sum()),
        "islands_with_positive_moth_tagged_photos_in_2_years":
            int(year_support.ge(2).sum()),
        "quarterly_proxy_light_island_years_NOT_annual": len(proxy),
        "quarterly_proxy_light_islands_NOT_annual": int(proxy.island_id.nunique()),
        "quarterly_proxy_light_moth_photo_island_years_NOT_annual":
            int(proxy.photos_with_positive_moth_visit_record.gt(0).sum()),
        "annual_VIIRS_eligible_island_years": None,
    }


def run(globi_zip, island_zip, photo_zip):
    g = locked_member(globi_zip, "GloBI")
    f = locked_member(island_zip, "island_flora")
    p = locked_member(photo_zip, "island_photo")
    validate(g, f, p)
    title = g.study_title.str.lower()
    is_inaturalist = title.str.contains("inaturalist", regex=False)
    is_trophiCH = title.str.contains("trophich", regex=False)
    scenarios = {
        "all_GloBI_bounded_stage2": g,
        "exclude_iNaturalist_interaction_sources": g.loc[~is_inaturalist],
        "exclude_iNaturalist_and_trophiCH": g.loc[~is_inaturalist & ~is_trophiCH],
        "non_iNaturalist_visitsFlowersOf_only":
            g.loc[~is_inaturalist & g.relation.eq("visitsFlowersOf")],
    }
    return {
        "schema": "alan_globi_primary_island_source_feasibility_20261111_v1",
        "provenance": {
            "GloBI": {"run": 41976552977, "artifact": 13391803415},
            "island": {"run": 34232450884, "artifact": 10058653212},
            "photo": {"run": 42171173575, "artifact": 13481850056},
            "members_SHA256": {k: v[1] for k, v in SOURCES.items()},
        },
        "source": {
            "all_GloBI_claim_rows": len(g),
            "moth_claim_rows": int(g.guild.eq(MOTH).sum()),
            "butterfly_claim_rows": int(g.guild.eq(BUTTERFLY).sum()),
            "iNaturalist_source_claim_rows": int(is_inaturalist.sum()),
            "trophiCH_source_claim_rows": int(is_trophiCH.sum()),
            "island_species_rows": len(f),
            "historical_island_IDs": int(f.island_id.nunique()),
        },
        "photo_source": photo_coverage(p),
        "scenarios": {name: scenario(table, f, p) for name, table in scenarios.items()},
        "source_limitation": {
            "GloBI_positive_visits_not_effective_nocturnal_pollination": True,
            "GloBI_recorded_only_not_true_pollinator_exclusivity": True,
            "global_plant_interaction_record_not_local_island_pollination": True,
            "introduced_status_not_naturalized": True,
            "photographic_plant_presence_not_insect_observation": True,
            "historical_island_ID_overlap_not_audited_geometry_equivalence": True,
            "annual_VIIRS_gecko_three_way_not_estimable": True,
            "no_protected_species_island_rows_published": True,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("globi-zip", "island-zip", "photo-zip", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    a = parser.parse_args()
    require(not a.output.exists(), "Refuse overwriting an existing result")
    report = run(a.globi_zip, a.island_zip, a.photo_zip)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, sort_keys=True,
                                   ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"source": report["source"], "photo_source": report["photo_source"],
                      "scenarios": report["scenarios"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
