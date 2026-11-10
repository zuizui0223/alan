"""Synthetic geographic and photo-documentation gates; no artificial biological observations."""
from __future__ import annotations
import unittest
import pandas as pd

from scripts.analysis.alan_island_functional_pollination_5species_exposure_census_20261110 import (
    archipelago_filter,visit_support,CONFIRMED)

class TrueIslandFunctionalGuildSourceGates(unittest.TestCase):
    def test_archipelago_excludes_cultivated_off_region_photo(self):
        # Canary source centroids ~28.5,-16; Mauritius ~-20.3,57.7
        sample=pd.DataFrame([
            {"island":"canary","year":2024,
             "latbin":int((28.5+90)/.05),"lonbin":int((-16+180)/.05)},
            {"island":"mauritius","year":2024,
             "latbin":int((-20.3+90)/.05),"lonbin":int((57.7+180)/.05)},
        ])
        canary=archipelago_filter(sample,"Canary")
        maur=archipelago_filter(sample,"Mauritius")
        self.assertEqual(list(canary.island),["canary"])
        self.assertEqual(list(maur.island),["mauritius"])

    def test_gecko_earlier_source_not_future_source(self):
        original=pd.DataFrame([
            {"island":"a","year":2019,"latbin":100,"lonbin":200},
            {"island":"a","year":2025,"latbin":100,"lonbin":200},
            {"island":"b","year":2025,"latbin":102,"lonbin":200}
        ])
        records=pd.DataFrame([
            {"island":"a","year":2020,"latbin":100,"lonbin":200,
             "actual_source_taxon_family":"Gekkonidae",
             "actual_GBIF_human_observation_records":3},
            {"island":"b","year":2023,"latbin":102,"lonbin":200,
             "actual_source_taxon_family":"Scincidae",
             "actual_GBIF_human_observation_records":2},
        ])
        genera=pd.DataFrame([
            {"island":"a","year":2020,"latbin":100,"lonbin":200,
             "original_genus":"Hemidactylus","original_source_records":3}
        ])
        a=visit_support(original,records,genera)
        self.assertEqual(a["original_islands_with_gecko_earlier_observation_years"],1)
        self.assertEqual(a["photo_sites_with_gecko_record_earlier_same_island"],1)
        self.assertEqual(a["photo_sites_with_Hemidactylus_genus_earlier_same_island"],1)
        self.assertEqual(a["photo_sites_with_gecko_record_same_cell_year"],0)

    def test_evidence_not_bloom_or_gbif_visits(self):
        self.assertEqual(len(CONFIRMED),8)
        self.assertEqual(sum("effective_native_day_gecko"==v["guild"] for v in CONFIRMED.values()),2)
        self.assertEqual(sum("effective_native_Gallotia_lizard_NOT_gecko"==v["guild"] for v in CONFIRMED.values()),1)
        self.assertEqual(sum("effective_native_Phelsuma_day_gecko"==v["guild"] for v in CONFIRMED.values()),1)
        self.assertTrue(all(v["evolutionary_ancestral_shift_proven"] is False for v in CONFIRMED.values()))

if __name__=="__main__":unittest.main()
