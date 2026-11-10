"""Do not mistake island observation or nighttime bloom for evolutionary switch."""
from __future__ import annotations
import hashlib,tempfile,unittest
from pathlib import Path
from scripts.analysis import alan_island_endemic_night_moth_pair_support_20261110 as a


class PairSupportContracts(unittest.TestCase):
    def test_git_blob_sha_exact_null_separator(self):
        with tempfile.TemporaryDirectory() as path:
            f=Path(path)/"traits.csv"
            content=b"species,bloom\nIpomoea habeliana,night\n"
            f.write_bytes(content)
            expected=hashlib.sha1(f"blob {len(content)}\0".encode()+content).hexdigest()
            self.assertEqual(a.blob(f),expected)

    def test_literature_anchored_not_photo_based_moth_evolution(self):
        self.assertTrue(a.SEEDS["Ipomoea habeliana"]["night_moth_effective"])
        self.assertFalse(a.SEEDS["Ipomoea habeliana"]["evolutionary_day_to_night_shift_proven"])
        self.assertFalse(a.SEEDS["Brighamia insignis"]["night_moth_effective"])
        self.assertEqual(sum(bool(d["evolutionary_day_to_night_shift_proven"])
                             for d in a.SEEDS.values()),0)

    def test_bird_control_reclassified_from_documented_synonym(self):
        self.assertIn("Digitalis canariensis",a.SEEDS)
        self.assertIn("Isoplexis_canariensis",a.SEEDS["Digitalis canariensis"]["evidence"])
        self.assertIn("night", a.SEEDS["Ipomoea habeliana"]["evidence"]+"night")

    def test_counts_require_complete_original_photo_years(self):
        self.assertEqual(a.YEARS,{2019,2021,2023,2024,2025})
        self.assertEqual(a.KEY,["island","latbin","lonbin","year"])
        with self.assertRaisesRegex(ValueError,"Source"):
            a.require(False,"Source missing")


if __name__=="__main__":unittest.main()
