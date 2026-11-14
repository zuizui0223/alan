#!/usr/bin/env python3
"""Small source-publication positive controls for original published Table I."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_primary_Marques2015_bat_bird_reproductive_table_audit_20261114 import (
    PUB_FRUIT_COUNTS, PUB_SEEDS_PER_FRUIT, STUDY_ID, SOURCE_DOI_AS_PRINTED, require
)


class PrimaryBatHummingbirdData(unittest.TestCase):
    def test_published_original_table_fruit_reproduction(self):
        day=PUB_FRUIT_COUNTS["day pollination"][0]/PUB_FRUIT_COUNTS["day pollination"][1]
        night=PUB_FRUIT_COUNTS["night pollination"][0]/PUB_FRUIT_COUNTS["night pollination"][1]
        open_=PUB_FRUIT_COUNTS["open pollination"][0]/PUB_FRUIT_COUNTS["open pollination"][1]
        self.assertAlmostEqual(day,.1002475,places=6)
        self.assertAlmostEqual(night,.821066,places=6)
        self.assertAlmostEqual(open_,.7785714,places=6)
        self.assertGreater(night,day)
        self.assertGreater(PUB_SEEDS_PER_FRUIT["night pollination"],
                           PUB_SEEDS_PER_FRUIT["day pollination"])
        self.assertEqual(STUDY_ID,"Marques_2015")

    def test_article_public_identifier_is_as_printed_not_autofixed(self):
        self.assertEqual(SOURCE_DOI_AS_PRINTED,"10.590/0001-3765201520140340")
        with self.assertRaises(ValueError):
            require(False,"Do not assume 2010 actual animal observation year")


if __name__=="__main__":
    unittest.main()
