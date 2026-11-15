"""Source-v1 diel+seasonal no-pseudoreplication and deterministic arithmetic guards."""
import tempfile
import unittest
from pathlib import Path
from scripts.analysis.alan_Grognuz2026_original_site_diel_seasonal_20261115 import (
  EVENTS,clock_offset,site_diff_bootstrap,check_source
)

class TestGrognuz2026(unittest.TestCase):
    def test_sunrise_normalization_and_midnight_wrap(self):
        self.assertEqual(clock_offset(8*60,6*60),120)
        self.assertEqual(clock_offset(5*60,6*60),-60)
        self.assertEqual(clock_offset(0*60+30,20*60+30),240)
        self.assertEqual(clock_offset(19*60+30,20*60+30),-60)
    def test_site_bootstrap_not_plant_day_rows(self):
        lit=[100,110,120,90,105,115]
        dark=[-10,-20,10,0,-5,5]
        a=site_diff_bootstrap(lit,dark,20261115,4000)
        b=site_diff_bootstrap(lit,dark,20261115,4000)
        self.assertEqual(a,b)
        self.assertGreater(a[0],0)
        self.assertGreater(a[1],a[0])
        self.assertIsNone(site_diff_bootstrap([1],[2],13,20))
    def test_no_fallback_on_wrong_original_hash(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"fake.csv";p.write_text("fabricated source")
            with self.assertRaises(ValueError):check_source(p,"0"*64)
    def test_daily_event_definitions_not_pollinator_mismatch(self):
        self.assertEqual(len(EVENTS),4)
        self.assertEqual(EVENTS[0][0],"Silene noctiflora")
        self.assertEqual(EVENTS[1][0],"Malva moschata")
        self.assertEqual(EVENTS[0][2],"Time.of.sunrise")
        self.assertEqual(EVENTS[1][2],"Time.of.sunrise")

if __name__=="__main__":unittest.main()
