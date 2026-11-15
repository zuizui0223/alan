"""Dawn co-availability is not proven pollinator mismatch."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.alan_Grognuz2026_same_date_flower_resource_window_20261115 import (
  circadian_offset,dawn_window,site_effect,BOUNDS
)

class DawnFlowerResourceTests(unittest.TestCase):
    def test_sunrise_clock_wrap_previous_evening(self):
        self.assertAlmostEqual(circadian_offset(21*60+50,6*60+45),-535)
        self.assertAlmostEqual(circadian_offset(8*60,6*60+30),90)
        self.assertAlmostEqual(circadian_offset(0*60+5,6*60+45),-400)
    def test_partial_vs_full_open_bounds(self):
        # Silene partially open until +90min, fully open until -60min;
        # Malva starts opening -90min but fully open +30min.
        partial=dawn_window(np.array([90]),np.array([-90]),150)
        full=dawn_window(np.array([-60]),np.array([30]),150)
        self.assertEqual(float(partial[0]),180)
        self.assertEqual(float(full[0]),0)
        self.assertGreater(partial[0],full[0])
        self.assertEqual(float(dawn_window(np.array([-500]),np.array([-20]),150)[0]),0)
    def test_site_bootstrap_counts_true_sites_not_46_dates(self):
        site_codes=["L"+str(i) for i in range(6)]+["D"+str(i) for i in range(6)]
        treatments=["LIGHT"]*6+["DARK"]*6
        z=pd.DataFrame({
          "Study.site":site_codes,"LIGHT_TRMT":treatments,
          "silene_begin_close":[90]*6+[0]*6,
          "malva_begin_open":[0]*12
        })
        result=site_effect(z,"silene_begin_close","malva_begin_open",150,seed=42,B=100)
        self.assertEqual((result["independent_lit_sites"],result["independent_dark_sites"]),(6,6))
        self.assertEqual(result["lit_minus_dark_minutes"],90)
        self.assertEqual(result["site_resample95_descriptive_minutes"],[90,90])
        self.assertTrue(result["site_labels_not_recovered_as_randomized_pairs"])
    def test_four_uncertain_window_stage_definitions(self):
        self.assertEqual(len(BOUNDS),4)
        self.assertEqual(BOUNDS["conservative_both_fully_open"],
                         ("silene_begin_close","malva_full_open"))
        self.assertEqual(BOUNDS["permissive_partially_open"],
                         ("silene_full_closed","malva_begin_open"))
if __name__=="__main__":unittest.main()
