"""Original six-meadow source name and exact ordinal null guard."""
import unittest
from scripts.analysis.alan_giavi2021_plantatt_moisture_dayvisits_20261120 import (
 source_name,exact,matches,add_doring_only_records
)
class TestGiaviNames(unittest.TestCase):
    def test_exact_name_keeps_sp_substring_in_real_binomial(self):
        self.assertEqual(source_name("Heracleum_sphondylium"),"Heracleum sphondylium")
        self.assertIsNone(source_name("Centaurea_sp"))
        self.assertIsNone(source_name("Valeriana_officinalis_aggr."))
        self.assertIsNone(source_name("Erigeron_annuus_s.l."))
        self.assertIsNone(source_name("Rubus_sp"))
    def test_exact_null_with_ties_not_insect_sample_size(self):
        row=[{"Kew_plantatt_F":float(x),
             "source_lit_minus_dark_day_visitation_log1p":float(y)}
              for x,y in ((4,2),(4,1),(8,-1),(8,-2))]
        z=exact(row,"Kew_plantatt_F")
        self.assertEqual(z["distinct_source_label_arrangements"],6)
        self.assertAlmostEqual(z["observed_Spearman_rho"],-0.8944271909999159)
        self.assertEqual(z["wet_F_gt7"],2)
        self.assertEqual(z["physical_meadow_light_pairs_total"],6)
    def test_true_missing_PLANTATT_not_replaced_by_Doring(self):
        raw=("taxon_name_binom,F_PLANTATT,F_Doring,taxonomic_status,aggregate\n"
             "Cirsium oleraceum,NA,7,Accepted,\n").encode("latin1")
        z=add_doring_only_records(raw,{})
        self.assertIsNone(z["Cirsium oleraceum"]["F_PLANTATT"])
        self.assertEqual(z["Cirsium oleraceum"]["F_Doring"],"7")
if __name__=="__main__":
    unittest.main()
