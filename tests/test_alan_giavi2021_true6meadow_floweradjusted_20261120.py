"""Synthetic original-meadow source counts and flower controls, no pseudo rep."""
import unittest
import pandas as pd
from scripts.analysis.alan_giavi2021_true6meadow_floweradjusted_moisture_20261120 import (
    estimate,site_signflip,source_support
)
def fixture():
    rows=[]
    for site in ["s1","s2","s3"]:
        for t in ["r1","r2"]:
            for plant,F in (("dry",4.),("medium",5.),("wet",8.)):
                flower={"dry":.1,"medium":-.3,"wet":.4}[plant]
                dr=-.25*F +.7*flower +(0.1 if site=="s1" else -.08)
                rows.append(dict(site=site,pair_key=site+"|"+t,plant=plant,
                    F_PLANTATT=F,F_Doring=F,delta_flower=flower,
                    delta_visits=dr,delta_v_Diptera=dr,delta_v_Hymenoptera=dr,
                    delta_v_Coleoptera=dr))
    return pd.DataFrame(rows)
class GiaviOriginalModel(unittest.TestCase):
    def test_original_source_pair_and_flower_slope(self):
        df=fixture()
        e=estimate(df,"F_PLANTATT","visits",True)
        self.assertEqual(e["source_true_site_pairs"],3)
        self.assertAlmostEqual(e["beta_F_per_ordinal_point"],-.25,places=8)
        self.assertAlmostEqual(e["beta_delta_log1p_flower"],.7,places=8)
    def test_one_real_meadow_flip_for_all_plant_repeats(self):
        df=fixture()
        e=site_signflip(df,"F_Doring","visits",True)
        self.assertEqual(e["exact_whole_meadow_label_flips"],8)
        self.assertTrue(0<=e["two_sided_abs_extreme_fraction"]<=1)
        self.assertTrue(e["not_original_randomized_unit_test"])
    def test_binary_wet_gate_requires_four_wet_plant_taxa(self):
        d=fixture()
        s=source_support(d,"F_PLANTATT")
        self.assertEqual(s["sites_with_wet_and_drier"],3)
        self.assertFalse(s["confirmatory_moisture_moderator_support"])
        self.assertEqual(s["status"],"SOURCE_WET_DRY_COMMON_SUPPORT_STOP")
    def test_no_inferential_rank_without_meadow_wet_dry_overlap(self):
        d=fixture()
        d.loc[d.site=="s1","F_PLANTATT"]=5.
        d.loc[d.site=="s2","F_PLANTATT"]=4.
        d.loc[d.site=="s3","F_PLANTATT"]=8.
        s=source_support(d,"F_PLANTATT")
        self.assertEqual(s["sites_with_wet_and_drier"],0)
        self.assertFalse(s["confirmatory_moisture_moderator_support"])
if __name__=="__main__":
    unittest.main()
