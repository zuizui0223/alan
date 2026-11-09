"""Synthetic tests; no hypothetical original RGB measurements are asserted as evidence."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_pooled_plus_tier08_910_original_photo_20261109 as m
from scripts.analysis import alan_pooled_historical_new36449_photoROI9_20261109 as base

class Tier08ExactPhotoIdPoolTests(unittest.TestCase):
    @staticmethod
    def photo(id,obs,arm,taxon,white,cell="2024|a|9"):
        return dict(photo_id=id,observation_id=obs,photo_source=arm,
                    target_id=cell,land_id="9",year=2024,
                    wcvp_accepted_plant_name_id=taxon,
                    real_numeric_ROI9=white is not None,
                    source_photo_white_fraction=white if white is not None else np.nan)

    def test_repeated_exact_photo_ids_removed_before_old_new_tier08_combination(self):
        prior=pd.DataFrame([
             self.photo("100","old","historical_main","A",.3),
             self.photo("101","other","new","B",.5)])
        tier08=pd.DataFrame([
             self.photo("100","newUUID","tier08","A",.3),
             self.photo("102","freshUUID","tier08","C",None)])
        kept,audit=m.same_photo_reconciliation(prior,tier08)
        self.assertEqual(audit["tier08_original_photo_ids_already_in_prior"],1)
        self.assertEqual(kept.photo_id.tolist(),["102"])
        self.assertEqual(len(set(prior.photo_id)&set(kept.photo_id)),0)

    def test_disagreeing_taxon_with_same_original_photo_id_rejected(self):
        prior=pd.DataFrame([self.photo("100","old","historical_main","A",.3)])
        tier08=pd.DataFrame([self.photo("100","new","tier08","B",.3)])
        with self.assertRaisesRegex(ValueError,"disagrees"):
            m.same_photo_reconciliation(prior,tier08)

    def test_true_source_rois_add_new_strict_cell_without_imputing_unknown(self):
        old=pd.DataFrame([self.photo("1","old","historical_main","A",.8),
                          self.photo("2","old2","historical_main","A",.6)])
        new=pd.DataFrame([self.photo("3","new","new-tier","B",.2),
                          self.photo("4","new2","new-tier","C",None)])
        tier08=pd.DataFrame([self.photo("5","tier8","tier08","C",.4)])
        earlier=pd.concat([old,new],ignore_index=True)
        self.assertEqual(len(base.photo_cell_stats(earlier)[1]),0)
        all_=pd.concat([earlier,tier08],ignore_index=True)
        trial,strict=base.photo_cell_stats(all_)
        self.assertEqual(len(strict),1)
        self.assertEqual(trial.photo_attempts.iloc[0],5)
        self.assertEqual(trial.strict9_success.iloc[0],4)
        self.assertAlmostEqual(strict.white_pp.iloc[0],100*(.7+.2+.4)/3)
        balanced,count=m.three_arm_equal_source_photo_white(old,new,tier08,strict)
        self.assertEqual(count["cells_with_3_measured_sources"],1)
        self.assertAlmostEqual(balanced.white_pp.iloc[0],100*(.7+.2+.4)/3)

    def test_unknown_not_used_as_coloured_species(self):
        old=pd.DataFrame([self.photo("1","o1","historical_main","A",.8)])
        new=pd.DataFrame([self.photo("2","o2","new-tier","B",None)])
        tier08=pd.DataFrame([self.photo("3","o3","tier08","C",.4),
                             self.photo("4","o4","tier08","D",None),
                             self.photo("5","o5","tier08","E",None)])
        strict=base.photo_cell_stats(pd.concat([old,new,tier08]))[1]
        self.assertEqual(len(strict),0)

if __name__=="__main__":
    unittest.main()
