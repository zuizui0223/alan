"""Synthetic-only tests: true original-photo pool, never fake ecological results."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis import alan_pooled_historical_new36449_photoROI9_20261109 as m

class ActualSourcePhotoPoolTests(unittest.TestCase):
    @staticmethod
    def row(photo,uuid,source,taxon,white,site="2024|10_10|9"):
        return dict(photo_id=str(photo),observation_id=str(uuid),
                    target_id=site,land_id="9",year=2024,
                    wcvp_accepted_plant_name_id=taxon,
                    real_numeric_ROI9=white is not None,
                    source_photo_white_fraction=(
                       float(white) if white is not None else np.nan),
                    photo_source=source)

    def test_full_pool_new_strict_site_from_two_individually_ineligible_arms(self):
        old=pd.DataFrame([
           self.row(1,"o1",m.SOURCE_OLD,"A",.9),
           self.row(2,"o2",m.SOURCE_OLD,"A",.7),
           self.row(3,"o3",m.SOURCE_OLD,"B",None)])
        new=pd.DataFrame([
           self.row(4,"n1",m.SOURCE_NEW,"B",.2),
           self.row(5,"n2",m.SOURCE_NEW,"C",.1)])
        self.assertEqual(len(m.photo_cell_stats(old)[1]),0)
        self.assertEqual(len(m.photo_cell_stats(new)[1]),0)
        both,eligible=m.photo_cell_stats(pd.concat([old,new]))
        self.assertEqual(len(both),len(eligible))
        self.assertEqual(int(eligible.photo_attempts.iloc[0]),5)
        self.assertEqual(int(eligible.strict9_success.iloc[0]),4)
        self.assertEqual(int(eligible.accepted_species_n.iloc[0]),3)
        self.assertAlmostEqual(eligible.white_pp.iloc[0],100*(.8+.2+.1)/3)
        e=eligible[["target_id","land_id","year","white_pp",
                    "photo_attempts","strict9_success","accepted_species_n"]]
        b,a=m.per_source_balanced_white(old,new,e)
        self.assertEqual(a["pooled_strict_cells_with_both_source_numeric_outcomes"],1)
        self.assertAlmostEqual(float(b.white_pp.iloc[0]),100*(.8+.15)/2)

    def test_observation_uuid_in_multiple_sources_is_invalid(self):
        photos=pd.DataFrame([
           self.row(1,"same",m.SOURCE_OLD,"A",.2),
           self.row(2,"same",m.SOURCE_NEW,"B",.6)])
        with self.assertRaisesRegex(ValueError,"Repeated independent"):
            m.photo_cell_stats(photos)

    def test_true_photo_unknown_not_converted_to_white(self):
        x=pd.DataFrame([
           self.row(1,"a",m.SOURCE_OLD,"A",.7),
           self.row(2,"b",m.SOURCE_OLD,"B",None),
           self.row(3,"c",m.SOURCE_NEW,"C",None),
           self.row(4,"d",m.SOURCE_NEW,"D",.1),
           self.row(5,"e",m.SOURCE_NEW,"E",.3)])
        allsites,eligible=m.photo_cell_stats(x)
        self.assertEqual(allsites.photo_attempts.iloc[0],5)
        self.assertEqual(allsites.strict9_success.iloc[0],3)
        self.assertEqual(eligible.accepted_species_n.iloc[0],3)
        self.assertAlmostEqual(eligible.white_pp.iloc[0],100*(.7+.1+.3)/3)

    def test_ambiguous_0p05_degree_landmass_is_never_assigned(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/"geo.csv"
            pd.DataFrame([
              dict(cell_id="1_1",target_id="2024|1_1|7",land_id="7",land_class="island"),
              dict(cell_id="1_1",target_id="2024|1_1|8",land_id="8",land_class="island"),
              dict(cell_id="2_2",target_id="2024|2_2|9",land_id="9",land_class="island"),
              dict(cell_id="3_3",target_id="2024|3_3|main",land_id="main",land_class="mainland")
            ]).to_csv(p,index=False)
            safe,audit=m.unique_original_gev(p,2024)
            self.assertEqual(safe.cell_id.tolist(),["2_2"])
            self.assertEqual(audit["ambiguous_original_geographic_rows"],2)

    def test_environment_missing_is_not_made_up(self):
        s=pd.DataFrame([dict(target_id="2024|x|1",land_id="1",year=2024,
                            white_pp=20.,photo_attempts=5)])
        env=pd.DataFrame([dict(target_id="2024|x|1",land_id="1",year=2024,
                               environment_complete=False,land_class="island",
                               log_light=np.nan,log_pop=np.nan,log_rain=np.nan,
                               bio5=np.nan,elevation=np.nan)])
        e,a=m.environment_qualified(s,env)
        self.assertEqual(len(e),0)
        self.assertEqual(a["original_strict_with_incomplete_environment"],1)

if __name__=="__main__":
    unittest.main()
