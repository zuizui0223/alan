"""Strict pairing of Page et al 2021 original SVE by true study and plant."""
import unittest
from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import (
    compare, QKEY, VKEY, PKEY
)
def row(study,plant,year,guild,metric="stigma pollen",
        unit="grains",visitmetric="visits/hour",q="10"):
    return {"studyID":study,"plant.species":plant,"study.year":year,
            "pollinator.group":guild,"pollinator.taxon":guild+"_sp",
            "effectiveness.metric":metric,"effectiveness.average.units":unit,
            "visit.metric":visitmetric,"effectiveness.average":q,"visit.value":"4"}
class PageSourceSamePlantPairing(unittest.TestCase):
    def test_different_studies_and_plants_never_pair(self):
        a=[row("a","Plant1","2020","moth"),row("b","Plant1","2020","fly"),
           row("a","Plant2","2020","fly")]
        self.assertEqual(compare(a,QKEY,True)["original_moth_fly_matched_cells"],0)
    def test_different_pollen_vs_seed_outcome_never_pair(self):
        a=[row("a","P","2020","moth",metric="pollen"),
           row("a","P","2020","fly",metric="seeds")]
        self.assertEqual(compare(a,PKEY,False)["original_moth_fly_matched_cells"],1)
        self.assertEqual(compare(a,QKEY,True)["original_moth_fly_matched_cells"],0)
    def test_same_study_plant_year_metric_actual_ratio_only(self):
        a=[row("a","P","2020","moth",q="10"),
           row("a","P","2020","fly",q="5")]
        r=compare(a,VKEY,True)
        self.assertEqual(r["original_moth_fly_matched_cells"],1)
        self.assertAlmostEqual(r["source_comparisons"][0]["source_fly_moth_mean_SVE_ratio"],.5)
    def test_visit_denominators_are_separate_gate(self):
        a=[row("a","P","2020","moth",visitmetric="visits/min"),
           row("a","P","2020","fly",visitmetric="visits/hour")]
        self.assertEqual(compare(a,QKEY,True)["original_moth_fly_matched_cells"],1)
        self.assertEqual(compare(a,VKEY,True)["original_moth_fly_matched_cells"],0)
if __name__=="__main__":unittest.main()
