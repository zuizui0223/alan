"""Original plant×visitor SVD functionality, false-absence and effective diversity gates."""
import unittest
from scripts.analysis.alan_ballantyne2015_heathland_visit_pollen_effectiveness_20261118 import (
    groups,plant_summary,effective_number,PLANTS
)
class TestHeathPollenVists(unittest.TestCase):
    def test_original_bees_flies_and_no_unobserved_moth(self):
        self.assertEqual(groups("Apis_mellifera"),"bees")
        self.assertEqual(groups("Episyrphus"),"flies")
        with self.assertRaises(ValueError):
            groups("imagined_night_moth")

    def test_same_plant_functional_vs_numeric_shares_not_equal(self):
        # SVD pollen delivery, rather than visitation counts, determines
        # function-weighted network links. Two 10-grain bees, two 1-grain flies.
        rows=[]
        for species,pollen in (("Apis_mellifera",10),("Apis_mellifera",10),
                               ("Episyrphus",1),("Episyrphus",1)):
            rows.append({"plant":"Erica_tetralix","date":"2013-01-01",
                         "taxon":species,"guild":groups(species),
                         "conspecific_pollen":pollen})
        x=plant_summary(rows,"Erica_tetralix")
        self.assertAlmostEqual(x["groups"]["flies"]["visit_share"],0.5)
        self.assertAlmostEqual(x["groups"]["flies"]["pollen_deposition_share"],2/22)
        self.assertAlmostEqual(x["fly_vs_bee_same_plant"]["mean_fly_per_visit_pollen_relative_to_bee"],0.1)

    def test_no_fly_visits_cannot_generate_true_zero_efficiency(self):
        rows=[{"plant":"Erica_tetralix","date":"day1",
               "taxon":"Apis_mellifera","guild":"bees","conspecific_pollen":6}]
        x=plant_summary(rows,"Erica_tetralix")
        self.assertIsNone(x["fly_vs_bee_same_plant"]["mean_fly_per_visit_pollen_relative_to_bee"])
        self.assertTrue(x["groups"]["flies"]["source_absence_in_nature_NOT_inferred"])

    def test_hill_effective_number(self):
        self.assertAlmostEqual(effective_number([10,10]),2.)
        self.assertAlmostEqual(effective_number([10,0]),1.)
        self.assertIsNone(effective_number([0,0]))

if __name__=="__main__":
    unittest.main()
