"""Mathematical/identification guards for Li 2026 original three-hayfield design."""
import unittest
from scripts.analysis.alan_Li2026_real_geometry_seed_hayfield_source_20261114 import (
    summary, contrast, field_comparisons, data_quality, complete_plant_field_checks
)

def controlled_fixture():
    out=[]
    for field in (1,2,3):
        for treatment in ("C","P","L"):
            for plant in (1,2):
                success=int(treatment=="L" or (treatment=="P" and field==1))
                out.append({"species":"ML","round":"1","Label":f"{field}-{treatment}-{plant}",
                            "replicate":str(field),"Treatment":treatment,"id":"NA",
                            "pollinated":str(success),"Seeds":str(success*(field+1)),
                            "Fruit":str(success),"Starting":"1","No.uf":"1"})
    return out

class TestOriginalGeometrySource(unittest.TestCase):
    def test_expected_plant_measures_without_pseudoreplication(self):
        x=controlled_fixture()
        z=summary(x)
        self.assertEqual(z["nonmissing_seedset_source_plants"],18)
        self.assertEqual(z["has_at_least_one_seed_plants"],8)
        self.assertEqual(contrast(x,"L")["raw_seed_presence_risk_difference"],1.0)
        f=field_comparisons(x,"L")
        self.assertEqual(f["hayfields_positive_difference"],3)
        self.assertEqual(f["independent_hayfields"],3)
        self.assertAlmostEqual(f["two_sided_sign_fraction_descriptive_not_original_assignment_p"],0.25)
        self.assertEqual(len(f["leave_one_whole_hayfield_out_pooled_plant_risk_differences"]),3)

    def test_unknown_not_zero_and_label_collision_reported(self):
        x=controlled_fixture()
        x[0]["pollinated"]="NA"
        x[0]["Seeds"]="NA"
        x[0]["Fruit"]="NA"
        q=summary(x)
        self.assertEqual(q["nonmissing_seedset_source_plants"],17)
        self.assertEqual(data_quality(x)["missing_binary_seedset_rows"],1)
        x[1]["Label"]=x[2]["Label"]
        # source plant label ID is only descriptive; no forced deduplication.
        self.assertEqual(summary(x)["all_original_source_rows"],18)

    def test_no_missing_species_field_silently_filled(self):
        x=controlled_fixture()
        z=complete_plant_field_checks(x)
        self.assertEqual(z["supported_strata"],3)
        self.assertEqual(z["positive_supported_strata"],3)
        self.assertEqual(z["negative_supported_strata"],0)
        self.assertEqual(len(z["strata"]),9)
        self.assertTrue(any(not t["complete_both_light_treatments"] for t in z["strata"]))

if __name__=="__main__":
    unittest.main()
