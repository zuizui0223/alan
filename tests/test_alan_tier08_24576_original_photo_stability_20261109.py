"""Synthetic parameter checks only, not observational ecological evidence."""
import unittest
from scripts.analysis import alan_tier08_24576_original_photo_stability_20261109 as m

class PrecommittedStabilityGatesTests(unittest.TestCase):
    def stage(self,beta,low,high,islands=121,cells=3300,equal=None):
        z={"island_equal_area_beta":beta,"island_equal_area_CI_low":low,
           "island_equal_area_CI_high":high,"multisite_islands":islands,
           "multisite_site_years":cells}
        if equal is not None:z["equal_source_area_beta"]=equal
        return z
    def test_stable_null_allowed_without_p_value_significance(self):
        x=[self.stage(.5,-.4,1.4),self.stage(.4,-.5,1.3),
           self.stage(.25,-.55,1.05),self.stage(.2,-.75,1.15,equal=.1)]
        y=m.convergence_gate(x)
        self.assertTrue(y["stable"])
        self.assertTrue(y["p_value_does_not_control_stopping"])
        self.assertFalse(y["needs_additional_genuine_photo_waves"])
    def test_significant_result_not_criterion_and_no_support_means_continue(self):
        x=[self.stage(4,3,5),self.stage(4.3,3.3,5.3),
           self.stage(5,4,6),self.stage(6,5,7,islands=40,cells=800,equal=6)]
        y=m.convergence_gate(x)
        self.assertFalse(y["stable"])
        self.assertTrue(y["needs_additional_genuine_photo_waves"])
        self.assertFalse(y["gates"]["at_least_120_spatially_estimable_islands"])
    def test_multiple_wave_changes_required(self):
        y=m.convergence_gate([self.stage(.5,-1,2)])
        self.assertFalse(y["stable"])
        self.assertEqual(y["reason"],"need baseline and three signed 8192-photo waves")
if __name__=="__main__":unittest.main()
