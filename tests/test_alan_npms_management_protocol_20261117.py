"""No fake absent species; source-reported management labels must truly match every year."""
import unittest
import pandas as pd
from scripts.analysis import alan_npms_management_protocol_20261117 as a

class NPMSManagementSourceGuardTests(unittest.TestCase):
    def test_stable_only_when_all_source_years_labeled(self):
        frame=pd.DataFrame([
            {"location_id":"p1","Management":"Mow"},
            {"location_id":"p1","Management":"Mow"},
            {"location_id":"p2","Management":"Grazing"},
            {"location_id":"p2","Management":None},
            {"location_id":"p3","Management":"Grazing"},
            {"location_id":"p3","Management":"Mow"},
            {"location_id":"p4","Management":None},
            {"location_id":"p4","Management":None},
        ])
        ids,stats=a.stable_plot_labels(frame,"Management")
        self.assertEqual(ids,{"p1"})
        self.assertEqual(stats["plots_stable_all_years_labeled"],1)
        self.assertEqual(stats["plots_with_any_missing_label"],2)
        self.assertEqual(stats["plots_with_two_or_more_observed_labels"],1)

    def test_single_year_does_not_prove_stability(self):
        frame=pd.DataFrame([{"location_id":"p1","Grazing":"Yes"}])
        ids,_=a.stable_plot_labels(frame,"Grazing")
        self.assertEqual(ids,set())

if __name__=="__main__":
    unittest.main()
