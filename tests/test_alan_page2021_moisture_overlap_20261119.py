"""No wet-plant moderator unless original studies independently support both traits."""
import unittest
from scripts.analysis.alan_page2021_moisture_function_overlap_20261119 import audit

def n(study,plant,missing=False):
    return {
        "source_study":study,"source_plant":plant,"source_top_visitor_group":"bee",
        "requires_more_than_equal_visits":True,
        "required_new_visits_per_lost_visit":2.,
        "taxa":5
    }

class WetLinkCoverage(unittest.TestCase):
    def test_two_studies_do_not_become_meta_evidence(self):
        records=[n("StudyA","Wet_a"), n("StudyA","Dry_b"),
                 n("StudyB","Wet_a"), n("StudyB","Dry_b"),
                 n("StudyC","No_record")]
        moisture={"Wet a":8.,"Dry b":4.}
        x=audit(records,moisture)
        self.assertEqual(x["matched_study_networks"],4)
        self.assertEqual(x["matched_wet_networks"],2)
        self.assertEqual(x["source_studies_with_both_wet_and_dry"],2)
        self.assertEqual(x["status"],"SOURCE_WITHIN_STUDY_MOISTURE_SUPPORT_STOP")
        self.assertIsNone(x["functional_redundancy_moisture_meta_moderator"])

    def test_missing_moisture_remains_missing(self):
        x=audit([n("Study","Unknown_plant")],{})
        self.assertEqual(x["matched_dry_networks"],0)
        self.assertEqual(x["matched_wet_networks"],0)
        self.assertTrue(x["missing_EIV_does_not_mean_dry"])

if __name__ == "__main__":
    unittest.main()
