"""Original full Lu network source: no silent missing-light reconstruction."""
import unittest

from scripts.analysis.alan_lu2026_original_fullnetwork_sameflower_20261118 import (
  parse_source,compare_views,SOURCE_ROWS,COMPLETE_RECORDS
)

class TestOriginalFullNetworkSource(unittest.TestCase):
    @staticmethod
    def fixture(two_invalid=False):
        names=["time","sampling_time","treatment","block","repeat","plant",
               "all","Lepidoptera","Diptera","Hymenoptera","Coleoptera","Others"]
        out=[names]
        for j in range(COMPLETE_RECORDS):
            out.append(["1","Nocturnal","Lit" if j%2 else "Unlit",
                        str(j%5+1),"1","Plant_A",1,1,0,0,0,0])
        out.append(["1","","","1","1","Plant_A",1,1,0,0,0,0])
        if two_invalid:out[-2][1]=""
        self_count=len(out)-1
        if self_count!=SOURCE_ROWS:raise AssertionError(self_count)
        return out

    def test_single_original_unlabelled_source_record_never_forward_filled(self):
        valid,support=parse_source(self.fixture())
        self.assertEqual(len(valid),COMPLETE_RECORDS)
        self.assertEqual(support["original_rows_with_unreported_light_or_diel"],1)
        self.assertEqual(support["original_distinct_plant_labels"],1)
        self.assertTrue(support["original_unlabelled_plant_not_assumed_0_visits"])

    def test_two_missing_source_light_labels_stop_without_guessing(self):
        with self.assertRaises(ValueError):
            parse_source(self.fixture(two_invalid=True))

    def test_three_views_dont_automatically_become_one_source_cohort(self):
        valid,_=parse_source(self.fixture())
        visit=[["sampling_time","abundance","D_abundance","L_abundance",
                "H_abundance","C_abundance","O_abundance"],
               ["Nocturnal",100,1,99,0,0,0]]
        focal=[["sampling_time","all","Diptera","Lepidoptera",
                "Hymenoptera","Coleoptera","Others"],
               ["Nocturnal",10,1,9,0,0,0]]
        result=compare_views(valid,visit,focal)
        self.assertFalse(result["source_three_views_exact_same_visitation_totals"])
        self.assertTrue(result["do_not_sum_across_three_views"])

if __name__=="__main__":
    unittest.main()
