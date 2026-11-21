"""Vial-level bootstrap ensures no individual organism pseudo-replication."""
import random,unittest
from scripts.analysis.alan_freshwater2026_three_consumer_vial_feeding import read_group, boot_logresponse

class TestSourceVialEffect(unittest.TestCase):
    def test_logratio_known_mean(self):
        x={0.:[1.,1.,1.],1.:[2.,2.,2.],10.:[.5,.5,.5]}
        stat,b=boot_logresponse(x,random.Random(33),199)
        self.assertAlmostEqual(stat["ratio_to_zero_by_lux"]["1"],2.)
        self.assertAlmostEqual(stat["ratio_to_zero_by_lux"]["10"],.5)
        self.assertEqual(len(b["1"]),199)

    def test_unreported_not_renamed_zero(self):
        spec={"species":"Synthetic","vials_expected_per_lux":2}
        sheet=[["ALAN","replicate","Feeding rate"],["0","1","1."],
               ["0","2",""],["1","1","2."],["1","2","3."],
               ["10","1","5."],["10","2","6."]]
        values,meta=read_group(sheet,spec)
        self.assertEqual(len(values[0.]),1)
        self.assertEqual(meta["0"]["missing_feeding_vials"],1)
        self.assertEqual(meta["0"]["original_vials"],2)

    def test_source_replicate_duplicate_rejected(self):
        spec={"species":"Synthetic","vials_expected_per_lux":2}
        sheet=[["ALAN","replicate","Feeding rate"],["0","1","1."],
               ["0","1","1."],["1","1","2."],["1","2","3."],
               ["10","1","5."],["10","2","6."]]
        with self.assertRaisesRegex(ValueError,"duplicated"):
            read_group(sheet,spec)

if __name__=="__main__":unittest.main()
