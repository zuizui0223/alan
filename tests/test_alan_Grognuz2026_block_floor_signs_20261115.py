import unittest
from scripts.analysis import alan_Grognuz2026_block_floor_signs_20261115 as m
class TestGrognuzFloor(unittest.TestCase):
    def test_source_structure(self):
        d=m.load()
        self.assertEqual(len(d),28)
        self.assertEqual(sum(x["block_id"].startswith("T") for x in d),20)
        self.assertEqual(sum(x["block_id"].startswith("C") for x in d),8)
    def test_recompute_known_original_block_balanced_means(self):
        d=m.load();a=m.analyze(d)
        t=a["Trifolium_pratense"];c=a["Centaurea_jacea"]
        self.assertAlmostEqual(t["mean_dark_seed_per_target"],4.925,places=4)
        self.assertAlmostEqual(c["mean_dark_seed_per_target"],1.47916667,places=4)
        self.assertAlmostEqual(t["mean_lit_seed_per_target"],3.51875,places=4)
    def test_nonnegativity_of_dark_zero(self):
        d=m.load()
        for z in d:
            if z["mean_dark_no_slug"]==0 and z["mean_dark_slug"]==0:
                self.assertTrue((z["mean_lit_no_slug"]+z["mean_lit_slug"])>=0)
    def test_author_repeated_lights_never_28_independent_experiments(self):
        a=m.analyze(m.load())
        self.assertTrue(a["source_is_derived_block_means_NOT_original_pot_rows"])
        self.assertNotIn("clustered_p_value",a["Trifolium_pratense"])
if __name__=="__main__":unittest.main()
