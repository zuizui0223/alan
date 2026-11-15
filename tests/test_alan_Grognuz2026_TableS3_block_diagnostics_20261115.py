import importlib.util
import math
from pathlib import Path
import unittest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_Grognuz2026_TableS3_block_diagnostics_20261115.py"
s=importlib.util.spec_from_file_location("grognuz_s3",SRC)
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

class TestTableS3Derived(unittest.TestCase):
    def test_n_and_missing(self):
        groups=m.load_source()
        self.assertEqual((len(groups["Trifolium"]),len(groups["Centaurea"])),(20,8))
        self.assertEqual(sum(sum(x["n"]) for g in groups.values() for x in g),425)
        self.assertEqual(28*16-425,23)
    def test_reproduces_original_published_means(self):
        groups=m.load_source()
        # This check fails in load_source if published means disagree.
        for name,g in groups.items():
            vals,n=m.extract(g)
            self.assertEqual(sum(n),311 if name=="Trifolium" else 114)
            self.assertEqual(len(vals),4)
    def test_critical_source_effect_splits(self):
        g=m.load_source()
        t=m.effects(m.extract(g["Trifolium"])[0])
        c=m.effects(m.extract(g["Centaurea"])[0])
        self.assertAlmostEqual(t["light_dark_ratio"],0.714467,places=5)
        self.assertAlmostEqual(c["light_dark_ratio"],0.714789,places=5)
        self.assertLess(c["no_slug_ratio"],1)
        self.assertGreater(c["slug_ratio"],1)
    def test_boot_resample_deterministic(self):
        group=m.load_source()["Centaurea"]
        self.assertEqual(m.block_boot(group,42,90),m.block_boot(group,42,90))
    def test_exploratory_only_cluster_limits(self):
        groups=m.load_source()
        for name,g in groups.items():
            x=m.one_species(name,g,100)
            self.assertTrue(x["source_has_plant_reuse_across_weeks"])
            self.assertIsNone(x["actual_independent_lamp_installations"])
            self.assertEqual(x["n_target_pots_total"],len(g)*16)

if __name__=="__main__":
    unittest.main()
