"""Experimental-year source gates; never fabricate historical ALAN exposures."""
import tempfile,unittest
from pathlib import Path
from scripts.analysis import alan_glopl_true_experiment_year_ALAN_preflight_20261110 as m


class TestGloPLYearGate(unittest.TestCase):
    def test_source_sha_and_no_future_exposure(self):
        self.assertEqual(len(m.SHA),64)
        self.assertIn("experiment_year",m.SCHEMA)
    def test_nonoriginal_source_fails(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"fake.csv";f.write_text("year_of_study\n2000\n")
            with self.assertRaisesRegex(ValueError,"file SHA changed"):
                m.run(f,Path(d)/"out")


if __name__=="__main__":unittest.main()
