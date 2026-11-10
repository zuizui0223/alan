"""Only synthetic transport/schema tests: no invented NPMS or GloPL results."""
import tempfile, unittest, zipfile
from pathlib import Path
from scripts.acquisition import alan_npms_glopl_source_overlap_preflight_20261110 as m


class TestSourceRequirements(unittest.TestCase):
    def test_real_like_npms_zip_csv_headers_are_inspected(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/"npms.zip"
            with zipfile.ZipFile(f,"w") as z:
                z.writestr("observations.csv","plot_id,taxon,year,cover\np1,Silene nutans,2021,1\n")
                z.writestr("metadata.txt","original file note")
            x=m.npms_inventory(f)
            self.assertEqual(x["csv_file_count"],1)
            self.assertEqual(x["members"][0]["parsed_rows"],1)
            self.assertIn("plot_id",x["members"][0]["columns"])

    def test_non_zip_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/"fake.zip"
            f.write_text("not a real data package")
            with self.assertRaisesRegex(ValueError,"not an actual ZIP"):
                m.npms_inventory(f)

    def test_effect_size_missing_is_unknown_not_zero(self):
        self.assertEqual(m.NPMS_DOI,"10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8")
        self.assertEqual(m.GLOPL_DOI,"10.5061/dryad.dt437")
        self.assertEqual(len(m.GLOPL_SHA),64)


if __name__=="__main__":
    unittest.main()
