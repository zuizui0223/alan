"""Fixed 30 DOI-source leads are not 30 new papers."""
import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
M=ROOT/"scripts/analysis/alan_crossrepo_p1_source_triage_20261124.py"
s=importlib.util.spec_from_file_location("p1_type",M)
mod=importlib.util.module_from_spec(s);s.loader.exec_module(mod)

DATA=[
"10.16904/envidat.628","10.5061/dryad.1rn8pk143",
"10.5061/dryad.8217n76","10.5061/dryad.95x69p8r8",
"10.5061/dryad.gb9534p","10.5061/dryad.gxd2547xh",
"10.5061/dryad.h7k26","10.5061/dryad.hhmgqnknt",
"10.5061/dryad.pq54gh5","10.5061/dryad.sxksn03gp",
"10.5281/zenodo.14516568","10.5281/zenodo.1743917",
"10.5281/zenodo.4036279","10.5281/zenodo.4540407",
"10.5285/7e7ce7c1-aa16-484f-86b1-90c5b088a64c",
"10.6084/m9.figshare.20304039","10.6084/m9.figshare.27273918.v1",
"10.6084/m9.figshare.c.3295991","10.6084/m9.figshare.c.3889075",
"10.6084/m9.figshare.c.6837642"
]

def row(doi):
    return {"doi_string":doi,
            "triage_priority":"P1_EXACT_UNMATCHED_WITH_ALAN_IWE_SOURCE",
            "source_repositories":"zuizui0223/alan",
            "source_file_locations":"zuizui0223/alan@abc:docs/foo.md:7"}

class TestSourceTypeReconciliation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources=[row(z) for z in (DATA+list(mod.ARTICLES))]
        cls.output,cls.meta=mod.triage(cls.sources)

    def test_count_distinguishes_dataset_deposits_from_journal_articles(self):
        self.assertEqual(len(self.output),30)
        self.assertEqual(self.meta["data_repository_identifiers"],20)
        self.assertEqual(self.meta["non_dataset_publication_identifiers"],10)
        self.assertEqual(self.meta["ALAN_biotic_manipulation_reports"],4)
        self.assertEqual(self.meta["Westhavelland_reused_reports"],2)
        self.assertEqual(self.meta["new_meta_eligible_programmes_added"],0)
        self.assertTrue(all(r["independently_adoptable_effect_added"] is False for r in self.output))

    def test_non_original_alan_not_mis_adopted(self):
        rec={r["source_doi_literal"]:r for r in self.output}
        self.assertEqual(rec["10.1098/rsbl.2019.0877"]["source_type_gate"],"NO_EXPERIMENTAL_ALAN_POLLEN_NETWORK")
        self.assertEqual(rec["10.1111/jvs.13168"]["source_type_gate"],"PLANT_EIV_TRAIT_REFERENCE")
        self.assertEqual(rec["10.1007/s00027-025-01161-7"]["known_or_candidate_physical_programme"],
                         rec["10.3389/fenvs.2017.00061"]["known_or_candidate_physical_programme"])
        self.assertEqual(rec["10.3389/fpls.2023.1263795/full"]["canonical_candidate_doi"],
                         "10.3389/fpls.2023.1263795")
        self.assertEqual(rec["10.5061/dryad.gb9534p"]["meta_admission_gate"],
                         "CROSSREFERENCE_TO_SOURCE_PROGRAMME_NOT_A_NEW_STUDY")

    def test_unclassified_or_mismatched_DOI_stops(self):
        with self.assertRaises(ValueError):
            mod.triage(self.sources+[row("10.9999/newpaper")])
        broken=[dict(x) for x in self.sources]
        broken[0]["doi_string"]="10.9999/newpaper"
        with self.assertRaises(ValueError):mod.triage(broken)
        broken=[dict(x) for x in self.sources]
        broken[0]["triage_priority"]="KNOWN_REFERENCE_NOT_NEW_STUDY"
        with self.assertRaises(ValueError):mod.triage(broken)

if __name__=="__main__":
    unittest.main()
