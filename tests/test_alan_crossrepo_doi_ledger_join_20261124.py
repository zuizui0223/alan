"""The cross-repo DOI join must never claim unverified papers are eligible."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analysis/alan_crossrepo_doi_ledger_join_20261124.py"
spec=importlib.util.spec_from_file_location("doi_join",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def write_csv(path,header,rows):
    with path.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=header)
        w.writeheader()
        w.writerows(rows)


class TestCrossRepoDOIJoin(unittest.TestCase):
    def test_doi_normalization_and_prefixes(self):
        self.assertEqual(mod.normalize_doi("HTTPS://DOI.ORG/10.1093/JPE/RTAG160."),
                         "10.1093/jpe/rtag160")
        self.assertEqual(mod.normalize_doi("doi:10.1038/nature23288;"),"10.1038/nature23288")
        self.assertFalse(mod.DOI_RE.fullmatch(mod.normalize_doi("MED:42184954")))

    def test_author_source_roster_contains_prior_alan_repro(self):
        status,meta=mod.load_existing()
        self.assertIn("ALAN_REPRODUCTIVE_ADMISSION",status["10.1093/jpe/rtag160"])
        self.assertIn("ALAN_MAIN_STUDY_REGISTER",status["10.1038/nature23288"])
        self.assertIn("ALAN_PENDING_55",status["10.1016/j.marenvres.2021.105285"])
        self.assertGreater(len(status),100)

    def test_676_source_ids_do_not_turn_into_676_studies(self):
        with tempfile.TemporaryDirectory() as dirname:
            home=Path(dirname)
            a=home/"inventory"
            a.mkdir()
            names=["10.1093/jpe/rtag160","10.9999/unknown-alan",
                   "10.9999/unknown-context","10.9999/unknown-unrelated"]
            names += [f"10.9999/unique-{i:04d}" for i in range(672)]
            self.assertEqual(len(set(names)),676)
            occ=[]
            for i,doi in enumerate(names):
                repo=("zuizui0223/alan" if i in (0,1)
                      else "zuizui0223/boids" if i==2
                      else "zuizui0223/fd")
                path="docs/alan_original_source.md" if i==2 else "docs/OTHER.md"
                occ.append({"repo":repo,"path":path,"git_blob_sha":"a"*40,
                            "line":str(i+1),"doi":doi})
            write_csv(a/"doi_occurrence_sources.csv",list(occ[0]),occ)
            rows=[{"repository":f"zuizui0223/repo_{i:02d}"} for i in range(62)]
            write_csv(a/"repository_coverage.csv",list(rows[0]),rows)
            info={**mod.EXPECTED_BASELINE,"DOI_source_occurrences_distinct_by_file":676,
                  "all_references_full_text_adjudicated":False,
                  "all_default_branch_text_files_scanned":False}
            (a/"summary.json").write_text(json.dumps(info),encoding="utf-8")
            result,report=mod.make_crosswalk(a)
            self.assertEqual(len(result),676)
            self.assertEqual(report["new_meta_programmes_added"],0)
            self.assertEqual(report["accepted_effects_added"],0)
            d={x["doi_string"]:x for x in result}
            self.assertEqual(d["10.1093/jpe/rtag160"]["triage_priority"],
                             "KNOWN_REFERENCE_NOT_NEW_STUDY")
            self.assertEqual(d["10.9999/unknown-alan"]["triage_priority"],
                             "P1_EXACT_UNMATCHED_WITH_ALAN_IWE_SOURCE")
            self.assertEqual(d["10.9999/unknown-context"]["triage_priority"],
                             "P2_EXACT_UNMATCHED_IN_RELEVANT_DOC_PATH")
            self.assertEqual(d["10.9999/unknown-unrelated"]["triage_priority"],
                             "P3_EXACT_UNMATCHED_OTHER_CONTEXT")
            self.assertFalse(d["10.9999/unknown-alan"]["novel_primary_article_confirmed"])
            self.assertFalse(d["10.9999/unknown-alan"]["meta_eligible_effect_confirmed"])

            # No silent source-coverage inflation on a changed artifact.
            info["repositories_with_full_tree"]=61
            (a/"summary.json").write_text(json.dumps(info),encoding="utf-8")
            with self.assertRaises(ValueError):mod.make_crosswalk(a)

if __name__=="__main__":
    unittest.main()
