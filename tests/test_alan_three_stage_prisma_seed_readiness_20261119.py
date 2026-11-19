import copy
import json
import unittest
from pathlib import Path
from scripts.analysis.alan_three_stage_prisma_seed_readiness_20261119 import audit

SOURCE=Path("config/alan_prisma_three_stage_seed_census_20261119.json")


class ThreeStagePRISMASourceReadinessTests(unittest.TestCase):
    def setUp(self):
        self.data=json.loads(SOURCE.read_text())

    def test_purposeful_literature_seed_is_not_prisma_meta_inclusion(self):
        d=audit(self.data)
        self.assertEqual(d["S1_known_intentionally_seeded_primary_paper_records"],11)
        self.assertEqual(d["S2_known_intentionally_seeded_primary_paper_records"],11)
        self.assertIsNone(d["S1_primary_META_ready_k"])
        self.assertIsNone(d["S2_primary_META_ready_k"])
        self.assertIsNone(d["PRISMA_records_identified_from_database"])
        self.assertIsNone(d["PRISMA_S1_included"])
        self.assertIsNone(d["PRISMA_S2_included"])
        self.assertFalse(d["WoS_Scopus_original_searches_executed"])
        self.assertTrue(d["no_lnRR_meta_or_causal_mediation_coefficient_estimated"])

    def test_preregistration_not_fabricated_after_prior_data_inspection(self):
        r=copy.deepcopy(self.data)
        r["review_status"]="PROSPECTIVELY_PREREGISTERED_WITH_OSF"
        with self.assertRaisesRegex(ValueError,"newly preregistered"):
            audit(r)

    def test_uvb_and_prior_review_endpoints_not_poolable_as_light_studies(self):
        r=copy.deepcopy(self.data)
        reviews=r["baseline_reviews"]
        guenat=next(x for x in reviews if x["doi"]=="10.1002/ece3.9990")
        guenat["source_ALAN_nonUVB_light_candidate_articles"]=22
        with self.assertRaisesRegex(ValueError,"UV-B"):
            audit(r)
        r=copy.deepcopy(self.data)
        r["seed_candidates"][0]["doi"]="10.1002/ece3.9990"
        with self.assertRaisesRegex(ValueError,"review paper"):
            audit(r)

    def test_light_experiment_not_multiple_papers_or_rows(self):
        r=copy.deepcopy(self.data)
        r["source_reference_id_concordance"]["penryn_2015_2018_2022_same_54_troughs"]=False
        with self.assertRaisesRegex(ValueError,"source independence"):
            audit(r)
        r=copy.deepcopy(self.data)
        r["seed_candidates"].append(copy.deepcopy(r["seed_candidates"][0]))
        with self.assertRaisesRegex(ValueError,"Duplicate stage"):
            audit(r)

    def test_westhavelland_same_physical_site_pair_is_not_two_independent_ALAN_trials(self):
        d=audit(self.data)
        self.assertTrue(d["Westhavelland_2017_2025_share_one_lit_one_dark_field_pair"])
        self.assertTrue(d["source_A_11_seed_papers_not_11_independent_light_programmes"])
        x=copy.deepcopy(self.data)
        x["seed_candidates"]=[r for r in x["seed_candidates"] if r["study_programme_id"]!="Manfrin2017Westhavelland"]
        x["counts"]["S1_purposeful_candidate_paper_records_not_search_inclusions"]=10
        with self.assertRaisesRegex(ValueError,"11 A/11 B"):
            audit(x)
        x=copy.deepcopy(self.data)
        x["source_reference_id_concordance"]["westhavelland_2017_2025_same_1lit_1dark_sitepair"]=False
        with self.assertRaisesRegex(ValueError,"source independence"):
            audit(x)

    def test_publish_bias_preserve_positive_and_negative_counterexamples(self):
        r=copy.deepcopy(self.data)
        r["seed_candidates"]=[x for x in r["seed_candidates"] if x["study_programme_id"]!="Macgregor2019"]
        # Removing any one pre-identified candidate is source registry drift.
        with self.assertRaisesRegex(ValueError,"11 A/11 B"):
            audit(r)
        r=copy.deepcopy(self.data)
        r["seed_candidates"][0]["screening_status"]="INCLUDED_IN_META"
        with self.assertRaisesRegex(ValueError,"systematically included"):
            audit(r)


if __name__=="__main__":
    unittest.main()
