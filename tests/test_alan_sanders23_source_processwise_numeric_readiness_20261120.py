import copy
import unittest
from collections import Counter

from scripts.analysis.alan_sanders23_source_processwise_numeric_readiness_20261120 import (
    process_quantitative_sufficiency, SOURCE_SHA, SCHEMA,
)

FAMILIES={
    "trophic_predation_parasitism":7,
    "aquatic_behaviour_not_emergence":5,
    "aerial_light_attraction_activity":5,
    "other_noncomparable_behaviour_reproduction":5,
    "aquatic_drift_source_NOT_equivalent_to_emergence":1,
}


def fake_original():
    rows=[]
    registry_entries=[]
    i=0
    for name,n in FAMILIES.items():
        for _ in range(n):
            sid="Study"+str(i)
            dois=["10.7777/source"+str(i)]
            if i==0:
                dois.append("10.7777/source0alt")
            registry_entries.append({
                "original_source_DOIs":dois,
                "preliminary_outcome_family_for_screen_only":name,
                "preliminary_original_habitat":"terrestrial_or_aquatic_unverified",
                "final_PRISMA_fulltext_decision":"NOT_SCREENED",
                "primary_illumination_programme_id_verified":False,
                "original_research_light_LRR_SE_verified":False
            })
            for di in dois:
                rows.append({
                    "Study":sid,"DOI":di,
                    "Title":"Source research "+sid,
                    "Delete":"Keep",
                    "Field.Lab":"Laboratory" if i==2 else "Field",
                    "Phylum":"Arthropoda","Class":"Insecta","Order":"Araneae",
                    "Sub.Category":"population abundance",
                    "Measure.description":"spider abundance at source night",
                    "Control_N":"6","Experimental_N":"6",
                    "Control_Mean":"1.4" if i!=1 else "0",
                    "Experimental_Mean":"1.9",
                    "Control_SD":"0.4","Experimental_SD":"0.3",
                })
            i+=1
    rows.append({
        "Study":"PlantOnly","DOI":"10.7777/plantA","Title":"Flowering study",
        "Delete":"Keep","Field.Lab":"Field",
        "Phylum":"Tracheophyta","Class":"Magnoliopsida","Order":"Lamiales",
        "Sub.Category":"plant reproduction","Measure.description":"seed count",
        "Control_N":"6","Experimental_N":"6","Control_Mean":"2.0",
        "Experimental_Mean":"1.7","Control_SD":"0.3","Experimental_SD":"0.4"
    })
    manual={
        "schema":"alan_sanders126_author_source_23_stageA_title_primary_process_review_v1",
        "source_original_CSV_SHA256":SOURCE_SHA,
        "distinct_source_study_records":23,
        "includes_title_screen_only_not_systematic_search":True,
        "no_automated_habitat_determination_from_plant_or_insect_taxonomy":True,
        "source_study_entries":registry_entries,
        "title_screen_only_families_counts_not_meta_ready":FAMILIES,
        "strict_aquatic_emergence_moderator_meta_eligible_k":"UNKNOWN",
        "strict_predator_abundance_moderator_meta_eligible_k":"UNKNOWN",
    }
    return rows,{"original_CSV_SHA256":SOURCE_SHA},manual


class SandersOriginalStageProcessNumericReadinessTests(unittest.TestCase):
    def test_source_process_classes_are_study_level_not_armlevel_meta_k(self):
        rows,base,reg=fake_original()
        x=process_quantitative_sufficiency(rows,base,reg,strict=False)
        self.assertEqual(x["schema"],SCHEMA)
        a=x["stageA_broad_original_source_sufficiency_totals"]
        self.assertEqual(a["source_original_Keep_StageA_Study_IDs"],23)
        self.assertEqual(a["source_original_Keep_StageA_broad_effects"],24)
        self.assertEqual(a["source_stageA_positive_original_arm_means_SD_n_rows"],23)
        self.assertEqual(a["source_stageA_Field_numerically_arithmetically_extractable_effect_rows"],22)
        self.assertEqual(x["stageA_original_23_papers_by_primary_biological_process"]["aquatic_drift_source_NOT_equivalent_to_emergence"]["source_Study_ID_clusters_from_original_author_keep"],1)
        self.assertEqual(x["stageB_original_authors_candidate_source_numeric_evidence"]["B_plant_seed_fruit_pollen_reproduction_CANDIDATE"]["source_original_Keep_effect_rows"],1)
        self.assertIsNone(x["confirmed_A_water_x_light_meta_regression_beta"])
        self.assertIsNone(x["confirmed_StageA_aquatic_emergence_and_predator_matched_effect_with_TRUE_independent_lamp_SE_k"])
        self.assertNotIn("10.7777",str(x))
        self.assertNotIn("Study1",str(x))

    def test_do_not_double_count_duplicate_author_study_dois_or_pseudoindependent_arms(self):
        rows,base,reg=fake_original()
        reg["source_study_entries"][0]["original_source_DOIs"]=["10.7777/source0"]
        with self.assertRaisesRegex(ValueError,"DOI list"):
            process_quantitative_sufficiency(rows,base,reg,strict=False)

    def test_wrong_source_digest_fails_before_biological_labeling(self):
        rows,base,reg=fake_original()
        base["original_CSV_SHA256"]="wrong"
        with self.assertRaisesRegex(ValueError,"SHA"):
            process_quantitative_sufficiency(rows,base,reg,strict=False)

    def test_do_not_convert_numerical_fields_to_truth_of_random_lighting(self):
        rows,base,reg=fake_original()
        x=process_quantitative_sufficiency(rows,base,reg,strict=False)
        for sub in x["stageA_original_23_papers_by_primary_biological_process"].values():
            self.assertTrue(sub["not_a_confirmed_independent_field_illumination_programme_count"])
            self.assertTrue(sub["not_an_original_verified_same_endpoint_meta_sample"])
        self.assertTrue(x["no_same_study_multiarm_shared_control_covariance_derived"])
        self.assertTrue(x["original_raw_meta_effect_rows_or_Study_ID_exported"])

if __name__=="__main__":
    unittest.main()
