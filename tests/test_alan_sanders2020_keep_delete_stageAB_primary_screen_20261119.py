import copy
import unittest

from scripts.analysis.alan_sanders2020_keep_delete_stageAB_primary_screen_20261119 import (
    audit_keep_delete,doi_norm,source_label_stage,original_arm_summaries_positive,
    SOURCE_SHA256,
)

def row(study,typ,status="Keep",design="Field",mean="5"):
    x={
       "Study":study,"ID":study+"effect1","DOI":"https://doi.org/10.1234/"+study,
       "Title":"Scientific published source study "+study,
       "Phylum":"Arthropoda" if typ=="predator" else "angiosperms",
       "Class":"Arachnida" if typ=="predator" else "monocots",
       "Order":"Araneae" if typ=="predator" else "Poales",
       "Species":"test species","Measure.description":
          "spider abundance" if typ=="predator" else "seed output",
       "Sub.Category":"Abundance" if typ=="predator" else "Reproductive.output",
       "Delete":status,"Field.Lab":design,
       "Control_N":"4","Experimental_N":"4",
       "Control_Mean":"3","Experimental_Mean":mean,
       "Control_SD":"1","Experimental_SD":"1"}
    return x

class SandersOriginalSourceKeepExclusionTests(unittest.TestCase):
    def setUp(self):
        self.base={"original_CSV_SHA256":SOURCE_SHA256}
        self.rows=[row("p1","predator"),row("p1","predator",status="Delete"),
                   row("p2","predator",design="Laboratory"),
                   row("p3","plant",mean="0")]
    def test_stage_keep_delete_and_independence_not_meta(self):
        x=audit_keep_delete(self.rows,self.base,strict=False)
        self.assertEqual(x["published_1304_original_rows"],4)
        a=x["stage_AB_source_screening_with_original_keep_delete_and_design"]["A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN"]
        b=x["stage_AB_source_screening_with_original_keep_delete_and_design"]["B_plant_seed_fruit_pollen_reproduction_CANDIDATE"]
        self.assertEqual(a["all_source_effect_rows_pre_source_exclusion"],3)
        self.assertEqual(a["author_original_Keep_effect_rows"],2)
        self.assertEqual(a["author_original_Delete_effect_rows"],1)
        self.assertEqual(a["author_original_Keep_unique_Study_IDs"],2)
        self.assertEqual(a["author_original_Keep_Field_effect_rows"],1)
        self.assertEqual(a["author_original_Keep_Laboratory_effect_rows"],1)
        self.assertEqual(b["author_original_Keep_positive_arm_means_SD_n_candidate_rows"],0)
        self.assertEqual(len(x["candidate_publication_bibliography_for_fulltext_screening"]["A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN"]),2)
        self.assertTrue(x["original_author_Delete_rows_not_assumed_biologically_null"])
        self.assertTrue(x["no_random_effect_meta_fitted"])
        self.assertNotIn("p1effect1",str(x))
    def test_doi_canonical_and_refuse_unverified_arms(self):
        self.assertEqual(doi_norm("https://doi.org/10.1234/X1"),"10.1234/x1")
        self.assertEqual(doi_norm("doi: 10.9999/example"),"10.9999/example")
        self.assertIsNone(doi_norm("https://example.com/example"))
        self.assertFalse(original_arm_summaries_positive(row("p2","predator",mean="-1")))
        self.assertTrue(original_arm_summaries_positive(row("p2","predator",mean="1")))
    def test_taxonomic_filter_wont_call_all_insects_predators(self):
        r=row("p2","predator")
        r["Measure.description"]="fly activity"
        r["Sub.Category"]="Activity"
        self.assertFalse(source_label_stage(r)["A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN"])
        r["Measure.description"]="spider abundance"
        self.assertTrue(source_label_stage(r)["A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN"])
    def test_do_not_claim_preregistration_from_author_delete_flag(self):
        bad=copy.deepcopy(self.base)
        bad["original_CSV_SHA256"]="incorrect"
        with self.assertRaisesRegex(ValueError,"SHA mismatch"):
            audit_keep_delete(self.rows,bad,strict=False)
        w=audit_keep_delete(self.rows,self.base,strict=False)
        self.assertTrue(all(r["not_an_included_PRISMA_study_yet"] for v in w["candidate_publication_bibliography_for_fulltext_screening"].values() for r in v))
if __name__=="__main__":
    unittest.main()
