import unittest
from scripts.analysis.alan_seymoure2023_studycluster_guild_time_space_screen_20261119 import (
    fisher_two_sided, paired_binomial_two_sided, source_role_audit
)


def row(i, typ, flags, study):
    dec, inc, tem, spa = flags
    return {
        "_original_excel_row":i,
        "interaction_type":typ,
        "method":"",
        "Consumer_taxon":"animal",
        "Resource_taxon":"resource",
        "interaction strength (IS) reduction":str(dec),
        "IS enhancement":str(inc),
        "strength sum":str(dec+inc),
        "temporal mismatch induced":str(tem),
        "spatial mismatch induced":str(spa),
        "effect sum":str(dec+inc+tem+spa),
        "Authors":study,"Article Title":"Study "+study,
        "Publication Year":"2021","DOI":"10.1234/"+study,"DOI Link":"",
    }


def original_like():
    x=[
        row(2,"predator-prey",(0,1,0,1),"P1"),
        row(3,"predator-prey",(1,0,0,0),"P2"),
        row(4,"herbivory",(1,0,0,1),"H1"),
        row(5,"herbivory",(0,1,1,0),"H2"),
        row(6,"pollination",(0,1,1,1),"P1"),
    ]
    foot=row(7,"",(0,0,0,0),"F")
    foot.update({"interaction_type":"","Authors":"","Article Title":"",
        "Publication Year":"","DOI":"",
        "interaction strength (IS) reduction":"29",
        "IS enhancement":"28","strength sum":"57",
        "temporal mismatch induced":"7","spatial mismatch induced":"27",
        "effect sum":""})
    return x+[foot]


class OriginalStudyClusterGuildTimeSpaceTests(unittest.TestCase):
    def test_probability_ordered_fisher_known_source_sign_table(self):
        self.assertAlmostEqual(fisher_two_sided(19,13,2,7),0.06695216451313872)
        self.assertAlmostEqual(fisher_two_sided(1,0,0,1),1.0)
        with self.assertRaises(ValueError):
            fisher_two_sided(-1,1,1,1)

    def test_exact_mcnemar_for_discordant_coding_only(self):
        self.assertAlmostEqual(paired_binomial_two_sided(24,4),0.00017999112606048584)
        self.assertEqual(paired_binomial_two_sided(0,0),1.0)
        self.assertEqual(paired_binomial_two_sided(1,0),1.0)

    def test_single_publication_multiple_interactions_cannot_be_pooled(self):
        ans=source_role_audit(original_like(),strict=False)
        self.assertEqual(ans["original_source_interactions"],5)
        self.assertEqual(ans["independent_publications_in_review"],4)
        self.assertEqual(ans["predator_prey_versus_herbivory_coding"]["predator_prey_publication_overlap_herbivory"],0)
        self.assertTrue(ans["inference_gate"]["not_a_random_effects_effect_size_meta_analysis"])
        self.assertTrue(ans["inference_gate"]["no_original_study_row_or_DOI_export"])
        self.assertNotIn("10.1234",str(ans))
        self.assertEqual(ans["within_publication_original_spatiotemporal_coding"]["both"],1)

    def test_cross_guild_overlap_requires_joint_cluster_model(self):
        x=original_like()
        x[3]["DOI"]="10.1234/P1"
        with self.assertRaisesRegex(ValueError,"share publications"):
            source_role_audit(x,strict=False)


if __name__=="__main__":
    unittest.main()
