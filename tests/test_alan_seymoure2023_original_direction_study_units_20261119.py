import copy
import unittest
from scripts.analysis.alan_seymoure2023_original_direction_study_units_20261119 import summarize_source


def dummy(i, typ, flags, study):
    dec,inc,tmp,spc=flags
    return {
        "_original_excel_row":i,
        "interaction_type":typ,
        "method":"field light experiment" if study=="A" else "",
        "Consumer_taxon":"insect",
        "Resource_taxon":"plant",
        "interaction strength (IS) reduction":str(dec),
        "IS enhancement":str(inc),
        "strength sum":str(dec+inc),
        "temporal mismatch induced":str(tmp),
        "spatial mismatch induced":str(spc),
        "effect sum":str(sum(flags)),
        "Authors":study,
        "Article Title":"Paper "+study,
        "Publication Year":"2022",
        "DOI":"10.1234/"+study,
        "DOI Link":"",
    }


def valid_source():
    footer=dummy(5,"",(0,0,0,0),"")
    footer.update({
       "interaction_type":"",
       "Authors":"","Article Title":"","Publication Year":"",
       "DOI":"",
       "interaction strength (IS) reduction":"29",
       "IS enhancement":"28",
       "strength sum":"57",
       "temporal mismatch induced":"7",
       "spatial mismatch induced":"27",
       "effect sum":"",
    })
    return [
        dummy(2,"predator-prey",(1,0,0,1),"A"),
        dummy(3,"pollination",(0,1,1,0),"A"),
        dummy(4,"herbivory",(0,0,0,0),"B"),
        footer,
    ]


class SeymoureOriginalSourceDirectionTests(unittest.TestCase):
    def test_study_clustered_interaction_counts(self):
        x=summarize_source(valid_source(),strict=False)
        i=x["interaction_level_original_coded_counts"]
        p=x["publication_level_DOI_aggregated_counts"]
        self.assertEqual(i["investigations"],3)
        self.assertEqual(i["coded_reduction"],1)
        self.assertEqual(i["coded_enhancement"],1)
        self.assertEqual(i["coded_temporal"],1)
        self.assertEqual(i["coded_spatial"],1)
        self.assertEqual(i["total_flagged_source_effect_entries_not_independent"],4)
        self.assertEqual(p["publications"],2)
        self.assertEqual(p["papers_with_multiple_interaction_rows"],1)
        self.assertEqual(p["coded_both_directions_in_same_article"],1)
        self.assertEqual(p["coded_any_strength"],1)
        self.assertTrue(x["source_limits"]["no_original_study_rows_exported"])
        self.assertNotIn("10.1234",str(x))
        self.assertNotIn("Paper A",str(x))

    def test_refuse_unidentified_footer_or_nonbinary_flags(self):
        rows=valid_source()
        rows[-1]["strength sum"]="49"
        with self.assertRaisesRegex(ValueError,"not authenticated"):
            summarize_source(rows,strict=False)
        rows=valid_source()
        rows[0]["interaction strength (IS) reduction"]="unknown"
        with self.assertRaisesRegex(ValueError,"code missing"):
            summarize_source(rows,strict=False)

    def test_refuse_broken_formula_or_moved_aggregate_row(self):
        rows=valid_source()
        rows[1]["effect sum"]="3"
        with self.assertRaisesRegex(ValueError,"not equal"):
            summarize_source(rows,strict=False)
        rows=valid_source()
        rows[-1]["_original_excel_row"]=1
        with self.assertRaisesRegex(ValueError,"not last"):
            summarize_source(rows,strict=False)


if __name__=="__main__":
    unittest.main()
