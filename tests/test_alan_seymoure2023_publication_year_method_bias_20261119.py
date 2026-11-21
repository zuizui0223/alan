import unittest
from scripts.analysis.alan_seymoure2023_publication_year_method_bias_20261119 import profile


def row(i, typ, flags, study, year, method=""):
    d,e,t,s=flags
    return {
        "_original_excel_row":i,"interaction_type":typ,
        "method":method,"Consumer_taxon":"animal","Resource_taxon":"resource",
        "interaction strength (IS) reduction":str(d),
        "IS enhancement":str(e),"strength sum":str(d+e),
        "temporal mismatch induced":str(t),
        "spatial mismatch induced":str(s),
        "effect sum":str(d+e+t+s),
        "Authors":study,"Article Title":"Study "+study,
        "Publication Year":str(year),"DOI":"10.54321/"+study,"DOI Link":""
    }


def source():
    records=[
        row(2,"predator-prey",(0,1,0,1),"A",2018,method="field"),
        row(3,"pollination",(1,0,0,0),"A",2018),
        row(4,"herbivory",(1,0,0,1),"B",2021),
    ]
    f=row(5,"",(0,0,0,0),"",2021)
    f.update({
        "interaction_type":"",
        "Authors":"","Article Title":"","Publication Year":"",
        "DOI":"",
        "interaction strength (IS) reduction":"29",
        "IS enhancement":"28","strength sum":"57",
        "temporal mismatch induced":"7",
        "spatial mismatch induced":"27",
        "effect sum":""
    })
    return records+[f]


class PublicationMethodYearCoverageTests(unittest.TestCase):
    def test_same_study_multiple_interaction_year_not_double_counted(self):
        result=profile(source(),strict=False)
        self.assertEqual(result["source_publication_count"],2)
        self.assertEqual(result["original_method_observed_studies"],1)
        early=result["year_stratified_source_coding"]["pre2020"]
        late=result["year_stratified_source_coding"]["2020plus"]
        self.assertEqual(early["publication_counts"]["n_studies"],1)
        self.assertEqual(early["publication_counts"]["coded_any_spatial"],1)
        self.assertEqual(late["publication_counts"]["n_studies"],1)
        self.assertEqual(late["interaction_type_counts"]["herbivory"]["coded_decrease_only"],1)
        self.assertTrue(result["source_caveats"]["method_blank_is_not_natural_light_control"])
        self.assertNotIn("10.54321",str(result))

    def test_year_drift_raises(self):
        cases=source()
        cases[1]["Publication Year"]="2021"
        with self.assertRaisesRegex(ValueError,"Conflicting paper"):
            profile(cases,strict=False)
        cases=source()
        cases[0]["Publication Year"]="2030"
        cases[1]["Publication Year"]="2030"
        with self.assertRaisesRegex(ValueError,"plausible interval"):
            profile(cases,strict=False)


if __name__=="__main__":
    unittest.main()
