"""Test no fake SVE pollen effects, source-study clustering and rare pollination roles."""
import unittest
from scripts.analysis.alan_page2021_original_macro_functional_network_mismatch import (
 candidate_groups,boot_study_rate,EFF,UNITS
)
def record(study,plant,taxon,guild,q,visits,metric=EFF,units=UNITS):
    return {"studyID":study,"plant.species":plant,"study.year":"2015",
        "effectiveness.metric":metric,"effectiveness.average.units":units,
        "visit.metric":"number of visits per unit time per unit area",
        "pollinator.taxon":taxon,"pollinator.group":guild,
        "effectiveness.average":str(q),"visit.value":str(visits)}
class PageMacroFunctionalSourceTests(unittest.TestCase):
    def test_most_visiting_not_most_pollen_and_rare_important(self):
        z=[record("a","P","Apis","bee",.1,100),
           record("a","P","Diptera","fly",10,4)]
        group,excl,cells=candidate_groups(z)
        self.assertEqual(len(group),1)
        self.assertFalse(group[0]["top_visitor_equals_top_functional_pollen"])
        self.assertEqual(len(group[0]["rare_under5pct_visitation_over10pct_functional_pollen"]),1)
        self.assertEqual(group[0]["top_visitor_taxon"],"Apis")
        self.assertEqual(group[0]["top_functional_taxon"],"Diptera")

    def test_different_Pollen_SVE_unit_and_study_not_pooled(self):
        z=[record("a","P","A","bee",2,10),
           record("b","P","B","fly",10,10),
           record("a","P","C","fly",10,10,metric="fruit set")]
        group,excluded,cells=candidate_groups(z)
        self.assertEqual(len(group),0)

    def test_original_Vanessa_source_only_taxonomy_exclusion(self):
        z=[record("a","P","Vanessa carye","moth",0,1),
           record("a","P","Sphingidae","moth",1,3),
           record("a","P","Diptera","fly",1,1)]
        allx,_,_=candidate_groups(z,False)
        fixed,_,_=candidate_groups(z,True)
        self.assertEqual(len(allx[0]["observed_source_pollinator_guilds"]),2)
        self.assertEqual(len(fixed[0]["observed_source_pollinator_guilds"]),2)
        self.assertEqual(allx[0]["original_source_pollinator_taxa"],3)
        self.assertEqual(fixed[0]["original_source_pollinator_taxa"],2)

    def test_boot_study_when_little_source_support_fails(self):
        x=[{"source_study":"a"}]*3
        self.assertEqual(boot_study_rate(x)["status"],"insufficient_source_independent_studies_STOP")
if __name__=="__main__":unittest.main()
