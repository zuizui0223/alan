"""Strict pairing of Page et al 2021 original SVE by true study and plant."""
import unittest
from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import (
    compare, within_original_sameplant_pollen_flux, taxon_verified_source_sensitivity, QKEY, VKEY, PKEY
)
def row(study,plant,year,guild,metric="stigma pollen",
        unit="grains",visitmetric="visits/hour",q="10"):
    return {"studyID":study,"plant.species":plant,"study.year":year,
            "pollinator.group":guild,"pollinator.taxon":guild+"_sp",
            "effectiveness.metric":metric,"effectiveness.average.units":unit,
            "visit.metric":visitmetric,"effectiveness.average":q,"visit.value":"4"}
class PageSourceSamePlantPairing(unittest.TestCase):
    def test_different_studies_and_plants_never_pair(self):
        a=[row("a","Plant1","2020","moth"),row("b","Plant1","2020","fly"),
           row("a","Plant2","2020","fly")]
        self.assertEqual(compare(a,QKEY,True)["original_moth_fly_matched_cells"],0)
    def test_different_pollen_vs_seed_outcome_never_pair(self):
        a=[row("a","P","2020","moth",metric="pollen"),
           row("a","P","2020","fly",metric="seeds")]
        self.assertEqual(compare(a,PKEY,False)["original_moth_fly_matched_cells"],1)
        self.assertEqual(compare(a,QKEY,True)["original_moth_fly_matched_cells"],0)
    def test_same_study_plant_year_metric_actual_ratio_only(self):
        a=[row("a","P","2020","moth",q="10"),
           row("a","P","2020","fly",q="5")]
        r=compare(a,VKEY,True)
        self.assertEqual(r["original_moth_fly_matched_cells"],1)
        self.assertAlmostEqual(r["source_comparisons"][0]["source_fly_moth_mean_SVE_ratio"],.5)
    def test_visit_denominators_are_separate_gate(self):
        a=[row("a","P","2020","moth",visitmetric="visits/min"),
           row("a","P","2020","fly",visitmetric="visits/hour")]
        self.assertEqual(compare(a,QKEY,True)["original_moth_fly_matched_cells"],1)
        self.assertEqual(compare(a,VKEY,True)["original_moth_fly_matched_cells"],0)
    def test_functional_stigma_pollen_flux_not_visit_rate(self):
        a=[row("a","Plant","2020","moth",metric="number of conspecific pollen grains",
               unit="mean",q="10"),
           row("a","Plant","2020","fly",metric="number of conspecific pollen grains",
               unit="mean",q="5")]
        a[0]["visit.value"]="2";a[1]["visit.value"]="4"
        z=within_original_sameplant_pollen_flux(a)
        d=z["source_paired_sameplant_records"][0]
        self.assertEqual(d["source_same_metric_pollen_times_visits_eligible"],True)
        self.assertAlmostEqual(d["visitation_ratio_fly_vs_moth"],2.)
        self.assertAlmostEqual(d["estimated_conspecific_pollen_flux_ratio_fly_vs_moth"],1.)
        self.assertTrue(d["source_functional_ratio_is_sve_times_visit_not_observed_fitness"])

    def test_no_infinite_ratio_when_moth_q_zero(self):
        a=[row("a","Plant","2020","moth",metric="number of conspecific pollen grains",
               unit="mean",q="0"),
           row("a","Plant","2020","fly",metric="number of conspecific pollen grains",
               unit="mean",q="7")]
        z=within_original_sameplant_pollen_flux(a)
        self.assertIsNone(z["source_paired_sameplant_records"][0][
            "estimated_conspecific_pollen_flux_ratio_fly_vs_moth"])

    def test_other_visitation_measure_cannot_make_flux(self):
        a=[row("a","Plant","2020","moth",metric="number of conspecific pollen grains",
               unit="mean",visitmetric="other"),
           row("a","Plant","2020","fly",metric="number of conspecific pollen grains",
               unit="mean",visitmetric="other")]
        z=within_original_sameplant_pollen_flux(a)
        self.assertFalse(z["source_paired_sameplant_records"][0][
            "source_same_metric_pollen_times_visits_eligible"])

    def test_original_Nymphalidae_Vanessa_not_moth_source_preserved(self):
        a=[row("study","Plant","2020","moth",metric="number of conspecific pollen grains",
               unit="mean",q="1.3"),
           row("study","Plant","2020","moth",metric="number of conspecific pollen grains",
               unit="mean",q="0"),
           row("study","Plant","2020","fly",metric="number of conspecific pollen grains",
               unit="mean",q="4")]
        a[0]["pollinator.taxon"]="Hyles annei"
        a[1]["pollinator.taxon"]="Vanessa carye"
        a[1]["visit.value"]="0.01"
        source_ratio=compare(a,VKEY,True)["source_comparisons"][0]["source_fly_moth_mean_SVE_ratio"]
        checked=taxon_verified_source_sensitivity(a)
        after=checked["corrected_taxon_reference_sensitivity_SVE"]["source_comparisons"][0]
        self.assertAlmostEqual(source_ratio,4/.65)
        self.assertAlmostEqual(after["source_fly_moth_mean_SVE_ratio"],4/1.3)
        self.assertTrue(checked["original_source_taxonomic_mismatch_not_silently_relabelled"])

if __name__=="__main__":unittest.main()
