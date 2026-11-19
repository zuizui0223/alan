import copy
import json
import unittest
from pathlib import Path

from scripts.analysis.alan_grognuz2026_figure4_phenology_reorganization_20261119 import analyze

SOURCE=Path("results/all_species_forest/biotic_interactions/grognuz2026_figure4_16plant_phenology_visual_approx_source_20261119.json")


class PublishedFigurePhenologyOpportunityTests(unittest.TestCase):
    def setUp(self):
        self.source=json.loads(SOURCE.read_text())

    def test_original_one_programme_14_sitepair_years_not_64_tests(self):
        z=analyze(self.source)
        self.assertEqual(z["evidence"]["original_programme_count"],1)
        self.assertEqual(z["evidence"]["independent_site_pair_years"],14)
        self.assertEqual(z["evidence"]["nested_original_inner_phenology_plots"],108)
        self.assertEqual(z["evidence"]["source_flowering_species"],16)
        self.assertEqual(z["evidence"]["author_reported_week_scale_shift_species"],12)
        self.assertTrue(z["gates"]["16_species_same_14_site_pair_year_source_not_16_independent_experiments"])
        self.assertTrue(z["gates"]["no_original_plot_level_SE_or_species_covariance"])

    def test_within_species_window_directional_heterogeneity(self):
        z=analyze(self.source)
        x=z["figure4_approx_bloom_window_changes"]
        self.assertEqual(x["observed_min_max_days"],[-24,24])
        self.assertAlmostEqual(x["species_mean_of_LAST_minus_FIRST_days"],1)
        self.assertAlmostEqual(x["species_median_of_LAST_minus_FIRST_days"],1.5)
        self.assertEqual(len(x["guaranteed_more_than_10d_extension_only_IF_visual_approx_within_declared_error_days"]),3)
        self.assertEqual(len(x["guaranteed_more_than_10d_contraction_only_IF_visual_approx_within_declared_error_days"]),2)
        self.assertTrue(x["not_a_fitted_duration_outcome"])
        self.assertTrue(x["not_actual_coflowering_length_or_fitness"])

    def test_inter_species_shift_differences_are_not_coflowering(self):
        z=analyze(self.source)
        p=z["figure4_approx_interspecies_PEAK_response_difference"]
        self.assertEqual(p["distinct_species_pair_differences"],120)
        self.assertEqual(p["threshold_absolute_differential_shift_days"],14)
        self.assertEqual(p["approx_pairs_differential_shift_at_least_threshold"],23)
        self.assertEqual(p["minimum_guaranteed_only_IF_fig4_readout_within_illustrative_perturbation"],7)
        self.assertEqual(p["maximum_possible_only_IF_fig4_readout_within_illustrative_perturbation"],60)
        self.assertTrue(p["NOT_true_overlap_or_relative_flowering_date"])
        self.assertTrue(p["120_pairs_not_independent_experiments"])

    def test_author_petals_diel_effect_is_not_seed_fitness(self):
        z=analyze(self.source)
        a=z["author_original_diel_petals"]
        self.assertEqual(a["Malva_moschata_day_flowering_morning_opening_min"],-147)
        self.assertEqual(a["Silene_noctiflora_morning_start_closing_min"],80)
        self.assertTrue(a["only_Malva_CI_excludes_zero"])
        self.assertTrue(a["diel_opening_not_pollen_deposition_or_seed_set"])

    def test_false_exact_original_raw_or_independent_site_design_rejected(self):
        s=copy.deepcopy(self.source)
        s["design"]["matched_light_dark_site_pairs_in_2022"]=50
        with self.assertRaisesRegex(ValueError,"pair-year"):
            analyze(s)
        s=copy.deepcopy(self.source)
        s["source_figure_4"]["not_original_data_or_published_exact_coefficients"]=False
        with self.assertRaisesRegex(ValueError,"original exact"):
            analyze(s)
        s=copy.deepcopy(self.source)
        s["design"]["actual_pollinator_visitation_or_seed_fitness_in_same_experiment_measured"]=True
        with self.assertRaisesRegex(ValueError,"DOES NOT"):
            analyze(s)
        s=copy.deepcopy(self.source)
        s["source_figure_4"]["source_species"].append(s["source_figure_4"]["source_species"][0])
        with self.assertRaisesRegex(ValueError,"16 original"):
            analyze(s)


if __name__=="__main__":
    unittest.main()
