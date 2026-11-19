import copy
import unittest

from scripts.analysis.alan_anic2022_Lotus_original_54trough_biomass_welch_20261119 import source_contrasts
from scripts.analysis.alan_anic2022_original_pdf_biomass_54trough_parse_20261119 import FILE_SHA256

# These are genuine author original 2022 PDF raw-trough-cell source summaries
# independently SHA-pinned from source reconstruction actions #44526503325.
# Test fixture ONLY stores source means, SD, n, not individual plant records.
LOTUS={
 "L1":[(59.41833333333333,42.358),(82.91666666666667,25.988),(76.33000000000000,48.799)],
 "L2":[(48.08666666666667,27.353),(52.94000000000000,30.251),(68.16166666666666,30.396)],
 "L3":[(34.48333333333333,21.660),(58.98166666666666,24.702),(47.35333333333333,29.523)]
}


def fixture():
    rows=[]
    for light,values in LOTUS.items():
        for ti,(mean,sd) in enumerate(values,1):
            rows.append({
                "species":"Lotus pedunculatus",
                "original_light_code":light,
                "original_trophic_code":f"T{ti}",
                "original_mesocosms_with_this_species":6,
                "original_biomass_total_g_arithmetic_mean":mean,
                "original_biomass_total_g_sample_SD":sd,
            })
    return {
        "source_original_figshare_sha256":FILE_SHA256,
        "experimental_light_trophic_full_factorial_verified":True,
        "original_unique_physical_mesocosms_in_biomass_table":54,
        "original_species_rows_parsed":665,
        "source_original_species_by_light_trophic_cell_descriptive_means":rows,
    }


class OriginalAnic2022IndependentBiomassTests(unittest.TestCase):
    def test_lotus_original_independent_trough_stratified_means_and_welch(self):
        result=source_contrasts(fixture())
        contrasts=result["source_original_exploratory_light_pair_contrasts"]
        self.assertEqual(len(contrasts),3)
        obs=next(c for c in contrasts if c["original_unmapped_light_code_contrast"]=="L3 minus L1")
        self.assertEqual(obs["original_physical_mesocosms_each_light_arm"],18)
        self.assertEqual(obs["source_original_trophic_strata"],3)
        self.assertAlmostEqual(obs["difference_g"],-25.9488888888889,delta=.01)
        self.assertAlmostEqual(obs["independent_trough_Welch_approx_SE_g"],11.223,delta=.05)
        self.assertLess(obs["exploratory_unadjusted_t_P"],.04)
        self.assertGreater(obs["exploratory_Holm_multiplicity_adjusted_P_three_light_pairs"],.05)
        self.assertLess(obs["exploratory_unadjusted_95pct_CI_g"][0],0)
        self.assertLess(obs["exploratory_unadjusted_95pct_CI_g"][1],0)
        self.assertFalse(result["source_original_light_code_mapping_verified_to_unlit_amber_white"])
        self.assertFalse(result["source_2018_Dryad_herbivore_rows_joined_by_original_trough"])
        self.assertTrue(result["no_original_plot_IDs_or_species_by_trough_values_exported"])

    def test_real_original_balanced_3x3_ANOVA_does_not_overstate_light_p(self):
        out=source_contrasts(fixture())
        anova=out["original_54trough_Lotus_balanced_3x3_factorial_ANOVA"]
        self.assertEqual(anova["independently_assigned_physical_trough_n"],54)
        self.assertEqual(anova["error_df"],45)
        self.assertEqual(anova["physical_troughs_per_cell"],6)
        f_light=anova["factors"]["LIGHT"]
        self.assertAlmostEqual(f_light["F"],2.9757,delta=.02)
        self.assertAlmostEqual(f_light["unadjusted_parametric_p_homoskedastic"],.06113,delta=.002)
        self.assertEqual(anova["factors"]["LIGHT_X_TROPHIC"]["df"],4)
        self.assertGreater(anova["factors"]["LIGHT_X_TROPHIC"]["unadjusted_parametric_p_homoskedastic"],.80)
        self.assertTrue(anova["not_author_original_AIC_selected_GLM"])
        self.assertTrue(anova["homoskedastic_assumption_unverified"])
        self.assertTrue(anova["not_corrected_across_13_species"])

    def test_changed_physical_plot_counts_or_source_SHA_rejected(self):
        f=fixture()
        f["original_source_figshare_sha256"]="fabricated-source"
        with self.assertRaisesRegex(ValueError,"SHA drift"):
            source_contrasts(f)
        f=fixture()
        f["source_original_species_by_light_trophic_cell_descriptive_means"][0]["original_mesocosms_with_this_species"]=500
        with self.assertRaisesRegex(ValueError,"n must be six"):
            source_contrasts(f)
        f=fixture()
        f["experimental_light_trophic_full_factorial_verified"]=False
        with self.assertRaisesRegex(ValueError,"54-trough factorial"):
            source_contrasts(f)

    def test_missing_source_trophic_cell_cannot_fake_lotus_representativity(self):
        f=fixture()
        f["source_original_species_by_light_trophic_cell_descriptive_means"].pop()
        with self.assertRaisesRegex(ValueError,"every randomization cell"):
            source_contrasts(f)


if __name__=="__main__":
    unittest.main()
