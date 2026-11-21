import unittest
from scripts.analysis.alan_anic2022_original_pdf_biomass_54trough_parse_20261119 import (
    parse_table_line,identify_species_prefix,BIOMASS_SPECIES,OUTPUT_SCHEMA
)


class AnicOriginal54TroughBiomassPDFTests(unittest.TestCase):
    def test_exact_species_prefix_and_6_original_biomass_columns(self):
        q=parse_table_line(
            "Achillea millefolium PLOT1 L1 T1 3.50 2.25 1.25")
        self.assertEqual(q["species"],"Achillea millefolium")
        self.assertEqual(q["meso"],"PLOT1")
        self.assertEqual(q["light"],"L1")
        self.assertEqual(q["trophic"],"T1")
        self.assertEqual(q["total"],3.50)
        self.assertAlmostEqual(q["total"],q["root"]+q["shoot"])
        self.assertEqual(len(BIOMASS_SPECIES),18)

    def test_non_target_header_and_nonbiomass_table_are_excluded(self):
        self.assertIsNone(parse_table_line("Species Mesocosm L T Total (g) root (g) shoot (g)"))
        self.assertIsNone(parse_table_line("L = Light treatment T = Trophic treatment"))
        self.assertIsNone(parse_table_line("This is a figure data archive"))
        self.assertEqual(OUTPUT_SCHEMA,
            "alan_anic2022_original_54trough_PDF_biomass_parsed_counts_and_cell_means_v1")

    def test_invalid_biomass_and_missing_columns_refuse_fabricated_data(self):
        with self.assertRaisesRegex(ValueError,"unexpected"):
            parse_table_line("Lotus pedunculatus PLOT1 L2 T3 1.0 2.0")
        with self.assertRaisesRegex(ValueError,"nonnumeric"):
            parse_table_line("Lotus pedunculatus PLOT1 L2 T3 NA 1.0 2.0")
        with self.assertRaisesRegex(ValueError,"impossible"):
            parse_table_line("Lotus pedunculatus PLOT1 L2 T3 -3.1 0.1 0.2")
        self.assertIsNone(parse_table_line("Unrelated taxon PLOT1 L2 T1 0.1 0.2 0.3"))

    def test_source_species_handling_preserves_biological_names(self):
        x=identify_species_prefix("Lotus corniculatus 41 LED P 2.10 0.4 1.7")
        self.assertEqual(x[0],"Lotus corniculatus")
        self.assertEqual(x[1],["41","LED","P","2.10","0.4","1.7"])
        y=identify_species_prefix("Lotus pedunculatus 12 AMB HERB 2.0 1.0 1.0")
        self.assertEqual(y[0],"Lotus pedunculatus")


if __name__=="__main__":
    unittest.main()
