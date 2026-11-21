import unittest
from scripts.acquisition.alan_deitsch2023_public_edi_source_audit_20261118 import source_table_summary


class EDISourceNonDisclosureAndGrainTests(unittest.TestCase):
    def test_csv_only_column_names_and_plot_counts_not_rows(self):
        raw = ("plot_id,treatment,group,attack\nPLOT_A,lit,spider,1\n"
               "PLOT_A,lit,spider,0\nPLOT_B,dark,beetle,1\n").encode()
        a = source_table_summary(raw)
        self.assertEqual(a["format"], "delimited_text")
        self.assertEqual(a["record_rows"], 3)
        self.assertEqual(a["unique_physical_plot_ID_counts_only"], [2])
        self.assertEqual(a["treatment_levels_counts_only"], [2])
        self.assertNotIn("PLOT_A", str(a))
        self.assertNotIn("PLOT_B", str(a))
        self.assertNotIn("spider", str(a))

    def test_tab_delimited_data(self):
        raw = b"Plot\tTreatment\tcount\nA\tL\t1\nA\tL\t2\nB\tD\t3\n"
        a = source_table_summary(raw)
        self.assertEqual(a["record_rows"], 3)
        self.assertEqual(a["unique_physical_plot_ID_counts_only"], [2])

    def test_not_csv_refuse_to_derive_effect(self):
        a = source_table_summary(b"an unrelated R script or file")
        self.assertEqual(a["format"], "not_structured_csv")


if __name__ == "__main__":
    unittest.main()
