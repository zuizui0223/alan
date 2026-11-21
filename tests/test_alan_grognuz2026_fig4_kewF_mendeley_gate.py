"""Ensure genuine UK F scores and Mendeley source gate cannot invent wet taxa."""
import json,unittest
from unittest.mock import patch
from scripts.acquisition import alan_grognuz2026_fig4_kewF_mendeley_source_gate as m

class Fig4OriginalSourceGates(unittest.TestCase):
    def test_fig4_source_label_counts_and_unresolved_names(self):
        self.assertEqual(len(m.PLANTS),16)
        self.assertEqual(len(set(m.PLANTS)),16)
        self.assertIn("Sonchus sp.",m.PLANTS)
        self.assertIn("Epilobium sp.",m.PLANTS)
        self.assertIn("Silene pratensis",m.PLANTS)
    def test_original_kew_csv_no_genus_imputation(self):
        h=("taxon_name_binom,F_PLANTATT,F_Doring,taxonomic_status,aggregate\n")
        rows=[f"Other{i},NA,NA,Accepted,\n" for i in range(3000)]
        rows.extend([
            "Achillea millefolium,5,5,Accepted,\n",
            "Centaurea cyanus,5,NA,Accepted,\n",
            "Trifolium incarnatum,2,2,Accepted,\n",
            "Centaurea jacea,NA,NA,Accepted,\n"
        ])
        raw=(h+"".join(rows)).encode("latin1")
        with patch.object(m,"K_GIT_BLOB",m.git_blob(raw)):
            s=m.kew_traits(raw)
        self.assertEqual(s["source_published_fig4_labels"],16)
        self.assertEqual(s["F_PLANTATT_scored"],3)
        self.assertEqual(s["F_PLANTATT_scored_wet_gt7"],0)
        a={z["paper_Figure4_label"]:z for z in
           s["figure4_source_name_and_exact_moisture"]}
        self.assertIsNone(a["Sonchus sp."]["Kew_F_PLANTATT"])
        self.assertTrue(a["Sonchus sp."]["not_a_species_level_label"])
        self.assertIsNone(a["Centaurea jacea"]["Kew_F_PLANTATT"])
        self.assertEqual(a["Trifolium incarnatum"]["Kew_F_PLANTATT"],2)
    def test_mendeley_unauthorized_never_becomes_original_file(self):
        def no_access(label,url):
            return {"label":label,"url":url,"status":"SOURCE_ACCESS_STOP",
                    "error_type":"HTTPError","http_code":401}
        with patch.object(m,"probe",no_access):
            source=m.mendeley_audit()
        self.assertEqual(source["status"],
                         "ORIGINAL_MENDELEY_PUBLISHED_FILE_SOURCE_STOP")
        self.assertFalse(source["source_json_metadata_readable"])
        self.assertFalse(source["source_raw_experimental_phenology_rows_downloaded"])
        self.assertFalse(source["plant_by_site_original_light_response_fitted"])
if __name__=="__main__":unittest.main()
