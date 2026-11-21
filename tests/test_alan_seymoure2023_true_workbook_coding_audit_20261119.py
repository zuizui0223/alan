import hashlib
import io
import unittest
from unittest.mock import patch

from openpyxl import Workbook
from scripts.analysis import alan_seymoure2023_true_workbook_coding_audit_20261119 as a


class OriginalSeymoureCodingAuditTests(unittest.TestCase):
    def source(self):
        wb=Workbook()
        s=wb.active
        s.title="database"
        headers=[
            "interaction_type","method","Consumer_taxon","Taxa1 (consumer) detailed",
            "Resource_taxon","Taxa2 (resource) detailed",
            "interaction strength (IS) reduction","IS enhancement","strength sum",
            "temporal mismatch induced","spatial mismatch induced","effect sum",
            "comments","Authors","Article Title","Source Title",
            "Publication Year","DOI","DOI Link","WoS Categories"]
        s.append(headers)
        s.append(["predation","experiment","spiders","M", "insects","Diptera",
                  1,0,1,0,0,1,"","A","paper A","Ecology",2020,
                  "10.1234/A","https://doi.org/10.1234/A","ecology"])
        s.append(["pollination","experiment","insects","bees","plants","Asteraceae",
                  0,1,1,0,1,2,"","A","paper A","Ecology",2020,
                  "https://doi.org/10.1234/a","","ecology"])
        s.append(["herbivory","observational","insects","moths","plants","seed",
                  0,0,0,1,0,1,"","B","paper B","Botany",2021,
                  "","","botany"])
        stream=io.BytesIO()
        wb.save(stream)
        wb.close()
        return stream.getvalue()

    def test_study_identity_is_not_row_identity(self):
        b=self.source()
        with patch.object(a,"SOURCE_SHA256",hashlib.sha256(b).hexdigest()):
            result=a.source_workbook_audit(b)
        s=result["source_observed"]
        self.assertEqual(s["original_excel_data_nonempty_rows"],3)
        self.assertEqual(s["source_distinct_publication_IDs_DOI_then_title"],2)
        self.assertEqual(s["source_distinct_DOIs_nonmissing"],1)
        self.assertEqual(s["source_DOI_missing_rows"],1)
        self.assertEqual(s["source_publications_by_number_of_rows"],{1:1,2:1})
        self.assertEqual(s["type_original_value_counts"],
                         {"herbivory":1,"pollination":1,"predation":1})
        self.assertFalse(result["gate"]["meta_effect_estimated"])
        self.assertFalse(result["gate"]["source_coding_direction_interpreted"])
        self.assertTrue(result["gate"]["source_excel_row_level_records_not_exported"])
        self.assertNotIn("10.1234",str(result))
        self.assertNotIn("paper A",str(result))

    def test_wrong_original_source_hash_must_fail(self):
        with self.assertRaisesRegex(ValueError,"Source SHA drift"):
            a.source_workbook_audit(self.source())

    def test_no_undocumented_codebook_relabeling(self):
        b=self.source()
        with patch.object(a,"SOURCE_SHA256",hashlib.sha256(b).hexdigest()):
            result=a.source_workbook_audit(b)
        values=result["literal_source_codebook"]["strength_reduction"]["observed_original_value_counts"]
        self.assertEqual(values,[{"original_value":"0","n":2},{"original_value":"1","n":1}])


if __name__=="__main__":
    unittest.main()
