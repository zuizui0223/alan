import io
import unittest
from openpyxl import Workbook

from scripts.acquisition.alan_perkin2014_four_streams_AppendixS2_source_20261120 import (
    audit, workbook_source_shape,DOI,WILEY_PUBLISHER_CANDIDATES
)


class Perkin2014FourStreamsOriginalAppendixTests(unittest.TestCase):
    def fixture(self):
        w=Workbook()
        sheet=w.active
        sheet.title="Original reach data"
        sheet.append(["Stream","Reach","Treatment","Drift","Emergence"])
        sheet.append(["PrivateRiver1","Reach1","Light",25,11])
        sheet.append(["PrivateRiver1","Reach2","Dark",50,9])
        out=io.BytesIO()
        w.save(out)
        w.close()
        return out.getvalue()

    def test_only_original_workbook_shapes_without_exposing_trout_or_drift(self):
        result=workbook_source_shape(self.fixture())
        self.assertEqual(len(result),1)
        self.assertEqual(result[0]["original_nonempty_rows_including_headers"],3)
        self.assertIn("Stream",result[0]["early_header_shapes_only"][0]["candidate_header_labels"])
        self.assertTrue(result[0]["original_stream_or_reach_ID_values_not_exported"])
        self.assertNotIn("PrivateRiver1",str(result))
        self.assertNotIn("Reach2",str(result))
        self.assertNotIn("50",str(result))

    def test_publisher_HTTP_denial_is_source_stop_not_zero_effect(self):
        def reject(url):
            raise ValueError("Blocked official supplement")
        x=audit(fetch=reject)
        self.assertEqual(len(x["official_publisher_source_attempts"]),len(WILEY_PUBLISHER_CANDIDATES))
        self.assertEqual(x["original_independent_stream_pair_count"],4)
        self.assertEqual(x["source_status"],"SOURCE_STOP")
        self.assertFalse(x["night_drift_original_field_meta_effect_fitted"])
        self.assertIsNone(x["source_original_workbook_sha256"])

    def test_successful_authenticated_original_bytes_are_only_source_shape(self):
        x=audit(fetch=lambda url:self.fixture())
        self.assertEqual(x["source_status"],"ORIGINAL_WILEY_APPX_S2_BYTES_AND_SHEETS_VERIFIED")
        self.assertEqual(x["original_independent_stream_pair_count"],4)
        self.assertIsNotNone(x["source_original_workbook_sha256"])
        self.assertFalse(x["aquatic_emergence_original_field_meta_effect_fitted"])
        self.assertNotIn("PrivateRiver1",str(x))

if __name__=="__main__":
    unittest.main()
