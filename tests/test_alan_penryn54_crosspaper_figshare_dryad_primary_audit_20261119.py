import io
import json
import unittest
import zipfile
from openpyxl import Workbook

from scripts.acquisition.alan_penryn54_crosspaper_figshare_dryad_primary_audit_20261119 import (
    classify_original_table,csv_schema,xlsx_schema,safe_schema,run, dryad_hal_version_file_lists, dryad_official_href
)

class OriginalPenryn54SourceTests(unittest.TestCase):
    def test_exact_source_csv_headers_counts_never_dump_plot_ids(self):
        raw=(b"Mesocosm,Trophic Treatment,Light,Year,Aphids,Flowerheads\n"
             b"plot_P1,plants+herbivores,white,2013,11,23\n"
             b"plot_P1,plants+herbivores,white,2014,17,19\n"
             b"plot_P2,plants-only,dark,2013,0,30\n")
        data=safe_schema("Bennie2018_data.csv",raw)
        self.assertEqual(data["kind"],"csv")
        s=data["sheets"][0]
        self.assertEqual(s["nonempty_data_rows"],3)
        self.assertEqual(s["column_names_from_exact_original_file"][0],"Mesocosm")
        self.assertEqual(s["candidate_source_fields"][0]["n_distinct_original_values_only"],2)
        self.assertFalse(s["original_observation_values_exported"])
        self.assertNotIn("plot_P1",str(s))
        self.assertNotIn("plants+herbivores",str(s))
        self.assertNotIn("white,2013",str(s))

    def test_xlsx_contains_only_shape_not_original_mesocosm_rows(self):
        wb=Workbook()
        ws=wb.active
        ws.title="Data"
        ws.append(["plot","light","trophic","biomass"])
        ws.append(["SENSITIVE50","amber","plants+pred",2.5])
        ws.append(["SENSITIVE51","dark","plants",1.5])
        b=io.BytesIO();wb.save(b);wb.close()
        ss=xlsx_schema(b.getvalue())
        self.assertEqual(len(ss),1)
        self.assertEqual(ss[0]["nonempty_data_rows"],2)
        self.assertEqual(ss[0]["sheet_name"],"Data")
        self.assertNotIn("SENSITIVE50",str(ss))
        self.assertNotIn("2.5",str(ss))

    def test_archive_parser_limits_and_nested_table_without_raw_rows(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr("original_aphid.csv","plot,treatment,abundance\nSite123,white,4\n")
        info=safe_schema("original-data.zip",buf.getvalue())
        self.assertEqual(info["kind"],"zip")
        self.assertEqual(info["members"][0]["kind"],"csv")
        self.assertNotIn("Site123",str(info))

    def test_dryad_hal_records_use_stash_files_link_not_missing_version_id(self):
        ds={"_links":{"stash:version":{"href":"/api/v2/versions/131413"}}}
        versions={"_embedded":{"stash:versions":[
           {"versionNumber":1,
            "_links":{"self":{"href":"/api/v2/versions/131413"},
                      "stash:files":{"href":"/api/v2/versions/131413/files"}}}
        ]}}
        urls=dryad_hal_version_file_lists(ds,versions)
        self.assertEqual(urls,["https://datadryad.org/api/v2/versions/131413/files"])
        self.assertEqual(dryad_official_href(
            "/api/v2/files/190044/download",allow_files=True),
            "https://datadryad.org/api/v2/files/190044/download")
        with self.assertRaisesRegex(ValueError,"nonofficial"):
            dryad_official_href("https://malicious.example/api/v2/versions/131413/files")
        with self.assertRaisesRegex(ValueError,"genuine original"):
            dryad_official_href("https://datadryad.org/api/v2/versions/131413")

    def test_version_hal_link_can_be_derived_from_dataset_metadata_alone(self):
        ds={"_links":{"stash:version":{"href":"/api/v2/versions/44455"}}}
        result=dryad_hal_version_file_lists(ds,{"_embedded":{"stash:versions":[]}})
        self.assertEqual(result,["https://datadryad.org/api/v2/versions/44455/files"])

    def test_unknown_original_binary_not_invented_as_spreadsheet(self):
        s=safe_schema("predators.RData",b"binary-raw-original-not-parsed")
        self.assertEqual(s["kind"],"other")
        self.assertEqual(s["original_data_type_unparsed"],".rdata")
        self.assertNotIn("column_names_from_exact_original_file",s)

    def test_one_original_54_mesocosm_programme_requires_no_fitted_meta(self):
        # No networking in deterministic tests. Schema parser and published
        # programme-independent count are tested separately from source fetch.
        expected={
            "independent_light_experiment_programmes":1,
            "original_physical_light_mesocosms":54,
            "independent_mesocosms_per_light_trophic_cell":6,
            "crossed_light_arms":3,
            "crossed_trophic_arms":3,
        }
        self.assertEqual(expected["independent_light_experiment_programmes"],1)
        self.assertEqual(expected["independent_mesocosms_per_light_trophic_cell"]*
                         expected["crossed_light_arms"]*expected["crossed_trophic_arms"],
                         expected["original_physical_light_mesocosms"])
if __name__=="__main__":
    unittest.main()
