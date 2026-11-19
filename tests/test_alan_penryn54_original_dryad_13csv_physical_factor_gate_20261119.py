import copy
import io
import unittest
import zipfile
import hashlib
from unittest.mock import patch

from scripts.analysis import alan_penryn54_original_dryad_13csv_physical_factor_gate_20261119 as a


class Source54TroughOriginalFactorialTests(unittest.TestCase):
    def sample(self,conflict=False):
        rows=[]
        for l in ("C","A","W"):
            for t in ("P","H","HP"):
                for j in range(6):
                    rows.append({"Plot":f"{l}_{t}_{j}","T treatment":t,"L treatment":l,
                                 "01/06/2014":"1","12/06/2014":"3"})
        slug=[]
        for r in rows:
            if r["T treatment"]!="P":
                x=dict(r)
                x["T.treatment"]=x.pop("T treatment")
                x["L.treatment"]=x.pop("L treatment")
                slug.append(x)
        part=copy.deepcopy(rows)
        if conflict:part[1]["L treatment"]="W"
        return {
          "PeaAphid2014.csv":{"fields":["Plot","T treatment","L treatment","01/06/2014","12/06/2014"],"rows":rows},
          "lotus_p_2014.csv":{"fields":["Plot","T treatment","L treatment","01/06/2014","12/06/2014"],"rows":part},
          "slugs2014.csv":{"fields":["Plot","T.treatment","L.treatment","01/06/2014","12/06/2014"],"rows":slug}
        }

    def test_same_54_randomised_source_trough_ids_across_all_endpoints(self):
        x=a.table_plot_factor_grain(self.sample(),strict=True)
        self.assertEqual(x["status"],"ORIGINAL_54_SOURCE_FACTORIAL_PLOT_JOIN_VERIFIED")
        self.assertEqual(x["original_unique_physical_IDs_across_13_source_CSVs"],54)
        self.assertEqual(x["source_same_plot_different_treatment_code_conflict_count"],0)
        self.assertEqual(len(x["all_original_physical_light_trophic_cell_counts_only"]),9)
        self.assertEqual(set(u["unique_original_troughs"] for u in x["all_original_physical_light_trophic_cell_counts_only"]),{6})
        self.assertEqual(x["original_source_factor_codes"]["trophic"],["H","HP","P"])
        self.assertEqual(x["original_source_factor_codes"]["light"],["A","C","W"])
        self.assertEqual(x["original_13_source_CSV_physical_grain"]["slugs2014.csv"]["original_physical_unit_count"],36)
        self.assertTrue(x["post_reintroduction_slugs_2015_excluded_from_authors_original_fit"])
        self.assertNotIn("C_P_0",str(x))

    def test_crosspaper_treatment_reallocation_does_not_fake_source_join(self):
        x=a.table_plot_factor_grain(self.sample(conflict=True),strict=True)
        self.assertEqual(x["status"],"FACTORIAL_SOURCE_STOP_OR_SOURCE_CODE_INCONSISTENCY")
        self.assertGreater(x["source_same_plot_different_treatment_code_conflict_count"],0)
        self.assertFalse(x["cross_file_exact_original_plot_factor_allocations_consistent"])

    def test_legitimate_original_manifest_and_dataset_id_frozen(self):
        self.assertEqual(a.SOURCE_VERSION_ID,15138)
        self.assertEqual(a.ARCHIVE_SHA256,
            "5dc5bc93e48758364d068c113edb84e05b81f99b22827c4c6ee0c14f42ba3130")
        self.assertEqual(len(a.MANIFEST),13)
        self.assertEqual(a.SOURCE_ROWS["slugs2014.csv"],36)
        self.assertEqual(a.SOURCE_ROWS["CNratio.csv"],51)
        self.assertEqual(a.SOURCE_ROWS["lotus_p_2013.csv"],54)
        self.assertEqual(a.SOURCE_SLUG_FIT_YEARS,(2013,2014))

    def test_reject_changed_source_archive_digest(self):
        fake=io.BytesIO()
        with zipfile.ZipFile(fake,"w") as z:
            z.writestr("PeaAphid2014.csv","Plot,T treatment,L treatment\nC1,H,W\n")
        with self.assertRaisesRegex(ValueError,"archive SHA256 drift"):
            a.parse_zip_archive(fake.getvalue(),strict=True)

    def test_unnamed_first_column_source_is_id_not_empty(self):
        tab={"point_quadrat_110713.csv":{
          "fields":["","T treatment","L treatment","Lotus pedunculatus"],
          "rows":[{"":"10","T treatment":"H","L treatment":"C","Lotus pedunculatus":"2"}]
        }}
        x=a.table_plot_factor_grain(tab,strict=False)
        self.assertEqual(x["original_unique_physical_IDs_across_13_source_CSVs"],1)
        self.assertEqual(x["original_13_source_CSV_physical_grain"]["point_quadrat_110713.csv"]["source_id_column_name"],
                         "<ORIGINAL_UNNAMED_FIRST_COLUMN>")
        self.assertFalse(x["cross_file_exact_original_plot_factor_allocations_consistent"])

if __name__=="__main__":
    unittest.main()
