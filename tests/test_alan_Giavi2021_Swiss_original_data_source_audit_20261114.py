#!/usr/bin/env python3
"""Swiss 2016 original light experiment stable source and schema guards."""
import hashlib,sys,unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_Giavi2021_Swiss_original_data_source_audit_20261114 import SRC, dataframe_info

class OriginalSwissData(unittest.TestCase):
    def test_zenodo_published_original_md5s_fixed(self):
        self.assertEqual(SRC["data.Rdata"],"f60b21d82c163cbb53520fda0bc61985")
        self.assertEqual(SRC["readme"],"3fccd78016a9cd41f122907e1561fa19")
        self.assertEqual(SRC["script.R"],"1ba8c516766a879b326d0b14bba69f1e")
    def test_no_biology_from_missing_plot_records(self):
        z=pd.DataFrame({"plant":["A","B",None],"treatment":["dark","lit","lit"],
                        "site":[1,2,3],"visits":[5,6,None]})
        info=dataframe_info(z)
        self.assertEqual(info["rows"],3)
        self.assertEqual(info["candidate_design_fields"]["treatment"]["distinct_count"],2)
        self.assertGreater(info["source_missing_fraction_per_col"]["plant"],0)
        self.assertFalse(info["individual_source_observation_rows_exported"])
        self.assertNotIn("source_rows",info)

if __name__=="__main__":
    unittest.main()
