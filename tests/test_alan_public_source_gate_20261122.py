"""Offline source-gate parser and access-safety tests; no internet dependency."""
import importlib.util
import unittest
from pathlib import Path

PATH=Path(__file__).resolve().parents[1]/'scripts'/'acquisition'/'alan_night_insect_public_source_gate_20261122.py'
spec=importlib.util.spec_from_file_location('alan_source_gate',PATH)
gate=importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)

class SafeSourceGateTests(unittest.TestCase):
    def test_hostname_allowlist(self):
        for u in ['https://dataverse.nl/api/datasets/1','https://api.figshare.com/v2/articles/29436527','https://ndownloader.figshare.com/files/1']:
            gate.check_url(u)
        for u in ['http://dataverse.nl/','https://evil.example/a','https://dataverse.nl.evil.test/a']:
            with self.assertRaises(ValueError):gate.check_url(u)

    def test_dataverse_manifest_file_id_and_license(self):
        p={'data':{'latestVersion':{'license':{'name':'CC0'},'files':[{'label':'count.csv','dataFile':{'id':19,'filesize':39,'restricted':False,'contentType':'text/csv'}}]}}}
        f, license=gate.dataverse_files(p)
        self.assertEqual(len(f),1)
        self.assertEqual(f[0]['id'],19)
        self.assertTrue(f[0]['download_url'].startswith('https://dataverse.nl/api/access/datafile/19'))
        self.assertEqual(license['name'],'CC0')

    def test_figshare_manifest(self):
        p={'files':[{'id':42,'name':'Roadlights_long.csv','size':35,'download_url':'https://ndownloader.figshare.com/files/42'}],'license':{'name':'CC BY'}}
        f,license=gate.figshare_files(p)
        self.assertEqual(f[0]['label'],'Roadlights_long.csv')
        self.assertEqual(license['name'],'CC BY')

    def test_csv_schema_not_effect_or_raw_redistribution(self):
        s=gate.tabular_summary('capture.csv',b'trap_id,date,light,caught\nT1,2023-07-11,dark,0\nT2,2023-07-11,red,3\n')
        self.assertEqual(s['row_count'],2)
        self.assertEqual(s['column_names'],['trap_id','date','light','caught'])
        self.assertEqual(gate.tabular_summary('untrusted.zip',b'abc')['schema_status'],'OTHER_FORMAT_NOT_PARSED')

if __name__=='__main__':unittest.main()
