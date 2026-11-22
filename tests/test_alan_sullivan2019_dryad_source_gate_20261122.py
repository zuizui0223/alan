import importlib.util
import unittest
from pathlib import Path

modpath=Path(__file__).resolve().parents[1]/'scripts/acquisition/alan_sullivan2019_dryad_experiment_source_gate.py'
spec=importlib.util.spec_from_file_location('sullivan_source',modpath)
S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)

class TestSullivanSourceGuards(unittest.TestCase):
    def test_two_independent_response_strata(self):
        self.assertIn('WetlandFieldExperiment_2013',S.EXPECTED[0]);self.assertEqual(len(S.EXPECTED),2)
        self.assertNotIn('UrbanStream',S.EXPECTED[0]);self.assertTrue(S.PAPER.startswith('10.1002/'))
    def test_disallow_untrusted_url(self):
        with self.assertRaises(ValueError):S.fetch('https://example.com/api/v2/files/1/download')
    def test_disallow_ftp(self):
        with self.assertRaises(ValueError):S.fetch('ftp://datadryad.org/api/v2/files/1/download')
    def test_readme_identifies_light_fields(self):
        x=S.summarize_text(b'Light treatment: lit or control\nPlot level: 10\nrepeated nights\n','README_for_Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt')
        self.assertEqual(x['format'],'README_text');self.assertEqual(len(x['relevant_lines_first_15']),2)
    def test_tsv_schema_and_groups(self):
        s=b'Plot\tLight\tYear\tCount\nA\tlit\t2013\t4\nA\tlit\t2013\t5\nB\tdark\t2013\t0\n'
        a=S.summarize_text(s,'Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt')
        self.assertEqual(a['nonempty_rows'],3);self.assertEqual(a['physical_unit_candidate_fields']['Light']['distinct'],2)
        self.assertEqual(a['physical_unit_candidate_fields']['Plot']['distinct'],2)
    def test_unavailable_source_not_interpreted_as_zero_effect(self):
        self.assertEqual(S.file_summaries({'_embedded':{'stash:files':[]}}),[])

if __name__=='__main__':unittest.main()
