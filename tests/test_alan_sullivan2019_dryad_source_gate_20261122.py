import importlib.util
import unittest
import io, zipfile
from pathlib import Path

modpath=Path(__file__).resolve().parents[1]/'scripts/acquisition/alan_sullivan2019_dryad_experiment_source_gate.py'
spec=importlib.util.spec_from_file_location('sullivan_source',modpath)
S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)

class TestSullivanSourceGuards(unittest.TestCase):
    def test_two_independent_response_strata(self):
        self.assertIn('WetlandFieldExperiment_2013',S.EXPECTED[0]);self.assertEqual(len(S.EXPECTED),2)
        self.assertNotIn('UrbanStream',S.EXPECTED[0]);self.assertTrue(S.PAPER.startswith('10.1002/'))
    def test_official_dryad_relative_hal_links_are_resolved_only_inside_host(self):
        source={'_links':{'stash:version':{'href':'/api/v2/versions/8747'},
                          'malicious':{'href':'https://evil.example/api/v2/files/1'},
                          'self':{'href':'/api/v2/datasets/doi%3A10.5061%2Fdryad.gb9534p'}}}
        out=S.links(source)
        self.assertIn('https://datadryad.org/api/v2/versions/8747',out)
        self.assertTrue(all(x.startswith('https://datadryad.org/api/v2/') for x in out))
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
    def test_only_original_2013_zip_members_enter_the_source_audit(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED) as z:
            z.writestr('README_for_Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt','Light: lit or dark\nPlot: 1\n')
            z.writestr('Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt','Plot\tLight\tCount\nA\tlit\t10\nB\tdark\t5\n')
            z.writestr('Sullivan_EAP18-0555_UrbanStreamGradient_2010-2011.txt','DO NOT MIX')
        files=S.unpack_expected_dataset_zip(b.getvalue())
        self.assertEqual(len(files),2)
        self.assertEqual(set(x['status'] for x in files),{'SOURCE_BYTES_VERIFIED'})
        self.assertEqual([x for x in files if x['name'].startswith('Sullivan_')][0]['physical_unit_candidate_fields']['Light']['distinct'],2)
    def test_grouped_experiment_only_matches_same_group_and_date_not_missing_zeros(self):
        data=('Station\tTreatment\tIntensity\tSpectra\tGroup\tSample\tDate\tTaxa\tCount\n'
              'X_DRK\tDRK\t0.1\tsky\tG1\tAquaticEmergent\t9/1/2013\tA\t2\n'
              'X_LED\tLEDI\t4\tLED\tG1\tAquaticEmergent\t9/1/2013\tA\t5\n'
              'X_LED\tLEDI\t4\tLED\tG1\tAquaticEmergent\t9/2/2013\tA\t0\n'
              'Y_DRK\tDRK\t0.1\tsky\tG2\tAquaticEmergent\t9/1/2013\tA\t1\n')
        z=S.audit_grouped_experiment_rows(data.encode())
        self.assertEqual(z['source_rows'],4)
        self.assertEqual(z['unique_original_station_names'],3)
        self.assertEqual(z['source_count_explicit_zero_rows'],1)
        self.assertEqual(z['by_sample']['AquaticEmergent']['matched_group_date_pairs']['LEDI'],1)
        self.assertEqual(z['by_sample']['AquaticEmergent']['fully_observed_five_treatment_group_dates'],0)
    def test_three_true_blocks_yield_four_correlated_dark_control_contrasts(self):
        rows=['Station\tTreatment\tGroup\tSample\tDate\tCount']
        for sm in ('AquaticEmergent','TerrestrialArthropod','WebCount'):
            for g in ('1','2','3'):
                for day in ('d1','d2'):
                    for t,n in (('DRK',1),('LEDI',2),('LEDH',3),('HPSI',4),('HPSH',5)):
                        rows.append(f'{g}_{t}\t{t}\t{g}\t{sm}\t{day}\t{n}')
        out=S.within_block_capture_contrasts(('\n'.join(rows)+'\n').encode())
        self.assertEqual(out['status'],'SOURCE_DERIVED_DESCRIPTIVE_BLOCK_ESTIMATES_ONLY')
        self.assertEqual(out['results']['AquaticEmergent']['independent_physical_block_n'],3)
        self.assertEqual(out['results']['AquaticEmergent']['station_n'],15)
        self.assertAlmostEqual(out['results']['AquaticEmergent']['estimands']['LEDI']['RR'],2.0)
        self.assertAlmostEqual(out['results']['AquaticEmergent']['shared_dark_covariance_log_scale']['LEDI']['HPSI'],0.0)
    def test_unavailable_source_not_interpreted_as_zero_effect(self):
        self.assertEqual(S.file_summaries({'_embedded':{'stash:files':[]}}),[])

if __name__=='__main__':unittest.main()
