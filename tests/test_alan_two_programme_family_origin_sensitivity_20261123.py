"""Unit and frozen-source tests; original author data fetched only in GitHub CI replay."""
import importlib.util
import math
import unittest
from pathlib import Path

SRC=(Path(__file__).resolve().parents[1] / 'scripts/analysis/alan_two_programme_family_origin_sensitivity_20261123.py')
S=importlib.util.spec_from_file_location('taxon_sensitivity',SRC)
A=importlib.util.module_from_spec(S);S.loader.exec_module(A)


def fixture():
    rows=[]
    for i in range(4):
        rows.append({'programme':'Parkinson2016','colour':'lit' if i<2 else 'dark',
                     'habitat':'littoral','unit':str(i),'aquatic':(18+i) if i<2 else (5+i),
                     'terrestrial':(5+i),'family':{'ephydridae':1,'caenidae':5,'leptoceridae':1}})
    for h in ('littoral','riparian'):
        for colour in ('green','purple','red','white','dark'):
            for i in range(3):
                rows.append({'programme':'SevenPonds2022','colour':colour,'habitat':h,
                            'unit':h+'_'+colour+str(i),'aquatic':(50+i) if colour!='dark' else (10+i),
                            'terrestrial':10+i,'family':{'ephydridae':1,'caenidae':5,'leptoceridae':1}})
    return rows

class OriginSensitivityTests(unittest.TestCase):
    def test_known_git_blob_sha1(self):
        self.assertEqual(A.git_blob_sha(b'test'), '30d74d258442c7c65512eafab474568dd706c430')
    def test_study_family_harmonization_preserves_physical_units(self):
        rows=fixture()
        adj=A.excluded_rows(rows,('ephydridae','caenidae','leptoceridae'))
        self.assertEqual(len(adj),len(rows))
        self.assertEqual(adj[0]['aquatic'],rows[0]['aquatic']-7)
        self.assertEqual(adj[0]['terrestrial'],rows[0]['terrestrial'])
        self.assertEqual(adj[4]['terrestrial'],rows[4]['terrestrial']-1)
        self.assertEqual(adj[4]['aquatic'],rows[4]['aquatic']-6)
        self.assertTrue(all(a['unit']==b['unit'] for a,b in zip(adj,rows)))
    def test_two_programmes_not_two_colour_cells(self):
        rows=fixture();p=A.one_programme(rows,'Parkinson2016');q=A.one_programme(rows,'SevenPonds2022')
        self.assertEqual(p['arms']['lit']['physical_light_units'],2)
        self.assertEqual(q['arms']['lit']['physical_light_units'],24)
        self.assertEqual(q['arms']['dark']['physical_light_units'],6)
        self.assertEqual(q['arms']['lit']['treatment_cells'],8)
        self.assertEqual(q['arms']['dark']['treatment_cells'],2)
    def test_colour_removal_preserves_dark_control_cells(self):
        d=fixture();p=A.one_programme(d,'SevenPonds2022')
        v=A.one_programme(d,'SevenPonds2022',omit_colour='purple')
        self.assertEqual(v['arms']['lit']['treatment_cells'],6)
        self.assertEqual(v['arms']['dark']['treatment_cells'],2)
        self.assertEqual(v['arms']['dark'],p['arms']['dark'])
    def test_sample_cov_of_identical_series(self):
        self.assertAlmostEqual(A.sample_cov([1,3,5],[1,3,5]),4)
        self.assertAlmostEqual(A.sample_cov([1,3,5],[5,3,1]),-4)
    def test_k_two_remle_closed_form_and_wide_t1(self):
        x={'log_ratio':math.log(5),'sampling_variance':0.1}
        y={'log_ratio':math.log(12),'sampling_variance':0.35}
        t=A.two_programme_mkh(x,y)
        self.assertEqual(t['programmes'],2)
        self.assertAlmostEqual(t['tau2_REML_two_programme_closed_form'],max(0,((math.log(5)-math.log(12))**2-0.45)/2))
        self.assertLess(t['lower'],1)
        self.assertGreater(t['upper'],1)
    def test_unavailable_taxon_is_not_silently_reclassified(self):
        with self.assertRaises(ValueError):A.excluded_rows(fixture(),('not_a_verified_family',))
    def test_nonlinear_two_study_pool_not_arithmetic_average(self):
        x={'log_ratio':math.log(5),'sampling_variance':0.1}
        y={'log_ratio':math.log(12),'sampling_variance':0.35}
        pooled=A.two_programme_mkh(x,y)['ratio']
        self.assertNotAlmostEqual(pooled,(5+12)/2)
    def test_author_source_commit_pinned_and_bounded(self):
        self.assertIn(A.SOURCE_COMMIT,A.AUTHOR_SOURCE_URL)
        self.assertEqual(A.GIT_OBJECTS['parkinson_2020_authors'],'4e2b7065d67e0604ca121ea26c96b5db6d64de94')
    def test_mkh_no_pretense_of_many_studies(self):
        self.assertTrue(all(len(fam)>=1 for _,fam in A.PRIMARY_VARIANTS))
        self.assertEqual(len(A.COLOURS),4)

if __name__=='__main__':unittest.main()
