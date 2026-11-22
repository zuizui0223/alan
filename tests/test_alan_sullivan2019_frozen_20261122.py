"""Check authentic 2013 wetland original-source derivatives and fail-closed counts."""
import csv, json, math, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'database/alan_experiments/arthropod_expansion_20261122/sullivan2019_source_reanalysis'
class TestSullivanFrozen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen=json.loads((BASE/'dryad_provenance_and_frozen_block_results.json').read_text())
        cls.source=next(f for f in cls.frozen['source_files'] if f['name']=='Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt')
        cls.effects=cls.source['provisional_block_contrasts']
    def test_authentic_source_identity(self):
        self.assertEqual(self.source['sha256'],'385444ab318ec4b15edd3161c740aa7f4b2d9b973fe4c2ddf6c29cc56cdec484')
        self.assertRegex(self.frozen['archive_source']['sha256'],r'^[0-9a-f]{64}$')
        self.assertEqual(self.source['experimental_count_grain']['source_rows'],1161)
        self.assertEqual(self.source['experimental_count_grain']['unique_original_station_names'],15)
    def test_3_blocks_one_programme_only(self):
        for s,d in self.effects['results'].items():
            self.assertEqual(d['station_n'],15)
            self.assertEqual(d['independent_physical_block_n'],3)
            self.assertEqual(d['matched_nights_by_block'],[4,4,4])
        self.assertEqual(self.frozen['new_meta_eligible_independent_programme_k'],0)
    def test_csv_effects_recompute_from_frozen_blocks(self):
        with (BASE/'experimental_endpoint_lnRR_candidates.csv').open() as f:rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),12)
        for row in rows:
            d=self.effects['results'][row['endpoint']]['estimands'][row['lamp_treatment']]
            self.assertAlmostEqual(float(row['lnRR']),sum(d['raw_block_lnRR'])/3,places=12)
            self.assertAlmostEqual(float(row['RR']),math.exp(float(row['lnRR'])),places=12)
            self.assertEqual(row['adoption_status'],'candidate_not_final_meta')
    def test_no_causal_mediation_misreporting(self):
        with (BASE/'within_experiment_cross_endpoint_ratios.csv').open() as f:rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),8)
        for r in rows:
            vals=self.effects['paired_cross_stratum_descriptive_ratios'][r['comparison_response']][r['lamp_treatment']]
            self.assertAlmostEqual(float(r['paired_ratio_of_RR']),vals['ratio_of_original_response_ratios'],places=12)
            self.assertEqual(r['status'],'posthoc_descriptive_not_causal_mediation')
    def test_current_run_fits_same_original_source_when_present(self):
        p=ROOT/'results/alan_sullivan2019_dryad_source_gate_20261122.json'
        if not p.exists():return
        run=json.loads(p.read_text())
        self.assertEqual(run['status'],'SOURCE_BYTES_READ_FOR_SCHEMA_AUDIT')
        rs=next(f for f in run['source_files'] if f['name']==self.source['name'])
        self.assertEqual(rs['sha256'],self.source['sha256'])
        for sm,d in self.effects['results'].items():
            for arm,v in d['estimands'].items():
                self.assertAlmostEqual(v['RR'],rs['provisional_block_contrasts']['results'][sm]['estimands'][arm]['RR'],places=12)
if __name__=='__main__':unittest.main()
