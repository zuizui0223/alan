import csv
import importlib.util
import math
from pathlib import Path
import unittest
import numpy as np

DB = Path(__file__).resolve().parents[1]/'database/alan_experiments'
spec = importlib.util.spec_from_file_location('effect_calculator', DB/'compute_effects.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ExperimentDatabase(unittest.TestCase):
    def setUp(self):
        with (DB/'effects.csv').open(encoding='utf-8') as f:
            self.rows = list(csv.DictReader(f))

    def test_hand_calculation_and_se_conversion(self):
        row = dict(self.rows[0],mean_lit='4',mean_dark='2',dispersion_lit='2',dispersion_dark='1',n_lit='10',n_dark='10')
        answer = module.compute(row)
        self.assertAlmostEqual(answer['lnRR'], math.log(2))
        self.assertAlmostEqual(answer['sampling_variance'], .05)
        row.update(dispersion_type='SE',dispersion_lit=str(2/math.sqrt(10)),dispersion_dark=str(1/math.sqrt(10)))
        self.assertAlmostEqual(module.compute(row)['sampling_variance'], .05)

    def test_nonestimable_designs_fail(self):
        for changes in ({'mean_dark':'0'},{'dispersion_unit':'individual'}, {'design':'paired'},{'n_lit':'1'},{'mean_lit':'nan'}):
            with self.assertRaises(ValueError):
                module.compute(dict(self.rows[0], **changes))

    def test_covariance_and_source_unit(self):
        ids = [r['effect_id'] for r in self.rows]
        matrix = np.zeros((len(ids),len(ids)))
        with (DB/'effect_covariance.csv').open(encoding='utf-8') as f:
            entries = list(csv.DictReader(f))
        self.assertEqual(len(entries),36)
        for r in entries:
            matrix[ids.index(r['effect_id_1']),ids.index(r['effect_id_2'])]=float(r['sampling_covariance'])
        np.testing.assert_allclose(matrix,matrix.T)
        self.assertGreaterEqual(np.linalg.eigvalsh(matrix).min(),-1e-12)
        np.testing.assert_allclose(matrix.diagonal(),[module.compute(r)['sampling_variance'] for r in self.rows])
        self.assertEqual(len(set(r['experiment_id'] for r in self.rows)),1)
        self.assertEqual(sum(r['analysis_role']=='primary' for r in self.rows),4)
        self.assertAlmostEqual(float(self.rows[0]['mean_lit'])/float(self.rows[0]['mean_dark']),163.7/24.7)


if __name__ == '__main__':
    unittest.main()
