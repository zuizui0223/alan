"""Guard physical-unit covariance, shared controls and source selection."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / 'database/alan_experiments/plant_sources/analyze_complete_species.py'


class PenrynCompleteSpeciesTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(MODULE.exists(), 'complete-species analysis is not implemented')
        spec = importlib.util.spec_from_file_location('penryn_complete', MODULE)
        self.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.m)
        self.base = ROOT / 'database/alan_experiments/plant_sources'

    def source(self):
        return (pd.read_csv(self.base/'species_biomass_source_rows.csv'),
                pd.read_csv(self.base/'species_survival_reconciliation.csv'),
                pd.read_csv(self.base/'anic2022_community_units.csv'))

    def test_complete_rows_do_not_override_conflicting_survival(self):
        units, selection = self.m.select_units(*self.source())
        self.assertEqual(set(units.species), {
            'Agrostis tenuis', 'Anthoxanthum odoratum',
            'Hypochaeris radicata', 'Lotus pedunculatus'})
        self.assertEqual(len(units), 216)
        self.assertEqual(units.mesocosm.nunique(), 54)
        self.assertFalse(selection.set_index('species').loc['Centaurea nigra', 'included'])

    def test_duplicate_unit_species_is_rejected(self):
        d, s, c = self.source()
        with self.assertRaises(ValueError):
            self.m.select_units(pd.concat([d, d.iloc[[0]]]), s, c)

    def test_missing_required_species_cell_is_rejected(self):
        d, s, c = self.source()
        i = d.index[d.species.eq('Lotus pedunculatus')][0]
        with self.assertRaises(ValueError):
            self.m.select_units(d.drop(i), s, c)

    def test_hand_computed_species_and_shared_dark_covariance(self):
        rows = []
        # Species b is always 3 times a; all species log responses coincide.
        for light, multiplier in [('L1', 1), ('L2', 2), ('L3', 4)]:
            for i, x in enumerate([1., 2., 3.]):
                for species, factor in [('a', 1), ('b', 3)]:
                    rows.append(dict(species=species, mesocosm=f'{light}_{i}',
                                     light=light, trophic='T1', total_g=x*multiplier*factor))
        effects, cov, scores = self.m.species_effects(pd.DataFrame(rows))
        lookup = {(r.species, r.light): i for i, r in enumerate(effects.itertuples())}
        a, b, w = lookup['a', 'L2'], lookup['b', 'L2'], lookup['a', 'L3']
        self.assertAlmostEqual(effects.iloc[a].lnRR, np.log(2))
        self.assertAlmostEqual(cov[a, a], 1/6)
        self.assertAlmostEqual(cov[a, b], 1/6)
        self.assertAlmostEqual(cov[a, w], 1/12)
        v = np.zeros(len(effects)); v[a] = 1; v[b] = -1
        self.assertAlmostEqual(float(v@cov@v), 0)
        np.testing.assert_allclose(scores.T@scores, cov, atol=1e-14)

    def test_row_order_and_species_mass_units_do_not_change_effects(self):
        units, _ = self.m.select_units(*self.source())
        e, v, _ = self.m.species_effects(units)
        altered = units.sample(frac=1, random_state=11).copy()
        altered.loc[altered.species.eq('Lotus pedunculatus'), 'total_g'] *= 1000
        f, w, _ = self.m.species_effects(altered)
        self.assertEqual(e.effect_id.tolist(), f.effect_id.tolist())
        np.testing.assert_allclose(e.lnRR, f.lnRR, atol=1e-13)
        np.testing.assert_allclose(v, w, atol=1e-13)

    def test_relative_response_cancels_identical_species_responses(self):
        units, _ = self.m.select_units(*self.source())
        values = units[units.species.eq('Lotus pedunculatus')].set_index('mesocosm').total_g
        units = units.copy()
        units['total_g'] = units.mesocosm.map(values)
        effects, cov, _ = self.m.species_effects(units)
        contrasts, matrix = self.m.contrast_matrix(effects)
        np.testing.assert_allclose(matrix@effects.lnRR.to_numpy(), 0, atol=1e-14)
        np.testing.assert_allclose(matrix@cov@matrix.T, 0, atol=1e-14)
        self.assertEqual(contrasts.kind.value_counts().to_dict(),
                         {'cell': 6, 'modifier': 4, 'marginal': 2})


if __name__ == '__main__':
    unittest.main()
