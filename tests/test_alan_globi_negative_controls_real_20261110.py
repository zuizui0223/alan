"""Source-locked negative-control semantics; synthetic values are ONLY code tests."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from scripts.analysis import alan_globi_negative_controls_real_20261110 as audit
from scripts.analysis import alan_full93m_moth_PCNLPrior_two_way_FE_20261109 as base


class TestNegativeControlSemantics(unittest.TestCase):
    def synthetic(self):
        rows=[]
        for island in range(40):
            for year in (2019,2023,2025):
                for cell in (0,1):
                    light=.01+.035*(island%7)+.4*cell+.012*(year-2019)
                    both=int(round(5000*(.3+.01*light)))
                    m_only=int(round(5000*(.04-.003*light)))
                    b_only=int(round(5000*(.2+.004*light)))
                    all_species=5000
                    d=dict(island=str(island),latbin=100+island,lonbin=200+cell,
                           year=year,all_accepted_recorded_species=all_species,
                           moth_visitor_positive_recorded_species=both+m_only,
                           all_accepted_observation_uuids=8+island%3+cell,
                           all_accepted_photo_ids=10+island%3+cell,
                           moth_positive_fraction_OF_RECORDED_WCVP_SPECIES=(both+m_only)/all_species,
                           pcnl_years_expected=year-1992,
                           pcnl_years_observed=year-1992,
                           historic_ntl_complete=True,
                           pcnl_1992_to_previous_year_mean_log1p=light,
                           both=both,moth=both+m_only,butterfly=both+b_only,
                           moth_record_no_butterfly=m_only,butterfly_record_no_moth=b_only,
                           neither_documented=all_species-both-m_only-b_only,
                           moth_pollinates=int(both/2),night_pollen=int(both/4))
                    rows.append(d)
        return pd.DataFrame(rows)

    def test_union_close_and_paired_difference(self):
        frame=self.synthetic()
        frame["any_moth_or_butterfly"]=frame.moth+frame.butterfly_record_no_moth
        frame["any_minus_moth"]=frame.any_moth_or_butterfly-frame.moth
        self.assertTrue(frame.any_minus_moth.eq(frame.butterfly_record_no_moth).all())
        x=base.prepare(frame)
        def effect(col):
            y=x.copy()
            y["frac_pp"]=100*y[col]/y.all_accepted_recorded_species
            z=base.fit(y,"island_equal")
            self.assertEqual(z["status"],"estimated_exploratory_observational")
            return z["terms"][base.LIGHT]["beta_pp"]
        self.assertAlmostEqual(effect("moth"),
            effect("both")+effect("moth_record_no_butterfly"),places=8)
        self.assertAlmostEqual(effect("any_moth_or_butterfly")-effect("moth"),
                               effect("any_minus_moth"),places=8)

    def test_nonempty_incomplete_guild_inputs_are_not_accepted(self):
        self.assertTrue(audit.CANONICAL_PCNL_SHA.startswith("440a"))
        self.assertIn("both",audit.SOURCE_NAMES.values())
        self.assertNotIn("bees",audit.SOURCE_NAMES.values())


if __name__=="__main__": unittest.main()
