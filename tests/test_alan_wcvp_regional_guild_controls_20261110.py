"""Synthetic source-guild closure tests; no ecological conclusions from mocks."""
from __future__ import annotations
import unittest
import pandas as pd
import numpy as np

from scripts.analysis import alan_wcvp_regional_guild_controls_20261110 as source


class RegionalGuildTests(unittest.TestCase):
    @staticmethod
    def sample():
        r=[]
        for island in range(40):
            for year in (2019,2021,2025):
                for site in (0,1):
                    light=.08*(island%7)+.23*site+.04*(year-2019)
                    n=100
                    moth=int(round(20+light*5))
                    butterfly=int(round(36+light*11))
                    night=int(round(8-light*1))
                    r.append({
                        "island":str(island),"year":year,"latbin":2000+island,
                        "lonbin":1000+site,
                        source.base.LIGHT:light,
                        "all_accepted_observation_uuids":10+island%3+site,
                        "all_accepted_recorded_species":100,
                        "recorded_species__native_at_TDWG3":n,
                        "recorded_moth__native_at_TDWG3":moth,
                        "recorded_butterfly__native_at_TDWG3":butterfly,
                        "recorded_night_pollen__native_at_TDWG3":night})
        return pd.DataFrame(r)

    def test_direct_moth_minus_butterfly_is_not_p_value_subtraction(self):
        d=self.sample()
        n=source.source_filter(d,"native_at_TDWG3")
        self.assertEqual(n.island.nunique(),40)
        a=source.fit_recorded_fraction(n,"recorded_moth__native_at_TDWG3",
                                     "recorded_species__native_at_TDWG3")
        b=source.fit_recorded_fraction(n,"recorded_butterfly__native_at_TDWG3",
                                     "recorded_species__native_at_TDWG3")
        n["moth_minus_butterfly"]=n.recorded_moth__native_at_TDWG3-n.recorded_butterfly__native_at_TDWG3
        c=source.fit_recorded_fraction(n,"moth_minus_butterfly",
                                     "recorded_species__native_at_TDWG3")
        self.assertAlmostEqual(a["beta_pp"]-b["beta_pp"],c["beta_pp"],places=8)
        self.assertEqual(a["n_islands"],40)

    def test_sparse_support_returns_not_estimable(self):
        d=self.sample().loc[lambda x:x.island.isin(["0","1","2"])]
        ans=source.fit_recorded_fraction(d,"recorded_moth__native_at_TDWG3",
                                       "recorded_species__native_at_TDWG3")
        self.assertEqual(ans["status"],"not_estimable")

    def test_wcvp_status_is_not_actual_wild_native(self):
        self.assertIn("native_at_TDWG3",source.STATES)
        self.assertNotIn("island_native",source.STATES)
        self.assertEqual(len(source.SOURCE_SHA),64)


if __name__=="__main__":unittest.main()
