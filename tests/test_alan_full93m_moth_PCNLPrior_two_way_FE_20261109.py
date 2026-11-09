"""Synthetic identification and provenance gates for ALAN source-flora PCNL model."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_full93m_moth_PCNLPrior_two_way_FE_20261109 as model


def synthetic(n_islands=40,reverse=False):
    rows=[]
    years=[2019,2023,2025]
    for island in range(n_islands):
        for year_idx,year in enumerate(years):
            for site in (0,1):
                obs=5+site+island%3
                light=0.2+0.03*(island%8)+.23*site+.07*year_idx
                frac=(36 + 2*(island%7)-5*light
                      + 1.8*np.log1p(obs)+.5*year_idx) /100
                if reverse:
                    frac=(100-frac*100)/100
                n_species=5000
                moth=int(round(n_species*frac))
                rows.append(dict(
                  island=f"island-{island}",latbin=1000+island*2+site,
                  lonbin=2000+island*2+site,year=year,
                  all_accepted_recorded_species=n_species,
                  moth_visitor_positive_recorded_species=moth,
                  all_accepted_observation_uuids=obs,
                  all_accepted_photo_ids=obs+3,
                  moth_positive_fraction_OF_RECORDED_WCVP_SPECIES=moth/n_species,
                  historic_ntl_complete=True,
                  pcnl_years_expected=year-1992,
                  pcnl_years_observed=year-1992,
                  pcnl_1992_to_previous_year_mean_log1p=light))
    return pd.DataFrame(rows)


class TestMacro(unittest.TestCase):
    def test_known_negative_source_coefficient_is_recovered(self):
        x=model.prepare(synthetic())
        self.assertEqual(x.island.nunique(),40)
        r=model.fit(x,"island_equal")
        self.assertEqual(r["status"],"estimated_exploratory_observational")
        self.assertAlmostEqual(r["terms"][model.LIGHT]["beta_pp"],-5.0,delta=.1)
        self.assertEqual(r["year_count"],3)

    def test_source_empty_and_nonidentifiable_are_not_fabricated(self):
        x=synthetic(n_islands=4)
        r=model.fit(model.prepare(x),"island_equal")
        self.assertEqual(r["status"],"insufficient_source_support")

    def test_missing_1992_to_prior_year_support_is_excluded(self):
        x=synthetic()
        x.loc[x.index==0,"pcnl_years_observed"]=10
        work=model.prepare(x)
        self.assertEqual(len(work),len(x)-1)
        self.assertTrue((work.pcnl_years_observed == work.year-1992).all())

    def test_double_counted_island_site_year_rejected(self):
        x=synthetic()
        with self.assertRaisesRegex(ValueError,"Double counted"):
            model.prepare(pd.concat([x,x.iloc[:1]],ignore_index=True))

    def test_genuine_sameyear_spatial_light_contrasts_required(self):
        x=synthetic()
        for island,group in x.groupby("island"):
            for year,sub in group.groupby("year"):
                x.loc[sub.index,model.LIGHT]=sub[model.LIGHT].iloc[0]
        filtered=model.prepare(x)
        self.assertEqual(len(filtered),0)

    def test_raw_photo_fraction_does_not_depend_on_FCP_ROI_pixels(self):
        x=synthetic()
        self.assertNotIn("white_pp",x)
        self.assertNotIn("strict9",x)
        r=model.fit(model.prepare(x),"site_equal")
        self.assertEqual(r["status"],"estimated_exploratory_observational")


if __name__=="__main__":
    unittest.main()
