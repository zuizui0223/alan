"""Synthetic checks only: ALAN island-area source identity, spatial site support, moderation."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_island_area_light_20261108 as a

class IslandAreaLightTests(unittest.TestCase):
    @staticmethod
    def sample():
        rng=np.random.default_rng(71)
        rows=[]
        areas={}
        for i in range(12):
            land=str(100+i)
            areas[land]=30.0*(1.65**i)
            log_area=np.log10(areas[land])
            for site in range(3):
                for year in (2021,2024):
                    log_light=0.3+0.08*site+0.018*i+0.013*(year==2024)
                    log_pop=1.0+0.02*i+0.03*site+rng.normal(0,.04)
                    bio5=25+0.08*i+0.25*site+rng.normal(0,.05)
                    log_rain=5.0-0.01*i+0.1*site+rng.normal(0,.05)
                    elev=80+3*i+25*site+rng.normal(0,3)
                    attempts=1.8+0.08*site+rng.normal(0,.04)
                    y=(26 + 2.0*log_light - 1.6*log_light*(log_area-2.5)
                       + .2*log_pop+.09*bio5+.1*log_rain+.01*elev
                       + .15*attempts+rng.normal(0,.015))
                    rows.append(dict(target_id=f"{year}|{100+i}_{50+site}|{land}",
                        year=year,land_id=land,is_island=1,white_pp=y,
                        log_light=log_light,log_pop=log_pop,bio5=bio5,
                        log_rain=log_rain,elevation=elev,log_attempts=attempts))
        area=pd.DataFrame([dict(land_id=k,land_class="island",area_km2=v)
                           for k,v in areas.items()])
        return pd.DataFrame(rows),area

    def test_genuine_sites_not_year_repeats(self):
        data,area=self.sample()
        all_island,within,g=a.join_area_and_select(data,area)
        self.assertEqual(len(all_island),72)
        self.assertEqual(within.land_id.nunique(),12)
        self.assertTrue(g.n_sites.eq(3).all())
        only_one=data.loc[data.target_id.str.contains("_50\\|")].copy()
        _,supported,g2=a.join_area_and_select(only_one,area)
        self.assertEqual(len(supported),0)
        self.assertTrue(g2.n_sites.eq(1).all())

    def test_area_source_identity_fail_closed(self):
        data,area=self.sample()
        area.loc[0,"land_class"]="mainland"
        with self.assertRaisesRegex(ValueError,"island/area"):
            a.join_area_and_select(data,area)
        area.loc[0,"land_class"]="island"
        area.loc[0,"area_km2"]=np.nan
        with self.assertRaisesRegex(ValueError,"island/area"):
            a.join_area_and_select(data,area)

    def test_source_moderator_is_estimable_and_numeric(self):
        d,area=self.sample()
        _,within,_=a.join_area_and_select(d,area)
        centre=float(np.log10(within.groupby("land_id").area_km2.first()).median())
        model=a.fit_area_moderator(within,centre)
        self.assertEqual(model["distinct_islands"],12)
        self.assertEqual(model["cell_years"],72)
        z=model["terms"]["log_light:log10_area_c"]
        self.assertTrue(np.isfinite(z["ci_low"]) and np.isfinite(z["ci_high"]))
        self.assertLess(z["beta"],0)
        self.assertAlmostEqual(z["beta"],-1.6,delta=0.8)
        self.assertNotIn("is_island",model["model"])

if __name__=="__main__":
    unittest.main()
