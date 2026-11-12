"""Synthetic numerical/source contract tests, NOT independent field observations."""
import sys,unittest
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import alan_13year_guild_photo_effort_20261112 as m

def example():
    rng=np.random.default_rng(500)
    data=[]
    for island in range(15):
        for year in (2019,2020,2021,2022):
            for cell in range(4):
                nrep=int(rng.integers(2,15))
                nins=int(rng.integers(6,45))
                ng=int(rng.integers(0,nrep+1))
                ni=int(rng.integers(1,nins))
                L=.3+.1*(year-2019)+.4*cell+float(rng.normal(0,.05))
                data.append(dict(
                    island=island,year=year,
                    islandyear=f"i{island}|{year}",site=f"i{island}|cell{cell}",
                    loglight=L,gecko_share=ng/nrep,
                    moth_share=100*ni/nins,
                    gecko_source_photo_count=ng,
                    original_total_reptile_photo_source_events=nrep,
                    all_four_family_records=nins,
                    log_insect=np.log1p(nins),
                    log_reptile=np.log1p(nrep),
                    latitude=22+.2*island+.1*cell))
    return pd.DataFrame(data)

class SourcePhotoEffortInvariants(unittest.TestCase):
    def test_count_weighted_spatial_FWL_matches_independent_WLS(self):
        d=example()
        for r,w in (("gecko","island_equal"),
                    ("gecko","reptile_capped_count"),
                    ("gecko_jeffreys","reptile_capped_count"),
                    ("moth","insect_capped_count")):
            result=m.fit_photo_share(d,r,w,repeat=False)
            self.assertTrue(result["estimable"])
            d["outcome"]=100*(d.gecko_share if r=="gecko" else (
               (d.gecko_source_photo_count+.5)/
               (d.original_total_reptile_photo_source_events+1)) if
               r=="gecko_jeffreys" else d.moth_share/100)
            counts=1/d.groupby("island").island.transform("size").to_numpy(float)
            if w=="reptile_capped_count":
                counts*=np.minimum(d.original_total_reptile_photo_source_events,10)
            if w=="insect_capped_count":
                counts*=np.minimum(d.all_four_family_records,40)
            reference=smf.wls("outcome ~ loglight + log_insect + log_reptile + latitude + C(islandyear)",
                               data=d,weights=counts).fit()
            self.assertAlmostEqual(result["effect"]["beta"],
                                   reference.params["loglight"],places=8)

    def test_two_way_FE_matches_explicit_source_site_WLS(self):
        d=example()
        result=m.fit_photo_share(d,"gecko",weight="reptile_capped_count",repeat=True)
        self.assertTrue(result["estimable"])
        self.assertEqual(result["fe_components"],15)
        self.assertEqual(result["fixed_effects_rank"],60+60-15)
        d["outcome"]=100*d.gecko_share
        weight=(1/d.groupby("island").island.transform("size").to_numpy(float)
                *np.minimum(d.original_total_reptile_photo_source_events,10))
        reference=smf.wls("outcome ~ loglight + log_insect + log_reptile + C(site) + C(islandyear)",
                           data=d,weights=weight).fit()
        self.assertAlmostEqual(result["effect"]["beta"],
                               reference.params["loglight"],places=7)

    def test_source_gate_requires_an_observed_light_range_in_the_same_islandyear(self):
        d=example()
        d.loc[d.islandyear.eq("i0|2019"),"loglight"]=.2
        result=m.gate(d,1,False)
        self.assertFalse((result.islandyear=="i0|2019").any())
        self.assertTrue((result.islandyear=="i1|2019").any())

    def test_reliability_does_not_infer_biological_absence(self):
        d=example()
        res=m.fit_photo_share(d,"gecko_jeffreys")
        self.assertTrue(res["estimable"])
        self.assertTrue(np.isfinite(res["effect"]["beta"]))
        self.assertFalse(m.fit_photo_share(d.loc[d.island.lt(10)],"gecko")["estimable"])

if __name__=="__main__":
    unittest.main()
