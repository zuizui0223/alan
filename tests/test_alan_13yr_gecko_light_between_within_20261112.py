"""Independent synthetic same-site PCNL within/between source regression tests."""
import sys,unittest
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
sys.path.insert(0,str(Path(__file__).parents[1]/"scripts"/"analysis"))
import alan_13yr_gecko_light_between_within_20261112 as m

def src():
    rng=np.random.default_rng(753)
    rows=[]
    for island in range(15):
        for cell in range(4):
            for year in (2019,2020,2021,2022):
                n=int(rng.integers(2,12))
                k=int(rng.integers(0,n+1))
                L=.4+.25*cell+.08*(year-2019)+rng.normal(0,.04)
                rows.append(dict(
                  island=island,year=year,
                  site=f"i{island}|cell{cell}",islandyear=f"i{island}|{year}",
                  loglight=L,
                  log_insect=float(rng.uniform(1,6)),
                  log_reptile=float(np.log1p(n)),
                  latitude=30+island/2+.1*cell,
                  gecko_source_photo_count=k,
                  original_total_reptile_photo_source_events=n,
                  gecko_share=k/n,
                  moth_share=float(rng.uniform(10,90))
                ))
    return pd.DataFrame(rows)

class Decomposition(unittest.TestCase):
    def test_correlated_pair_equals_independent_statsmodels_WLS(self):
        d=src()
        d["meanlight"]=d.groupby("site").loglight.transform("mean")
        d["withinlight"]=d.loglight-d.meanlight
        for q in ("insect","reptile"):
            d["mean_"+q+"_effort"]=d.groupby("site")["log_"+q].transform("mean")
            d["within_"+q+"_effort"]=d["log_"+q]-d["mean_"+q+"_effort"]
        for r in ("gecko","gecko_smoothed","moth"):
            o=m.fit(d,r)
            self.assertEqual(o["original_GSHHG_islands"],15)
            d["response"]=100*d.gecko_share if r=="gecko" else (
                100*(d.gecko_source_photo_count+.5)/(d.original_total_reptile_photo_source_events+1)
                if r=="gecko_smoothed" else d.moth_share)
            w=1/d.groupby("island").island.transform("size").to_numpy(float)
            fit=smf.wls("response ~ meanlight + withinlight + mean_insect_effort + within_insect_effort + mean_reptile_effort + within_reptile_effort + latitude + C(islandyear)",data=d,weights=w).fit()
            self.assertAlmostEqual(o["effects"]["between_site_meanlight"]["beta"],fit.params["meanlight"],places=8)
            self.assertAlmostEqual(o["effects"]["within_site_light_change"]["beta"],fit.params["withinlight"],places=8)
            diff=o["effects"]["within_minus_between"]["beta"]
            self.assertAlmostEqual(diff,fit.params["withinlight"]-fit.params["meanlight"],places=8)

    def test_count_weight_still_true_original_observation_unit(self):
        d=src()
        o=m.fit(d,"gecko","reptile_count")
        self.assertEqual(o["source_cell_years"],240)
        self.assertEqual(o["sites"],60)
        self.assertTrue(np.isfinite(o["effects"]["within_minus_between"]["beta"]))

if __name__=="__main__":unittest.main()
