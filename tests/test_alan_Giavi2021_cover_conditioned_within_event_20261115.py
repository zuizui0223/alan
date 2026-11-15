"""Paired floral-exposure allocation avoids fake 193 independent lighting plots."""
import unittest
import pandas as pd
from scripts.analysis.alan_Giavi2021_cover_conditioned_within_event_20261115 import (
 conditional_exposure_rows,per_site,summarize
)

def sample():
    records=[]
    for site in range(6):
        for plant in ("A","B"):
            for pair in (1,2,3):
                x={"site":str(site),"plant":plant,"pair":pair,
                  "floral_cover_T":2., "floral_cover_C":1.,
                  "all3_visits_T":2,"all3_visits_C":2}
                if site==0 and pair==1 and plant=="A":
                    x["all3_visits_T"]=0
                    x["all3_visits_C"]=0
                records.append(x)
    return pd.DataFrame(records)

class PairedOriginalExposureChecks(unittest.TestCase):
    def test_zero_visits_not_missing_and_site_total(self):
        a=conditional_exposure_rows(sample(),"all3_visits")
        self.assertEqual(len(a),36)
        self.assertEqual(int(a.zero_visit_pair.sum()),1)
        self.assertAlmostEqual(a.iloc[1].exposure_balance_lit_share,2/3)
        self.assertAlmostEqual(a.iloc[1].lit_excess_visits_beyond_flower_share,2-4*2/3)
        r,q=per_site(a,3)
        self.assertEqual(len(r),6)
        self.assertEqual(len(q),36)
        self.assertTrue(all(x["observed_minus_expected_lit_visit_fraction"]<0 for x in r))
        s=summarize(r,q,13)
        self.assertEqual(s["sites"],6)
        self.assertEqual(s["source_matched_zero_visitation_rows"],1)
        self.assertTrue(s["flower_cover_is_post_treatment_proxy_and_may_mediate_total_ALAN_effect"])
        self.assertTrue(s["visitor_count_conditions_on_recorded_matched_plant_events"])
    def test_eligibility_requires_matched_plant_occasions(self):
        a=conditional_exposure_rows(sample(),"all3_visits")
        r,q=per_site(a,4)
        self.assertEqual(len(q),0)
        self.assertEqual(len(r),0)
        s=summarize(r,q,13)
        self.assertEqual(s["sites"],0)
    def test_no_bad_physical_values(self):
        a=sample()
        a.loc[0,"floral_cover_C"]=0
        with self.assertRaises(ValueError):conditional_exposure_rows(a,"all3_visits")
        a=sample()
        a.loc[0,"all3_visits_T"]=-1
        with self.assertRaises(ValueError):conditional_exposure_rows(a,"all3_visits")
    def test_observation_missing_not_imputed_as_zero(self):
        a=sample()
        a.loc[0,"floral_cover_T"]=float("nan")
        with self.assertRaises(ValueError):conditional_exposure_rows(a,"all3_visits")
if __name__=="__main__":unittest.main()
