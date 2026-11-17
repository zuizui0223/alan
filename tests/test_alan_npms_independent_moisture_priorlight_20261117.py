"""Synthetic numerical tests, no biological conclusions from synthetic fixtures."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_npms_independent_moisture_priorlight_20261117 as a

class TestIndependentNPMSWetness(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(1726)
        taxa=["plant_"+str(x) for x in range(4)]
        moisture=[2.,5.,8.,10.]
        rec=[]
        for site in range(40):
            for yr in range(2018,2023):
                light=.3+site*.003+(yr-2018)*(.04+site*.0011)
                year_src=.01*(yr-2018)
                plotyear_src=np.sin(site*.13+yr*.27)*.15
                for j,taxon in enumerate(taxa):
                    doy=(int(rng.integers(140,231))-180)/30
                    signal=36 + site*.2 + j*.8 + year_src + plotyear_src + .73*light*(moisture[j]-7)+.19*doy
                    rec.append({"location_id":str(site),"year":yr,
                                "survey_method_id":"154","season":"spring","species":taxon,
                                "prior_light":light,"wet7":float(moisture[j]>7),
                                "moist_centered":moisture[j]-7,
                                "light_moist":light*(moisture[j]-7),
                                "light_wet7":light*float(moisture[j]>7),
                                "doy30":doy,"cover_pp":signal,
                                "panel_key":str(site)+"|154|spring|"+taxon,
                                "species_year":taxon+"|"+str(yr),
                                "plotyear_method_season":str(site)+"|154|spring|"+str(yr)})
        cls.z=pd.DataFrame(rec)

    def test_synthetic_full_plotyear_fixed_effect_resolves_interaction(self):
        r=a.fit(self.z,["light_moist","doy30"],
                ["panel_key","species_year","plotyear_method_season"],"plot_equal")
        self.assertEqual(r["status"],"observational_positive_cover_only")
        self.assertAlmostEqual(r["coefficients"]["light_moist"]["beta"],.73,delta=.003)
        self.assertAlmostEqual(r["coefficients"]["doy30"]["beta"],.19,delta=.003)

    def test_synthetic_two_way_residualization_has_finite_coefficients(self):
        r=a.fit(self.z,["prior_light","light_moist","doy30"],
                ["panel_key","species_year"],"plot_equal")
        self.assertEqual(r["status"],"observational_positive_cover_only")
        self.assertTrue(np.isfinite(r["coefficients"]["light_moist"]["beta"]))
        self.assertEqual(r["plots"],40)

    def test_wet_and_dry_common_support_definition(self):
        m,c=a.eligibility(self.z)
        self.assertEqual(c["mixed_source_plots"],40)
        self.assertEqual(c["mixed_source_plotyear_methodseason_units"],200)
        self.assertEqual(c["wet_species"],2)
        self.assertEqual(c["drier_species"],2)
        self.assertFalse(c["pass_primary"]) # 4, not >=40 taxa

if __name__=="__main__":
    unittest.main()
