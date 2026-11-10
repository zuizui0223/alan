"""Synthetic regression identities for real NPMS recorded-cover fixed-effects model."""
from __future__ import annotations
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.alan_npms_recorded_cover_vs_priorPCNL_plot_FE_20261110 import (
    choose,fit_source,weighted_panel_partial_out)

class NPMSTruePlantCoverFETests(unittest.TestCase):
    @staticmethod
    def source():
        rng=np.random.default_rng(41)
        rows=[]
        for plot in range(55):
            light_base=0.2+0.025*(plot%10)
            for species in ("A alba","B beta","C chrysis"):
                for year in range(2018,2024):
                    # Genuine source light changes within physical sites,
                    # but within plot the exposure is identical across taxa.
                    light=float(light_base+.05*(year-2018)+.015*(plot%3)*(year-2018))
                    site_baseline=18+(plot%7)+(ord(species[0])-65)*2
                    year_effect=(year-2018)*(ord(species[0])-65)*0.22
                    doy=float(120+plot%10+year%5+((plot*year)%7))
                    cov=float(site_baseline+1.8*light+year_effect+
                              .06*(doy-120)+rng.normal(0,.025))
                    rows.append(dict(location_id=str(plot),
                       survey_method_id="154",season="spring",species=species,
                       year=year,median_cover=cov/100.,
                       median_source_day_of_year=doy,n_visits_with_cover=1,
                       source_domin_encodings="percent_text_range",
                       source_record_status="V",original_annual_priorPCNL=light,
                       source_prior_unique_values=1,
                       recorded_cover_pp=cov,
                       recorded_cover_log=np.log1p(cov),
                       site_light=light,
                       doy_30=(doy-180)/30.,
                       log_recorded_visits=np.log(2),
                       panel_key=f"{plot}|154|spring|{species}",
                       species_year_key=f"{species}|{year}",
                    ))
        return pd.DataFrame(rows)

    def test_subset_requires_original_light_variability_and_repeat_species_positive(self):
        d=self.source()
        x,a=choose(d,source_min_range=.02)
        self.assertEqual(len(x),len(d))
        self.assertEqual(a["source_true_independent_NPMS_plots"],55)
        self.assertEqual(a["source_true_source_panel_plot_species_method_season_groups"],165)
        d.loc[d.location_id.eq("0"),"site_light"]=1.
        sub,info=choose(d)
        self.assertEqual(info["source_true_independent_NPMS_plots"],54)
        self.assertFalse((sub.location_id=="0").any())

    def test_exact_two_set_original_FE_recovers_within_site_light_slope(self):
        d=self.source()
        m=fit_source(d,"recorded_cover_pp","plot_equal","species_year",True)
        self.assertEqual(m["status"],
                         "estimated_exploratory_within_positive_recorded_plant_cover")
        self.assertEqual(m["actual_independent_source_plots"],55)
        self.assertAlmostEqual(m["terms"]["site_light"]["beta"],1.8,delta=.50)

    def test_semantics_conditioned_positive_plant_observation_not_absence(self):
        d=self.source()
        miss=d.drop(d.index[:2])
        q,a=choose(miss)
        self.assertGreater(a["source_true_source_panel_plot_species_method_season_groups"],100)
        self.assertTrue(q.median_cover.between(0,1).all())
        self.assertTrue((q.n_visits_with_cover>0).all())

if __name__=="__main__":
    unittest.main()
