"""Unit guards for P1 real Lilium night/day contrast; no pseudo-replication."""
import math,unittest
from scripts.analysis.alan_p1_lilium_night_day_crossstudy_20261115 import summarize

def row(study,year,ratio,log_sd=.2,day=.2):
    return dict(study_ID=study,year=year,
        daytime_seed_ovule=day,nighttime_seed_ovule=day*ratio,
        night_day_ratio=ratio,log_night_day=math.log(ratio),
        approximate_independent_arm_sampling_se=log_sd)

class TestLiliumNightDay(unittest.TestCase):
    def test_actual_ratio_sign_and_two_source_papers_only(self):
        v=[row("Morinaga_2009",2003,.8),
           row("Morinaga_2009",2004,1),
           row("Morinaga_2009",2005,1.25),
           row("Nakajima_2018",2016,16)]
        x=summarize(v)
        self.assertEqual(x["independent_publications"],2)
        self.assertTrue(x["years_are_not_three_independent_mainland_studies"])
        self.assertAlmostEqual(x["night_day_ratio_of_ratios_island_over_mainland"],16)
        self.assertFalse(x["causal_island_adaptation_established"])
        self.assertFalse(x["ALAN_exposure_measured"])
        self.assertFalse(x["night_dependence_from_day_over_open_estimable_in_mainland"])

    def test_swapped_island_day_and_night_changes_sign(self):
        main=[row("Morinaga_2009",k,1) for k in (2003,2004,2005)]
        pos=summarize(main+[row("Nakajima_2018",2016,5)])
        neg=summarize(main+[row("Nakajima_2018",2016,1/5)])
        self.assertGreater(pos["island_minus_mainland_log_night_day_equal_year"],0)
        self.assertLess(neg["island_minus_mainland_log_night_day_equal_year"],0)

    def test_author_years_not_independent_publications(self):
        with self.assertRaises(ValueError):
            summarize([row("Morinaga_2009",2003,1),row("Nakajima_2018",2016,2)])

    def test_sampling_only_interval_not_causal(self):
        v=[row("Morinaga_2009",k,1) for k in (2003,2004,2005)]
        v.append(row("Nakajima_2018",2016,10,day=.03))
        x=summarize(v)
        lo,hi=x["sampling_only_approximate_ratio95_WALD_not_full_CI"]
        self.assertLess(lo,10);self.assertGreater(hi,10)
        self.assertFalse(x["source_figure_extraction_error_quantified"])
        self.assertFalse(x["same_study_island_mainland_paired"])
        self.assertEqual(x["island_day_vs_night_access_hours"],9)
        self.assertEqual(x["mainland_day_vs_night_access_hours"],12)

if __name__=="__main__":unittest.main()
