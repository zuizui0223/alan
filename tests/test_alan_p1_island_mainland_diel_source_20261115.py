"""Test source island scope, sign orientation and refusal to invent evidence."""
import unittest
import pandas as pd
from scripts.analysis.alan_p1_island_mainland_diel_source_20261115 import classify,compare

def row(study,country,place,subtype='fruit set',unit='proportion',day=0.5,openmean=1,lat='0.0'):
    return {'study_ID':study,'country':country,'locality':place,
            'coordinates_lat':lat,'plant_species':'Test plant',
            'treatment_effectiveness_metric':subtype,'effectiveness_units':unit,
            'day_mean':day,'open_mean':openmean,'night_loss':1-day/openmean}

class TestP1Scope(unittest.TestCase):
    def test_exact_locality_not_just_country(self):
        p=classify(pd.DataFrame([
          row('McMullen_2009','Ecuador','Galapagos',lat='0.554473'),
          row('Wolff_2003','Ecuador','Podocarpus National Park',lat='-3.97'),
          row('Macgregor_2019','United Kingdom','Yorkshire',lat='53.84'),
          row('Chapurlat_2015','Sweden','Skane',lat='56.75027778')]))
        self.assertEqual(list(p.geo_class),['oceanic','continent','large_island','shelf_small'])

    def test_geocoding_warning_fails_closed(self):
        with self.assertRaises(ValueError):
            classify(pd.DataFrame([row('McMullen_2007','Ecuador','Pinta Island',lat='0.5858')]))
        a=classify(pd.DataFrame([row('McMullen_2007','Ecuador','Pinta Island',lat='-0.542222')]))
        self.assertTrue(bool(a.geo_warning.iloc[0]))
        with self.assertRaises(ValueError):
            classify(pd.DataFrame([row('Mystery_2026','Japan','Unknown island')]))

    def test_true_zero_day_not_pseudocount(self):
        a=classify(pd.DataFrame([row('Barrios_2012','Cuba','Pan de Matanzas',day=0,lat='23.03')]))
        self.assertTrue(pd.isna(a.log_day_open.iloc[0]))
        self.assertEqual(a.night_loss.iloc[0],1.0)

    def test_island_minus_mainland_direction(self):
        d=pd.DataFrame([
          row('Jaca_2020','Canary Islands','Tenerife',day=0.7),
          row('McMullen_2009','Ecuador','Galapagos',day=0.8,lat='0.554473'),
          row('Cutler_2012','Canada','Debert, Nova Scotia',day=0.2),
          row('Wolff_2003','Ecuador','Podocarpus National Park',day=0.3)])
        p=classify(d)
        z=compare(p,'fruit_proportion',seed=13)
        self.assertEqual(z['island_studies'],2)
        self.assertEqual(z['mainland_studies'],2)
        # Island night loss mean .25 vs mainland .75 -> signed delta -.5
        self.assertAlmostEqual(z['exploratory_island_minus_mainland_night_loss'],-.5)
        self.assertFalse(z['confirmatory'])
        self.assertEqual(z['descriptive_paper_bootstrap95'][0]<0,True)

if __name__=='__main__':unittest.main()
