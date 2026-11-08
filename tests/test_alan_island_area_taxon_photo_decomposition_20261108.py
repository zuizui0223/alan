"""Synthetic contracts for the ALAN original-photo community vs photo deviation split."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import alan_island_area_taxon_photo_decomposition_20261108 as m

SPEC=m.taxa.SPECIES_COL

def sample():
    rng=np.random.default_rng(9301)
    rows=[]
    for i in range(12):
        land=str(700+i)
        lga=np.log10(30*1.8**i)
        lgd=np.log10(1+2+(i*7)%55)
        for year in (2019,2024):
            for site in range(3):
                light=0.27+0.08*site+.012*i+.025*(year==2024)
                pop=.7+.05*site+.02*i+rng.normal(0,.01)
                temp=18+.3*i+.4*site+rng.normal(0,.05)
                rain=5.7+.08*i-.04*site+rng.normal(0,.05)
                elev=150+7*i+5*site+rng.normal(0,1)
                attempt=1.8+.08*site+rng.normal(0,.02)
                ref=28+1.4*light*lga +.8*pop+.1*temp+rng.normal(0,.04)
                deviation=1.5+.9*light*lgd +.1*rain+rng.normal(0,.04)
                rows.append({
                  "land_id":land,"area_km2":10**lga,
                  "distance_mainland_vertex_km":10**lgd-1,
                  "year":year,"target_id":f"{year}|{i}_{site}|{land}",
                  "white_pp":ref+deviation,"log_light":light,"log_pop":pop,
                  "bio5":temp,"log_rain":rain,"elevation":elev,
                  "log_attempts":attempt,
                  "reference":ref,"deviation":deviation,
                })
    return pd.DataFrame(rows)

class DecompositionTests(unittest.TestCase):
    def test_source_reference_join_exact_real_species_gates(self):
        d=pd.DataFrame([
          {"target_id":"2019|a|9","white_pp":40,"accepted_species_n":2},
          {"target_id":"2019|b|9","white_pp":25,"accepted_species_n":1}])
        s=pd.DataFrame([
          {"target_id":"2019|a|9",SPEC:"A"},
          {"target_id":"2019|a|9",SPEC:"B"},
          {"target_id":"2019|b|9",SPEC:"C"}])
        t=pd.DataFrame([
          {"target_id":"2019|a|9","taxon_average_white_fraction":.3,
           "source_singleton_taxa":1,"taxon_species_count":2},
          {"target_id":"2019|b|9","taxon_average_white_fraction":.4,
           "source_singleton_taxa":0,"taxon_species_count":1}])
        x=m.source_components(d,s,t)
        self.assertTrue(np.allclose(x.white_pp,x.taxon_reference_white_pp+x.taxon_relative_photo_pp))
        self.assertEqual(list(x.taxon_relative_photo_pp),[10.,-15.])
        t.loc[t.target_id.eq("2019|a|9"),"taxon_species_count"]=3
        with self.assertRaisesRegex(ValueError,"species"):
            m.source_components(d,s,t)

    def test_exact_same_rows_area_and_distance_ols_and_equal_island(self):
        a=sample()
        c1=np.log10(a.groupby("land_id").area_km2.first()).median()
        c2=np.log10(1+a.groupby("land_id").distance_mainland_vertex_km.first()).median()
        for eq in (False,True):
            r=m.exact_decomposition_model(
                a,"white_pp","reference","deviation",c1,c2,eq)
            self.assertTrue(r["both_interactions_reconstruct_exactly"])
            for term in m.TERMS:
                o=r["results"]["observed"]["terms"][term]["beta"]
                u=r["results"]["taxon_photo_reference"]["terms"][term]["beta"]
                v=r["results"]["taxon_relative_photo_deviation"]["terms"][term]["beta"]
                self.assertAlmostEqual(o,u+v,places=7)
                self.assertTrue(np.isfinite(o))
            self.assertEqual(r["results"]["observed"]["unique_islands"],12)

    def test_prior_photo_reference_never_uses_holdout_year(self):
        s=pd.DataFrame([
          {"target_id":"2019|old|9","year":2019,SPEC:"A","photo_white_fraction":.2},
          {"target_id":"2021|old|9","year":2021,SPEC:"B","photo_white_fraction":.8},
          {"target_id":"2024|a|9","year":2024,SPEC:"A","photo_white_fraction":.9},
          {"target_id":"2024|a|9","year":2024,SPEC:"B","photo_white_fraction":.5},
          {"target_id":"2025|b|9","year":2025,SPEC:"A","photo_white_fraction":.7},
          {"target_id":"2025|b|9","year":2025,SPEC:"B","photo_white_fraction":.5},
        ])
        cells=pd.DataFrame([
          {"target_id":"2024|a|9","year":2024,"land_id":"9",
           "accepted_species_n":2,"log_light":0.4},
          {"target_id":"2025|b|9","year":2025,"land_id":"9",
           "accepted_species_n":2,"log_light":0.8}])
        actual,receipt=m.earlier_year_reference_holdout(s,cells)
        self.assertEqual(receipt["test_cell_years_having_2plus_genuinely_prior_photographed_species"],2)
        self.assertEqual(receipt["test_islands_having_2plus_true_different_sites"],1)
        self.assertEqual(len(actual),2)
        self.assertAlmostEqual(actual.loc[actual.year.eq(2024),"prior_photo_reference_pp"].iloc[0],50)
        self.assertAlmostEqual(actual.loc[actual.year.eq(2024),"current_photo_white_pp"].iloc[0],70)
        self.assertTrue(receipt["holdout_photo_years_not_used_to_define_historical_taxon_reference"])
        s.loc[s.year.eq(2024)&s[SPEC].eq("A"),"photo_white_fraction"]=.1
        other,_=m.earlier_year_reference_holdout(s,cells)
        self.assertAlmostEqual(other.loc[other.year.eq(2024),"prior_photo_reference_pp"].iloc[0],50)

if __name__=="__main__":
    unittest.main()
