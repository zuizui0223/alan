"""Original measured-source support only; no fictitious botanical labels."""
import unittest
import pandas as pd
from scripts.analysis import alan_extra_original1181146_photo_source_audit_20261108 as m

def frame(year, target, land_class, status_column):
    row=dict(target_id=target,land_class=land_class,
             photo_attempts=3,strict9_success=1,accepted_species_n=1,
             viirs_valid_months=4,viirs_median=.5,
             bio5=21.0,bio12=1134.0,elevation=82.0)
    if status_column=="population_status":
        row.update(population_status="OK",population_density_per_km2=230.)
    else:
        row.update(population_source_status="OK",population_density=189.)
    return pd.DataFrame([row]),pd.DataFrame([{"target_id":target}])

class SourceIntakeTests(unittest.TestCase):
    def test_differently_named_year_population_source_is_coalesced(self):
        baseline={
            2019:frame(2019,"2019|3_4|100","island","population_status"),
            2024:frame(2024,"2024|3_4|100","island","population_source_status"),
        }
        photos=pd.DataFrame([
          dict(unit_id="199",observation_id="a",
               target_id="2019|3_4|100",land_id="100",year=2019,
               land_class="island",source_design="bulk128_131072",
               real_strict9=True,wcvp_accepted_plant_name_id="A"),
          dict(unit_id="249",observation_id="b",
               target_id="2024|3_4|100",land_id="100",year=2024,
               land_class="island",source_design="bulk128_131072",
               real_strict9=True,wcvp_accepted_plant_name_id="B")])
        r,parts=m.summarize(photos,baseline)
        self.assertEqual(r["new_independent_observation_photos_true_strict9_and_WCVP"],2)
        self.assertEqual(r["strict9_WCVP_new_observations_with_existing_source_ALAN_exposure"]["island_photo_records"],2)
        self.assertEqual(r["new_strict9_accepted_photos_on_existing_frozen_strict_sites"],2)
        self.assertEqual(len(parts),2)

    def test_missing_exposure_not_claimed_new_measured_ALAN_result(self):
        base={2019:frame(2019,"2019|3_4|100","island","population_status")}
        q=pd.DataFrame([
          dict(unit_id="199",observation_id="a",target_id="2019|3_4|100",
               land_id="100",year=2019,land_class="island",
               source_design="bulk128_131072",real_strict9=True,
               wcvp_accepted_plant_name_id="A"),
          dict(unit_id="1991",observation_id="b",target_id="2019|9_4|111",
               land_id="111",year=2019,land_class="island",
               source_design="bias_refill_2048",real_strict9=True,
               wcvp_accepted_plant_name_id="B"),
          dict(unit_id="1992",observation_id="c",target_id="2019|9_4|111",
               land_id="111",year=2019,land_class="island",
               source_design="bias_refill_2048",real_strict9=True,
               wcvp_accepted_plant_name_id=None)])
        r,_=m.summarize(q,base)
        self.assertEqual(r["strict9_WCVP_new_observations_with_existing_source_ALAN_exposure"]["island_photo_records"],1)
        self.assertEqual(r["strict9_WCVP_new_observations_requiring_new_ALAN_exposure_or_other_env"]["island_photo_records"],1)
        self.assertTrue(r["no_supplementary_photo_retreated_as_original_historical_main"])

    def test_same_observation_extra_photo_rescues_real_unknown_not_new_observation(self):
        canonical=pd.DataFrame([
          {"observation_id":"A","photo_id":"1","inat_taxon_id":"20",
           "source_subarm":"historical_main","observed_on":"2024-05-20",
           "latitude":35.0,"longitude":139.0,"strict9":False,
           "strict9_white":float("nan")},
          {"observation_id":"B","photo_id":"2","inat_taxon_id":"21",
           "source_subarm":"historical_main","observed_on":"2024-05-20",
           "latitude":36.0,"longitude":140.0,"strict9":True,
           "strict9_white":.3},
        ])
        extra=pd.DataFrame([
          {"observation_id":"A","photo_id":"3","inat_taxon_id":"20",
           "real_strict9":True,"strict9_white":.8},
          {"observation_id":"A","photo_id":"4","inat_taxon_id":"20",
           "real_strict9":True,"strict9_white":.6},
          {"observation_id":"B","photo_id":"5","inat_taxon_id":"21",
           "real_strict9":True,"strict9_white":.4},
        ])
        rescue,report=m.candidate_observation_rescues(extra,canonical)
        self.assertEqual(report["extra_genuine_strict9_sibling_photos"],3)
        self.assertEqual(report["technical_unknown_observation_UUIDs_with_1plus_measured_sibling_photos"],1)
        self.assertEqual(report["strict9_extra_photos_for_canonical_already_classifiable_observations"],1)
        self.assertEqual(report["historical_main_observations_rescuable_original_roi9"],1)
        self.assertEqual(len(rescue),1)
        self.assertEqual(int(rescue.extra_strict9_photos_on_same_observation.iloc[0]),2)
        self.assertAlmostEqual(float(rescue.rescued_white_fraction.iloc[0]),.7)

if __name__=="__main__":
    unittest.main()
