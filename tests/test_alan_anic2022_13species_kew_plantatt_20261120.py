"""Source pin and no post-hoc pseudo-WET promotion for 13 meadow species."""
import io,unittest
from scripts.analysis.alan_anic2022_13species_kew_plantatt_20261120 import (
  git_blob_hash, trait_lookup, describe, infer, rho
)
class TestPLANTATT(unittest.TestCase):
  def test_git_blob_digest_matches_known_definition(self):
    import hashlib
    b=b"this is a check\n"
    self.assertEqual(git_blob_hash(b),hashlib.sha1(b"blob 16\0"+b).hexdigest())
  def test_read_original_F_not_plant_nitrogen(self):
    h="taxon_name_binom,F_PLANTATT,F_Doring,taxonomic_status,aggregate\n"
    rows="Lotus pedunculatus,8,8,Accepted,\nLotus corniculatus,4,4,Accepted,\n"
    out,n=trait_lookup((h+rows).encode('utf-8'),min_rows=2)
    self.assertEqual(out["Lotus pedunculatus"]["F_PLANTATT"],8)
    self.assertEqual(out["Lotus corniculatus"]["F_PLANTATT"],4)
    self.assertEqual(n,2)
  def test_full_cohort_wet_one_only_no_false_confirmatory(self):
    source=[]
    for i in range(13):
      name=f"Plant {i}"
      source.append({"source_name":name,
        "name_match_method":"original_exact_binomial",
        "F_PLANTATT":8. if i==12 else 4.+(i%3),
        "Kew_taxon_aggregate_flag":"",
        "source_lightcode_log_biomass_ratio":{
          "L2/L1":i/25.,"L3/L1":i/30.,"L3/L2":i/140.}})
    a=infer(source,99)
    self.assertEqual(a["source_F_nonmissing"],13)
    self.assertFalse(a["binary_wet_moisture_test_supported"])
    self.assertEqual(a["moisture_above_7_original_plant_names"],1)
    self.assertEqual(a["permutations"],99)
    self.assertTrue(a["no_ALAN_causal_moisture_moderator_identified"])
  def test_source_synonym_explicit_only(self):
    row=[{"species":"Agrostis tenuis","codes":{
      "L2/L1":0.,"L3/L1":0.,"L3/L2":0.}}]
    info={"Agrostis capillaris":{
      "F_PLANTATT":5.,"F_Doring":"5","Kew_taxonomic_status":"Accepted",
      "Kew_aggregate_label":""}}
    z=describe(row,info)[0]
    self.assertEqual(z["independent_Kew_name"],"Agrostis capillaris")
    self.assertEqual(z["name_match_method"],"authoritative_Kew_POWO_synonym")
if __name__=="__main__":unittest.main()
