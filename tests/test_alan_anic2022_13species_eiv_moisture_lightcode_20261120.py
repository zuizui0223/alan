"""Moisture × original light codes: no fabricated zeros or wet taxa."""
import unittest
from scripts.analysis.alan_anic2022_13species_eiv_moisture_lightcode_source_gate_20261120 import phenotype_modulation

def species(i,trait):
    return {"species":f"Plant_{i}",
       "MOIST_exact":trait,"wet7":trait>7 if trait is not None else None,
       "logratio_light_codes":{"L2/L1":.1*(trait or 0),
                                "L3/L1":-.05*(trait or 0),
                                "L3/L2":-.15*(trait or 0)}}

class MOISTGates(unittest.TestCase):
    def test_sparse_wet_species_stop_even_if_rho_seems_strong(self):
        x=[species(i,4.+i/10) for i in range(10)]
        x += [species(10,8.),species(11,8.5),species(12,9.)]
        y=phenotype_modulation(x)
        self.assertEqual(y["status"],"SOURCE_MOISTURE_MODERATOR_COMMON_SUPPORT_STOP")
        self.assertIsNone(y["multiple_contrast_confirmatory_inference"])

    def test_unknown_moisture_never_called_dry(self):
        x=[species(i,float(5+i/2)) for i in range(10)]
        x+=[species(10,None),species(11,None),species(12,None)]
        y=phenotype_modulation(x)
        self.assertEqual(y["matched_original_plants"],10)
        self.assertNotEqual(y.get("matched_original_plants"),13)

    def test_three_source_light_codes_are_not_relabelled_light_intensity(self):
        x=[species(i,2.5+i*.7) for i in range(13)]
        y=phenotype_modulation(x)
        self.assertIn("descriptive_by_original_light_code",y)
        self.assertTrue(all(k in y["descriptive_by_original_light_code"]
                            for k in ("L2/L1","L3/L1","L3/L2")))
        if y["status"].startswith("EXPLORATORY"):
            self.assertTrue(all(0<=p<=1 for p in y["max3_absolute_rho_adjusted_permutation_p"].values()))

if __name__=="__main__":unittest.main()
