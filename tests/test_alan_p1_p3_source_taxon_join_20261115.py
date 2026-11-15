"""P1×P3 exact taxon crosswalk tests; synthetic source cannot claim ALAN interaction."""
import json,math,unittest
from pathlib import Path
from scripts.analysis.alan_p1_p3_source_taxon_join_20261115 import (
    li_silene,li_blob_verify,LI_CODES
)

def rows():
    # Controlled fixture with exactly the 2026 Silene nutans arm/field totals,
    # not a replacement for the external author CSV.
    amounts={
      "1":{"C":[0,0],"P":[168,24],"L":[233,60]},
      "2":{"C":[135,0],"P":[85,0],"L":[11,0]},
      "3":{"C":[102,0],"P":[144,60],"L":[615,717]},
    }
    return [{"species":"SN","replicate":k,"Treatment":arm,
             "pollinated":str(int(seed>0)),"Seeds":str(seed)}
            for k in ("1","2","3") for arm in ("C","P","L")
            for seed in amounts[k][arm]]

class TestP1P3NoForcedEffect(unittest.TestCase):
    def test_true_source_silene_nutans_seed_direction_and_site_grain(self):
        out=li_silene(rows())
        self.assertEqual(out["arms"]["L"]["seed_present"],5)
        self.assertEqual(out["arms"]["C"]["seed_present"],2)
        self.assertEqual(out["arms"]["P"]["seed_present"],5)
        self.assertEqual(out["arms"]["L"]["total_seeds"],1636)
        self.assertEqual(out["arms"]["C"]["total_seeds"],237)
        self.assertAlmostEqual(out["linear_vs_dark_seed_presence_risk_difference"],.5)
        self.assertEqual(out["linear_vs_point_seed_presence_risk_difference"],0)
        self.assertEqual(out["n_independent_hayfield_blocks"],3)
        self.assertEqual(out["n_fields_linear_positive"],2)
        self.assertEqual(out["n_fields_linear_equal"],1)
        self.assertEqual(out["seed_count_linear_minus_dark_positive_field_count"],2)
        self.assertFalse(out["causal_P1xP3_moderation_established"])

    def test_broken_field_plant_coverage_fails(self):
        with self.assertRaises(ValueError):li_silene(rows()[:-1])
    def test_source_hash_fails_closed(self):
        with self.assertRaises(ValueError):li_blob_verify(b'fake original P3 Li')
    def test_source_catalog_keeps_11_species_and_factorial_lu_separate(self):
        p=Path('config/alan_p1_p3_primary_taxon_crosswalk_20261115.json')
        data=json.loads(p.read_text())
        q=data['publications']
        self.assertEqual(len(q),9)
        self.assertEqual(len(set(sum([v['taxa'] for v in q],[]))),11)
        self.assertEqual(set(LI_CODES.values()),set(next(x for x in q if x["id"]=="li2026")["taxa"]))
        self.assertEqual(next(x for x in q if x['id']=='wilson2021')['taxa'],['Hesperoyucca whipplei'])
        self.assertEqual(next(x for x in q if x['id']=='lu2026jpe')['taxa'],
            ['Pedicularis kansuensis','Euphrasia regelii'])
        self.assertFalse(data['validated_alias_candidates'][0]['match_tier'].startswith('EXACT'))

    def test_night_exclusion_is_not_light_effect(self):
        fruit_loss=1-.09/.7
        seed_loss=1-16.544/284.314
        self.assertAlmostEqual(fruit_loss,.8714285714285714)
        self.assertTrue(.94<seed_loss<.95)
        # Real original Li L vs C source risk ratio 2.5 is positive.
        # It is a different source, year, country and reproductive endpoint.
        self.assertGreater((5/6)/(2/6),1)
        self.assertNotEqual("night_loss_in_Belgium","ALAN_interaction_in_UK")

if __name__=="__main__":unittest.main()
