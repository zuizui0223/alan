"""Test only physical block and taxonomic complement logic; synthetic fixture."""
import csv,importlib.util,io,json,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/"scripts/analysis/alan_sullivan2019_other_family_decomposition_20261124.py"
spec=importlib.util.spec_from_file_location("sullivan_taxa_decomp",FILE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def create_source():
    # Test fixture is not biological observation data.
    chiro={
      "DRK":[120,118,118],"LEDI":[400,400,362],
      "LEDH":[500,400,395],"HPSI":[500,500,457],"HPSH":[200,230,260]
    }
    other={
      "DRK":[17,17,17],"LEDI":[20,15,15],
      "LEDH":[25,20,22],"HPSI":[36,36,36],"HPSH":[45,45,46]
    }
    cols=["Station","Treatment","Group","Sample","Date","Taxa","Count"]
    d=[]
    for i in range(3):
        for tr in ("DRK",)+m.LIGHTS:
            for tax,n in (("Chironomidae",chiro[tr][i]),("Culicidae",other[tr][i])):
                d.append({"Station":f"{tr}-{i+1}","Treatment":tr,
                    "Group":str(i+1),"Sample":"AquaticEmergent",
                    "Date":"10/01/2013","Taxa":tax,"Count":str(n)})
    while len(d)<198:
        d.append({"Station":"DRK-1","Treatment":"DRK","Group":"1",
            "Sample":"AquaticEmergent","Date":"10/01/2013",
            "Taxa":"Culicidae","Count":"0"})
    for _ in range(963):
        d.append({"Station":"DRK-1","Treatment":"DRK","Group":"1",
            "Sample":"LeafFall","Date":"10/01/2013",
            "Taxa":"leaf","Count":""})
    assert len(d)==1161
    fh=io.StringIO()
    writer=csv.DictWriter(fh,fieldnames=cols,delimiter="\t")
    writer.writeheader();writer.writerows(d)
    return fh.getvalue().encode()

class TaxonDecompositionTest(unittest.TestCase):
    def test_block_effort_not_independent(self):
        a=m.inspect(create_source())
        self.assertEqual(a["source_total_non_chironomidae_catches"],412)
        self.assertEqual(a["true_wetland_habitat_blocks"],3)
        self.assertEqual(len(a["physical_block_family_contributions"]),3)
        self.assertEqual(a["source_family_label_richness"],2)
        self.assertEqual(a["source_non_chironomidae_family_counts"]["Culicidae"],412)
        self.assertEqual(a["new_meta_programmes_added"],0)
        alt=a["explicit_two_aquatic_larval_family_sensitivity"]
        self.assertEqual(alt["physical_habitat_block_n"],3)
        self.assertTrue(alt["source_selection_posthoc"])
        self.assertEqual(alt["source_selected_families"],["Chironomidae","Culicidae"])
        self.assertEqual(len(alt["three_original_block_counts_and_ratios"]),3)
        self.assertFalse(not alt["not_independent_experimental_meta_effect"])
        self.assertTrue(alt["descriptive_t2_95_RRR"][0]<alt["equal_block_RRR"]<alt["descriptive_t2_95_RRR"][1])
        self.assertGreaterEqual(alt["three_block_signflip_reference_p"],0.25)
        self.assertTrue(all(x["chiro_vs_all_other_response_RRR"]>0 for x in a["physical_block_family_contributions"]))

    def test_source_frozen_receipt_rejects_synthetic_fixture(self):
        pinned=json.loads((ROOT/"database/alan_experiments/arthropod_expansion_20261122/chironomid_emergence_20261124/sullivan_family_composition_frozen_20261124.json").read_text(encoding="utf-8"))
        self.assertEqual(pinned["original_physical_habitat_groups"],3)
        self.assertAlmostEqual(pinned["equal_block_chiro_vs_other_RRR"],1.59110636626693,places=12)
        self.assertAlmostEqual(pinned["equal_block_chiro_vs_culicidae_RRR"],2.705677495636492,places=12)
        self.assertEqual(pinned["new_independent_programmes"],0)
        with self.assertRaises(ValueError):
            m.verify_frozen(m.inspect(create_source()),pinned)

    def test_corrupted_taxon_source_blocks_stop(self):
        raw=create_source().decode()
        with self.assertRaises(ValueError):
            m.inspect(raw.replace("Chironomidae","DipteraX",1).encode())
        with self.assertRaises(ValueError):
            m.inspect(raw.replace("5372","5370",1).encode() if "5372" in raw else raw.replace("\t120\r\n","\t1200\r\n",1).encode())
        with self.assertRaises(ValueError):
            m.inspect(raw.replace("AquaticEmergent","Not_Aquatic",1).encode())

if __name__=="__main__":
    unittest.main()
