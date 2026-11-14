#!/usr/bin/env python3
"""Synthetic XML unit test; CI must fetch and reparse original publisher JATS."""
import sys, unittest, xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_Boom2020_nursery_pollinator_seedpredation_spectra_20261114 import (
  ARTICLE_DOI, SOURCE_NAMES,EXPECTED, AVERAGE_OVULES,get_table,decomposition
)

def synthetic_xml():
    a=ET.Element("article")
    front=ET.SubElement(a,"front")
    t=ET.SubElement(front,"article-id",{"pub-id-type":"doi"})
    t.text=ARTICLE_DOI
    tb=ET.SubElement(a,"table-wrap")
    cap=ET.SubElement(tb,"caption")
    cap.text="Proportions of flowers fertilized, fruits infested and ovules fertilized"
    table=ET.SubElement(tb,"table")
    for k in SOURCE_NAMES:
        p,c,o,e,n1,n2,n3=EXPECTED[k]
        tr=ET.SubElement(table,"tr")
        for v in (k,f"{p} ± .01 ({n1:,})",f"{c} ± .01 ({n2:,})",
                  f"{o} ± .01 ({n3:,})",f"{e} ± 20"):
            td=ET.SubElement(tr,"td")
            td.text=v
    return a

class OriginalNurseryLight(unittest.TestCase):
    def test_exact_publisher_table_values_and_source_nested_counts(self):
        rows=get_table(synthetic_xml())
        self.assertEqual(set(rows),set(SOURCE_NAMES))
        self.assertEqual(rows["Dark"]["fruit_set_source_nested_flowers"],689)
        self.assertEqual(rows["Green"]["infestation_source_nested_fruits"],429)
        self.assertEqual(rows["Green"]["ovule_sample_nested_in_fruit_plant_transect_site"],29127)

    def test_true_reproductive_cost_factors_not_causal(self):
        z=decomposition(get_table(synthetic_xml()))
        self.assertAlmostEqual(z["Green"]["algebraic_source_product_ratio_factors"]["fertilized_flower_factor"],.75/.73)
        self.assertAlmostEqual(z["Green"]["algebraic_source_product_ratio_factors"]["escape_seed_predation_factor"],.47/.72)
        self.assertLess(z["Green"]["ratio_to_dark_recomputed"],.53)
        self.assertGreater(z["Green"]["ratio_to_dark_recomputed"],.50)
        self.assertTrue(z["Green"]["not_causal_mediator_effects"])
        self.assertAlmostEqual(z["Dark"]["recomputed_from_rounded_source_means"],.73*.72*.56*AVERAGE_OVULES)

    def test_fail_wrong_doi_and_missing_arm(self):
        a=synthetic_xml()
        node=next(a.iter("article-id"))
        node.text="10.9999/wrong"
        with self.assertRaisesRegex(ValueError,"wrong DOI"):get_table(a)
        b=synthetic_xml()
        tr=next(t for t in b.iter("tr") if t[0].text=="Red")
        tr.clear()
        with self.assertRaisesRegex(ValueError,"4 light-spectrum arms"):get_table(b)

if __name__=="__main__":unittest.main()
