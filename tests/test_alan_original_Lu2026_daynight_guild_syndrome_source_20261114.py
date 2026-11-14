#!/usr/bin/env python3
"""Field block replication and source visit conservation mathematical guards."""
import sys,unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_original_Lu2026_daynight_guild_syndrome_source_20261114 import (
 exact_paired_block_sign,original_design_model,species_descriptive,XLSX_SHA
)

class SourceVerifiedFiveBlock(unittest.TestCase):
 def test_sha_pinned(self):
  self.assertEqual(len(XLSX_SHA),64)
  self.assertEqual(XLSX_SHA,"ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f")
 def test_signflip_support_not_pseudoreplication(self):
  x=exact_paired_block_sign([-2,-2,-2,-2,-2])
  self.assertEqual(x["two_sided_complete_32_signflip_p_DESCRIPTIVE"],.0625)
  self.assertTrue(x["five_block_values_all_same_sign"])
 def test_conservation_paired_original_light_design(self):
  a=[]
  for time in range(1,5):
   for p in ("Diurnal","Nocturnal"):
    for t in ("Lit","Unlit"):
     for b in range(1,6):
      for rep in range(1,6):
       n=1 if t=="Unlit" else 2
       a.append({"time":time,"sampling_time":p,"treatment":t,
          "block":b,"repeat":rep,"abundance":5*n,
          "D_abundance":2*n,"H_abundance":n,"L_abundance":n,
          "C_abundance":n,"O_abundance":0})
  q=original_design_model(pd.DataFrame(a))
  self.assertEqual(q["Nocturnal"]["independent_field_blocks"],5)
  self.assertEqual(q["Nocturnal"]["experimental_survey_pair_group_effects"]["Total"]["n_paired_block_rounds"],20)
  self.assertAlmostEqual(q["Nocturnal"]["experimental_survey_pair_group_effects"]["Total"]["lit_dark_visit_ratio"],2)
 def test_species_source_unknown_without_plant_flowering_denominator(self):
  z=pd.DataFrame([{"plant":"X","time":1,"sampling_time":"Nocturnal",
     "treatment":"Lit","block":1,"repeat":1,"all":2,
     "Lepidoptera":1,"Coleoptera":0,"Diptera":1,"Hymenoptera":0,"Others":0}] * 1038)
  a=species_descriptive(z)
  self.assertTrue(a["source"]["no_missing_plant_survey_row_coded_zero"])
  self.assertEqual(a["by_period"]["Diurnal"]["status"],"NOT_ESTIMABLE_NO_POSITIVE_SOURCE_COMMON_SUPPORT")

if __name__=="__main__":
 unittest.main()
