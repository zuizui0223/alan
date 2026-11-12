#!/usr/bin/env python3
"""Strict species-direct morphology classification tests; no synthetic ecology claims."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from alan_direct_flower_syndrome_light_20261112 import (
    classify, taxon_support, stratum_test, site_count, regression,
)

class TestSourceTraits(unittest.TestCase):
    def test_single_species_direct_unknown_not_negative(self):
        species=pd.DataFrame({
            "wcvp_accepted_plant_name_id":["1","2","3","4","5","6"],
            "wcvp_accepted_name":["A a","B b","C c","D d","E e","F f"],
            "class":["old_only","old_only","common","common","new_only","neither"]
        })
        traits=pd.DataFrame({
            "flower_primary_color":["white","white","blue_purple",pd.NA,"white","white"],
            "floral_form":["tubular","funnel_trumpet","tubular","salverform","open_radial",pd.NA],
            "tube_depth_class":["deep",pd.NA,"deep","deep","deep","deep"]
        },index=["A a","B b","C c","D d","E e","F f"])
        result=classify(species,traits)
        self.assertEqual(result.broad_candidate.tolist(),[True,True,False,False,False,False])
        self.assertEqual(result.strict_candidate.tolist(),[True,False,False,False,False,False])
        self.assertEqual(result.eligible_broad.tolist(),[True,True,True,False,True,False])
        self.assertEqual(result.broad_other_known.tolist(),[False,False,True,False,True,False])
        self.assertEqual(result.eligible_strict.tolist(),[True,False,True,False,True,False])
        tally=taxon_support(result)
        self.assertEqual(tally["old_only"]["broad_candidate"],2)
        self.assertEqual(tally["old_only"]["strict_candidate"],1)

    def test_posthoc_fisher_source_trait_presence(self):
        x=pd.DataFrame({"class":["old_only"]*10+["common"]*10,
                        "eligible_broad":[True]*20,
                        "broad_candidate":[True]*4+[False]*6+[True]*2+[False]*8})
        out=stratum_test(x,"broad_candidate",gate="eligible_broad")
        self.assertEqual((out["old_only_positive"],out["common_positive"]),(4,2))
        self.assertEqual((out["old_only_denominator"],out["common_denominator"]),(10,10))
        self.assertGreater(out["odds_ratio"],1)

    def test_photo_site_counts_count_each_taxon_once(self):
        sites=pd.DataFrame({"key":["k1","k2","k3"]})
        obs=pd.DataFrame({"key":["k1","k1","k2"],"white":[True,False,True]})
        a=site_count(sites,obs,obs.white)
        self.assertEqual(a.tolist(),[1,1,0])
        b=regression(pd.DataFrame({"key":["k1","k2","k3"],
            "n_species":[2,2,2],"land_class":["island","island","mainland"],
            "region5":["a","b","c"]}),a,"synthetic_support_only")
        self.assertEqual(b["model_status"],"insufficient_independent_geographic_or_photo_support")
        self.assertIsNone(b["beta"])

if __name__=="__main__":
    unittest.main()
