#!/usr/bin/env python3
"""Numerical and provenance guardrails for original GloPL × direct autonomous selfing."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.analysis.alan_glopl_direct_autonomous_selfing_PCNLPrior_20261113 import (
    fit_source, source_traits, EXPOSURE, TARGET, PRIMARY
)

class OriginalGloPLAutonomousSourceTests(unittest.TestCase):
    def test_refuse_untrusted_direct_plant_trait_file(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"fake.csv"
            p.write_text("accepted_species,trait_name,normalized_value\nFake alba,autonomous_selfing_capacity,autonomous\n")
            with self.assertRaises(ValueError):
                source_traits(p)

    def test_not_estimated_for_pseudoreplicated_small_study(self):
        a=pd.DataFrame({"autonomy":[0,0,1,1],
                        "study_key":["a","a","a","a"],
                        "species":["x","x","y","y"]})
        r=fit_source(a,PRIMARY)
        self.assertEqual(r["status"],"NOT_ESTIMABLE_INDEPENDENT_STUDY_TRAIT_SUPPORT")

    def test_synthetic_original_publication_cluster_interaction_sign(self):
        rng=np.random.default_rng(1013)
        d=[]
        for study in range(64):
            for plant in range(4):
                light=float(rng.uniform(.1,3.5))
                a=int((study+plant)%2)
                d.append({
                    "study_key":f"s{study}", "species":f"sp{(study*4+plant)%104}",
                    "pcnl_prior_mean_log1p":light,
                    "autonomy":a,
                    "PL_log_ratio":.15 +.18*light-.32*light*a +.001*rng.normal(),
                    "exp_year":2004, "lat":float(rng.uniform(-10,10)),
                    "lat2":0.0,"lon":float(rng.uniform(-10,10)),"lon2":0.0,
                    "PL_Effect_Size_Type1":"s", "PL_Effect_Size_Type2":"Sup",
                    "Constant_added":"False","Level_of_Supplementation":"flower",
                    "geo10":"0_0"
                })
        x=pd.DataFrame(d)
        # Both lat2 and lon2 are zero, so the full original-formula design is
        # deliberately rank deficient. The source model MUST report that,
        # never produce a misleading small p value.
        r=fit_source(x,PRIMARY)
        self.assertEqual(r["status"],"NOT_ESTIMABLE_DESIGN_COLLINEARITY")
        # Remove collinear columns, validate interaction calculation.
        f=f"PL_log_ratio ~ {EXPOSURE}*autonomy + lat + lon"
        r=fit_source(x,f)
        self.assertEqual(r["status"],
                         "EXPLORATORY_ORIGINAL_EXPERIMENT_PL_PRIOR_PCNL_TRAIT_MODERATION")
        self.assertAlmostEqual(r["beta"],-.32,delta=.01)
        self.assertLess(r["p_nominal_study_cluster"],.05)
        self.assertIsNotNone(r["two_way_study_and_species_cluster"].get("CI95_student_t"))

if __name__=="__main__":
    unittest.main()
