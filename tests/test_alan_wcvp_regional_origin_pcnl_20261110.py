"""Synthetic correctness gates for source-verified WCVP botanical-country status.

Artificial tests do NOT establish actual island native plant ecology.
"""
from __future__ import annotations

import unittest
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import box

from scripts.analysis import alan_wcvp_regional_origin_pcnl_20261110 as m


class RegionalWCVPContractTest(unittest.TestCase):
    def test_unique_tdwg_region_and_boundary_unknown(self):
        geo = gpd.GeoDataFrame(
            {"tdwg":["AAA","BBB"]},
            geometry=[box(0,0,1,1), box(1,0,2,1)], crs="EPSG:4326")
        sites=pd.DataFrame([
            dict(island="1",latbin=20,lonbin=40,year=2025,lat=.5,lon=.5),
            dict(island="1",latbin=20,lonbin=41,year=2025,lat=.5,lon=1.5),
            dict(island="1",latbin=20,lonbin=42,year=2025,lat=.5,lon=1.0),
        ])
        mapped=m.overlay_sites(sites,geo)
        self.assertEqual(mapped.tdwg_level3_code.tolist(),["AAA","BBB",""])
        self.assertEqual(mapped.tdwg_region_state.iloc[-1],"outside_or_boundary")

    def test_assessed_native_introduced_unknown_conserve_species_and_moth(self):
        rows=[]
        for j,name,yes in [(0,"A species",1),(1,"B species",0),(2,"C species",1),
                           (3,"D species",1)]:
            rows.append(dict(island="1",latbin=10,lonbin=20,year=2025,
                species=name,moth_positive_photo_ids=bool(yes),
                butterfly_positive_photo_ids=False,
                moth_pollinates_photo_ids=False,night_pollen_photo_ids=False))
        species=pd.DataFrame(rows)
        sites=pd.DataFrame([dict(island="1",latbin=10,lonbin=20,year=2025,
            all_accepted_recorded_species=4,
            moth_visitor_positive_recorded_species=3,
            all_accepted_observation_uuids=6,
            lon=.5,lat=.5)])
        place=pd.DataFrame([dict(island="1",latbin=10,lonbin=20,year=2025,
            tdwg_level3_code="AAA",tdwg_region_state="unique_pinned_TDWG3")])
        native=pd.DataFrame([dict(species="A species",tdwg_level3_code="AAA",
            native_introduced_state="native_at_TDWG3"),
            dict(species="B species",tdwg_level3_code="AAA",
            native_introduced_state="introduced_at_TDWG3"),
            dict(species="C species",tdwg_level3_code="AAA",
            native_introduced_state="doubtful")])
        out=m.summarise_status(species,place,native,sites)
        q=out.iloc[0]
        self.assertEqual(q.recorded_species__native_at_TDWG3,1)
        self.assertEqual(q.recorded_species__introduced_at_TDWG3,1)
        self.assertEqual(q.recorded_species__unknown,2)
        self.assertEqual(q.recorded_moth__native_at_TDWG3,1)
        self.assertEqual(q.recorded_moth__introduced_at_TDWG3,0)
        self.assertEqual(q.recorded_moth__unknown,2)

    def test_species_restricted_island_contrast_support_does_not_fake_beta(self):
        rows=[]
        for island in range(4):
            for year in (2019,2025):
                for site in (1,2):
                    rows.append({"island":str(island),"year":year,"latbin":site,
                        "lonbin":site,"pcnl_1992_to_previous_year_mean_log1p":0.1*site,
                        "all_accepted_observation_uuids":8,
                        "recorded_species__native_at_TDWG3":5,
                        "recorded_moth__native_at_TDWG3":2})
        result=m.model_group(pd.DataFrame(rows),"native_at_TDWG3")
        self.assertEqual(result["status"],"not_estimable")
        self.assertEqual(result["eligible_islands"],4)

    def test_botanical_country_native_not_island_or_individual_wild(self):
        self.assertEqual(m.STATUS,("native_at_TDWG3","introduced_at_TDWG3"))
        self.assertNotIn("wild",m.STATUS)


if __name__=="__main__":unittest.main()
