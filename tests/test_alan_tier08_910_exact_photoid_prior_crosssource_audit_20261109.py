"""Synthetic-only source/photo-site integrity, not original RGB evidence."""
import unittest
import pandas as pd
from scripts.analysis import alan_tier08_910_exact_photoid_prior_crosssource_audit_20261109 as a

class Tier08PhotoDedupTests(unittest.TestCase):
    def test_empty_source_is_not_a_strict_cell(self):
        z=pd.DataFrame(columns=[
           "target_id","photo_id","observation_id","acquisition_arm",
           "year","wcvp_accepted_plant_name_id","strict9"])
        self.assertEqual(a.strict_sites(z)["strict_cell_years"],0)

    def test_two_distinct_sameyear_sites_pass_only_with_three_real_measured_taxa(self):
        f=[]
        for site in ("a","b"):
            for k in range(5):
                f.append(dict(target_id=f"2024|{site}|1",
                   photo_id=f"{site}{k}",observation_id=f"{site}{k}",
                   wcvp_accepted_plant_name_id=f"sp{k}",
                   strict9="true" if k<=2 else "false",
                   year="2024",
                   acquisition_arm="area_stratified_sameyear_two_site_light_pair"))
        x=pd.DataFrame(f)
        self.assertEqual(a.strict_sites(x)["strict_cell_years"],2)
        self.assertEqual(a.strict_sites(x)["strict_pair_island_years"],1)
        t=x.loc[x.photo_id.ne("a0")].copy()
        self.assertEqual(a.strict_sites(t)["strict_cell_years"],1)
        self.assertEqual(a.strict_sites(t)["strict_pair_island_years"],0)

    def test_legacy_source_photo_checksum_frozen(self):
        self.assertEqual(a.N["tier08"],910)
        self.assertEqual(len(a.SHA["tier08"]),64)
        self.assertEqual(a.N["canonical"],1032882)

if __name__=="__main__":
    unittest.main()
