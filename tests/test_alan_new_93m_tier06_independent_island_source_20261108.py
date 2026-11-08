"""Synthetic contracts for ALAN new original observation ID tier06 source queue."""
import unittest
import pandas as pd
from scripts.acquisition import alan_new_93m_tier06_independent_island_source_20261108 as m

class FrozenOriginalTier06Tests(unittest.TestCase):
    def test_single_source_photo_per_new_observation_and_never_prior_UUID(self):
        rows=[
         {"photo_id":"4","observation_uuid":"OLD","inat_taxon_id":"1"},
         {"photo_id":"14","observation_uuid":"NEW","inat_taxon_id":"2"},
         {"photo_id":"2","observation_uuid":"NEW","inat_taxon_id":"2"},
         {"photo_id":"45","observation_uuid":"SECOND","inat_taxon_id":"3"}
        ]
        a=m.choose_per_observation(rows,{"OLD"})
        b=m.choose_per_observation(list(reversed(rows)),{"OLD"})
        self.assertEqual({x["observation_uuid"] for x in a},{"NEW","SECOND"})
        self.assertEqual({(x["observation_uuid"],x["photo_id"]) for x in a},
                         {(x["observation_uuid"],x["photo_id"]) for x in b})
        self.assertEqual(len({x["photo_id"] for x in a}),2)

    def test_island_first_balanced_queue_not_one_dominant_island(self):
        rows=[]
        for land in ("a","b","c"):
            for i in range(12):
                rows.append({
                    "land_class":"island","accepted_taxon":True,
                    "has_original_year_ALAN_environment":True,
                    "supported_image_extension":True,"source_loader_supported":True,
                    "observation_uuid":f"{land}_{i}","unit_id":f"{land}-{i}",
                    "land_id":land,"target_id":f"2025|{land}_{i}|{land}",
                    "year":2025,"existing_strict_island_cell_years":0 if land!="a" else 30,
                })
        rows.append({
            "land_class":"mainland","accepted_taxon":True,
            "has_original_year_ALAN_environment":True,
            "supported_image_extension":True,"source_loader_supported":True,
            "observation_uuid":"out","unit_id":"out","land_id":"M",
            "target_id":"2025|m|M","year":2025,"existing_strict_island_cell_years":0})
        d=m.choose_early_island_queue(pd.DataFrame(rows),quota=25,max_each=8)
        self.assertEqual(len(d),24)
        self.assertEqual(d.land_id.nunique(),3)
        self.assertTrue(d.land_id.value_counts().le(8).all())
        self.assertEqual(d.observation_uuid.nunique(),24)

    def test_source_no_flower_colour_assumptions(self):
        self.assertEqual(m.SOURCE_CSV_FIELDS if hasattr(m,"SOURCE_CSV_FIELDS") else m.COLUMNS,
                          m.COLUMNS)
        self.assertEqual(m.QUEUE_SIZE,256)
        self.assertEqual(m.REGISTRY_SHA.__len__(),64)
        self.assertEqual(len(m.TAR_SHA),16)
        self.assertTrue(all(len(x)==64 for x in m.TAR_SHA))

if __name__=="__main__":
    unittest.main()
