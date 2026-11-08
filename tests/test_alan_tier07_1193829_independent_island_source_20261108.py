"""Synthetic only: signed tier07 metadata are never assumed to be floral pixels."""
import unittest
from scripts.acquisition import alan_tier07_1193829_independent_island_source_20261108 as s

class ALANTier07OriginalSourceTests(unittest.TestCase):
    def test_all_16_original_tier07_source_groups_have_explicit_SHA256(self):
        self.assertEqual(len(s.TAR_SHA),16)
        self.assertTrue(all(len(x)==64 and all(c in "0123456789abcdef" for c in x) for x in s.TAR_SHA))
        self.assertEqual(len(s.REGISTRY_SHA),64)
        self.assertEqual(s.QUEUE_SIZE,256)
        self.assertEqual(s.YEARS,(2019,2021,2023,2024,2025))

    def test_one_photo_per_new_UUID_colour_blind_and_excludes_prior(self):
        rows=[
          {"photo_id":"2","observation_uuid":"new","inat_taxon_id":"5"},
          {"photo_id":"8","observation_uuid":"prior","inat_taxon_id":"5"},
          {"photo_id":"3","observation_uuid":"new","inat_taxon_id":"5"},
          {"photo_id":"15","observation_uuid":"also","inat_taxon_id":"6"},
        ]
        q=s.choose_per_observation(rows,{"prior"})
        t=s.choose_per_observation(list(reversed(rows)),{"prior"})
        self.assertEqual({x["observation_uuid"] for x in q},{"new","also"})
        self.assertEqual({(x["observation_uuid"],x["photo_id"]) for x in q},
                         {(x["observation_uuid"],x["photo_id"]) for x in t})
        self.assertEqual(len(q),2)

if __name__=="__main__":
    unittest.main()
