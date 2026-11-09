"""Synthetic validation only; never interpreted as measured 93m source totals."""
import unittest
from collections import Counter
from scripts.analysis import alan_full93m_globi_moth_original_photo_census_20261109 as m

class Global93mMothSourceUnitTests(unittest.TestCase):
    @staticmethod
    def row(photo,uuid,taxon,tier=0,good=True,year=2024):
        d={key:"" for key in m.COLUMNS}
        d.update({"photo_id":str(photo),"observation_uuid":str(uuid),
                  "inat_taxon_id":str(taxon),"tier":str(tier),
                  "valid_observed_date":"true" if good else "false",
                  "valid_geopoint":"true","frozen_loader_supported":"true",
                  "observed_on":str(year)+"-06-01"})
        return d

    def test_moth_visitor_pollen_walton_and_butterfly_not_exclusive(self):
        # 1= moth visit, 2=GloBI pollinates, 4=Walton nocturnal pollen,
        # 8= butterfly-positive. An obs with two separate original photos
        # counts 2 photos but only one biological observation UUID.
        d=[self.row(1,"u1",100),self.row(2,"u1",100),
           self.row(3,"u2",200),self.row(4,"u3",300,good=False),
           self.row(5,"u4",400)]
        stats=Counter()
        m.classify_original_photo_rows(d,{"100":1|2|4|8,"200":1|8,"300":8},stats,{},0)
        self.assertEqual(stats["all_source_photos"],5)
        self.assertEqual(stats["original_photo_IDs_moth_visit_positive"],3)
        self.assertEqual(stats["original_photo_IDs_moth_database_pollinates_positive"],2)
        self.assertEqual(stats["original_photo_IDs_walton2020_nocturnal_pollen_transport_positive"],2)
        self.assertEqual(stats["original_photo_IDs_butterfly_visit_positive"],4)
        self.assertEqual(stats["source_RGB_eligible_original_photo_IDs_moth_visit_positive"],3)
        self.assertEqual(stats["independent_source_observation_UUIDs_moth_visit_positive"],2)
        self.assertEqual(stats["original_photos_unknown_WCVP_accepted_name"],1)

    def test_unrecorded_GloBI_stays_unknown_not_negative(self):
        d=[self.row(1,"u1",100),self.row(2,"u2",200)]
        stat=Counter()
        m.classify_original_photo_rows(d,{"100":0},stat,{},0)
        self.assertEqual(stat["accepted_taxon_photos_without_positive_record"],1)
        self.assertEqual(stat["original_photos_unknown_WCVP_accepted_name"],1)
        self.assertEqual(stat["original_photo_IDs_moth_visit_positive"],0)

    def test_no_malformed_or_cross_tier_original_photo_source(self):
        s=Counter()
        with self.assertRaisesRegex(ValueError,"identity invalid"):
            m.classify_original_photo_rows([self.row(1,"u",10,tier=9)],
                                           {"10":1},s,{},0)

    def test_conflicting_one_observation_taxon_raises(self):
        s=Counter()
        with self.assertRaisesRegex(ValueError,"contradictory botanical taxa"):
            m.classify_original_photo_rows([
                self.row(1,"sameUUID",1),self.row(2,"sameUUID",2)],
                {"1":1,"2":8},s,{},0)
if __name__=="__main__":
    unittest.main()
