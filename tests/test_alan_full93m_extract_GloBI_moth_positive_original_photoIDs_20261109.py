"""Tests use SYNTHETIC source photo IDs, never substitute for source 93m output."""
import csv
import io
import unittest
from collections import Counter
from scripts.acquisition import alan_full93m_extract_GloBI_moth_positive_original_photoIDs_20261109 as m

class OriginalFullMothPhotoIDExtractionTests(unittest.TestCase):
    def row(self,pid,uuid,taxon,tier=7,valid=True):
        d={k:"" for k in m.SOURCE_TIER_FIELDS}
        d.update(photo_id=str(pid),observation_uuid=str(uuid),
                 inat_taxon_id=str(taxon),observed_on="2024-09-21",
                 valid_observed_date="true" if valid else "false",
                 valid_geopoint="true",frozen_loader_supported="true",
                 photo_license="cc-by",image_extension="jpeg",
                 longitude="139.9",latitude="35.5",tier=str(tier))
        return d

    def test_proven_guild_source_only_preserves_independent_UUID_and_photo_ID(self):
        src=[self.row(1,"u1",100),self.row(2,"u1",100),
             self.row(3,"u2",200),self.row(4,"u3",300)]
        out=io.StringIO()
        writer=csv.DictWriter(out,fieldnames=m.OUTPUT_FIELDS)
        writer.writeheader()
        stat=Counter();obs=set();taxa=Counter()
        seen=m.extract_rows(src,{"100":1|2|4|8,"200":1|8,"300":0},
                            writer,7,stat,obs,taxa)
        self.assertEqual(seen,4)
        self.assertEqual(stat["moth_positive_original_photo_IDs"],3)
        self.assertEqual(stat["moth_database_pollinates_photo_IDs"],2)
        self.assertEqual(stat["walton_nocturnal_pollen_carriage_photo_IDs"],2)
        self.assertEqual(stat["also_butterfly_positive_photo_IDs"],3)
        self.assertEqual(stat["true_source_RGB_eligible_moth_photo_IDs"],3)
        self.assertEqual(len(obs),2)
        emitted=list(csv.DictReader(io.StringIO(out.getvalue())))
        self.assertEqual([x["photo_id"] for x in emitted],["1","2","3"])
        self.assertEqual(emitted[0]["moth_database_pollinates"],"true")
        self.assertEqual(emitted[2]["moth_database_pollinates"],"false")

    def test_technical_unknown_is_not_promoted_to_actual_RGB(self):
        src=[self.row(1,"u1",10,valid=False)]
        stat=Counter()
        writer=csv.DictWriter(io.StringIO(),fieldnames=m.OUTPUT_FIELDS)
        m.extract_rows(src,{"10":1},writer,7,stat,set(),Counter())
        self.assertEqual(stat["moth_positive_original_photo_IDs"],1)
        self.assertEqual(stat["true_source_RGB_eligible_moth_photo_IDs"],0)

    def test_invalid_tier_rejected_before_emitting_photographs(self):
        writer=csv.DictWriter(io.StringIO(),fieldnames=m.OUTPUT_FIELDS)
        with self.assertRaisesRegex(ValueError,"immutable photo"):
            m.extract_rows([self.row(1,"u1",10,tier=8)],
                           {"10":1},writer,7,Counter(),set(),Counter())

    def test_source_evidence_flags_are_positive_not_daytime_control(self):
        stat=Counter()
        writer=csv.DictWriter(io.StringIO(),fieldnames=m.OUTPUT_FIELDS)
        m.extract_rows([self.row(5,"u5",10),self.row(6,"u6",11)],
                       {"10":8},writer,7,stat,set(),Counter())
        self.assertEqual(stat["moth_positive_original_photo_IDs"],0)
if __name__=="__main__":
    unittest.main()
