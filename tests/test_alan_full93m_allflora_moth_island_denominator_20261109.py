"""Synthetic tests for one-denominator GloBI-positive island plant census."""
from __future__ import annotations
import sqlite3
import unittest
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

from scripts.analysis import alan_full93m_allflora_moth_island_denominator_20261109 as mod
from scripts.analysis import alan_full93m_globi_moth_original_photo_census_20261109 as source


def fake(photo, uuid, plant, x=0.25, y=0.25, yr="2024"):
    vals = {k: "" for k in source.COLUMNS}
    vals.update(photo_id=photo, observation_uuid=uuid, inat_taxon_id=plant,
                observed_on=yr+"-06-05", latitude=str(y),longitude=str(x),
                valid_observed_date="true",valid_geopoint="true",
                frozen_loader_supported="true",tier="0")
    return vals


class TestAllFloraDenominator(unittest.TestCase):
    def setUp(self):
        self.islands=gpd.GeoDataFrame({"island":["Test-A"]},
            geometry=[box(-1,-1,1,1)],crs="EPSG:4326")
        self.taxon={"10":"Moth-reported Plant","20":"Unreported Plant"}
        self.evidence={"10":1|8,"20":0}

    def dataframe(self,*rows):
        return pd.DataFrame(rows,columns=list(source.COLUMNS))

    def test_moth_positive_and_unrecorded_both_remain_in_denominator(self):
        r=[fake("1","ob1","10"),fake("2","ob1","10"),
           fake("3","ob2","20"),fake("4","ob3","30")]
        normalized=mod.normalize_rows(self.dataframe(*r),self.taxon,self.evidence)
        self.assertEqual(len(normalized),3)
        tuples,n,amb=mod.grouped_island_rows(normalized,self.islands,set())
        self.assertEqual(n,3)
        self.assertEqual(amb,0)
        self.assertEqual(len(tuples),2)
        w={x[4]:x for x in tuples}
        self.assertEqual(w["Moth-reported Plant"][5:10],(2,1,2,1,0))
        self.assertEqual(w["Unreported Plant"][5:10],(1,1,0,0,0))

    def test_unknown_guild_not_biological_negative(self):
        x=mod.normalize_rows(self.dataframe(fake("1","ob1","20")),self.taxon,self.evidence)
        self.assertEqual(len(x),1)
        self.assertEqual(int(x.iloc[0].moth_mask),0)

    def test_original_observation_uuid_not_double_counted_across_chunks(self):
        seen=set()
        a=mod.normalize_rows(self.dataframe(fake("1","same","10")),self.taxon,self.evidence)
        b=mod.normalize_rows(self.dataframe(fake("2","same","10")),self.taxon,self.evidence)
        r1,_,_=mod.grouped_island_rows(a,self.islands,seen)
        r2,_,_=mod.grouped_island_rows(b,self.islands,seen)
        self.assertEqual(r1[0][6],1)
        self.assertEqual(r2[0][6],0)
        self.assertEqual(r2[0][5],1)

    def test_original_photo_year_and_landgate(self):
        x=mod.normalize_rows(self.dataframe(fake("1","a","10",yr="2022"),
                   fake("2","b","10",x=100),fake("3","c","10")),self.taxon,self.evidence)
        self.assertEqual(len(x),2)
        _,n,_=mod.grouped_island_rows(x,self.islands,set())
        self.assertEqual(n,1)

    def test_sql_aggregation_preserves_species_grain(self):
        db=sqlite3.connect(":memory:")
        db.execute(mod.TABLE)
        base=("Island",2,3,2024,"A")
        db.execute(mod.UPSERT,(*base,2,1,2,1,0,0,2,0.2,0.2))
        db.execute(mod.UPSERT,(*base,1,1,1,1,0,0,1,0.2,0.2))
        db.execute(mod.UPSERT,("Island",2,3,2024,"B",2,2,0,0,0,0,0,0.3,0.3))
        r=db.execute("SELECT COUNT(DISTINCT species),SUM(photo_ids),SUM(observation_uuids),"
                     "SUM(CASE WHEN moth_positive_photo_ids>0 THEN 1 ELSE 0 END) FROM flora").fetchone()
        self.assertEqual(r,(2,5,4,1))

    def test_bad_photo_id_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"identity malformed"):
            mod.normalize_rows(self.dataframe(fake("invalid","ob","10")),self.taxon,self.evidence)


if __name__=="__main__":
    unittest.main()
