"""Source provenance and non-iNaturalist claim separation: synthetic code tests."""
from __future__ import annotations
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
from scripts.analysis import alan_globi_inaturalist_source_independence_20261110 as p


class TestIndependence(unittest.TestCase):
    def test_one_exact_source_accepted_species_is_not_duplicated(self):
        with tempfile.TemporaryDirectory() as tmp:
            with sqlite3.connect(Path(tmp)/"allflora-island-group-00.sqlite") as conn:
                conn.execute("""CREATE TABLE flora (
                    island TEXT,latbin INTEGER,lonbin INTEGER,year INTEGER,
                    species TEXT,moth_positive_photo_ids INTEGER,
                    butterfly_positive_photo_ids INTEGER)""")
                conn.executemany("INSERT INTO flora VALUES (?,?,?,?,?,?,?)", [
                    ("I",12,34,2023,"Plant alpha",2,1),
                    ("I",12,34,2023,"Plant beta",0,0)])
            s=pd.DataFrame([dict(island="I",latbin=12,lonbin=34,year=2023)])
            g=pd.DataFrame([
                dict(species="Plant alpha",moth_global=True,butterfly_global=True,
                     moth_external=True,butterfly_external=False,
                     any_external=True,both_external=False,
                     moth_inat_only=False,butterfly_inat_only=True,
                     external_moth_no_external_butterfly=True)])
            with patch.object(p,"SOURCE_GROUPS",1):
                d=p.original_species(tmp,s,g)
            self.assertEqual(len(d),2)
            self.assertEqual(int(d.moth_global.sum()),1)
            self.assertEqual(int(d.moth_external.sum()),1)
            self.assertEqual(int(d.moth_inat_only.sum()),0)
            self.assertEqual(int(d.butterfly_inat_only.sum()),1)

    def test_native_and_missing_guild_never_imputed(self):
        self.assertEqual(p.FIELDS[0],"moth_external")
        self.assertNotIn("native",p.FIELDS)
        self.assertNotIn("cultivated",p.FIELDS)


if __name__=="__main__": unittest.main()
