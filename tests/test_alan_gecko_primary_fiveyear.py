"""Ensure coarse coexistence requires the recorded gecko that actually has the prey link."""
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/analysis"))
from alan_gecko_primary_fiveyear import ID, linked_prey


class TestPredatorVisitorLink(unittest.TestCase):
    def test_specific_gecko_and_named_prey(self):
        sites = pd.DataFrame(dict(key=["i","m"], year=[2024,2024], land_id=["1","2"],
                                  land_class=["island","mainland"], x1=[0,0], y1=[0,0]))
        land = sites.assign(species=["Hemidactylus frenatus","Hemidactylus mabouia"])
        members = sites[["key"]].assign(wcvp_accepted_name="Plant example", **{ID: "1"})
        visits = pd.DataFrame(dict(guild=["moth_family_visit"],insect_name=["Moth example"],plant_name=["Plant example"]))
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            prey = folder / "species_prey.csv"
            prey.write_text("gecko,moth,source_stage,non_inat_trophich,relationship_rows\n"
                            "Hemidactylus frenatus,Moth example,unknown,False,2\n"
                            "Hemidactylus mabouia,,adult,False,1\n")
            (folder / "receipt.json").write_text(json.dumps(dict(scanned=24453954,
                source_sha256="670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3",
                prey_table_sha256=hashlib.sha256(prey.read_bytes()).hexdigest())))
            r, _ = linked_prey(sites, members, visits, land, folder)
        rows = {(v["gecko"],v["prey_stage"]):v for v in r}
        self.assertEqual(rows[("Hemidactylus frenatus","any")]["site_years_by_land"], {"island":1})
        self.assertFalse(rows[("Hemidactylus frenatus","adult")]["site_years_by_land"])
        self.assertFalse(rows[("Hemidactylus mabouia","any")]["site_years_by_land"])

    def test_local_chain_needs_same_moth_species_and_stage(self):
        sites = pd.DataFrame(dict(key=["a","b"], year=[2024]*2, land_id=["7"]*2, land_class=["island"]*2,
                                  x1=[0]*2, y1=[0]*2, gx=[3600,3601], gy=[1800]*2))
        geckos = sites.assign(species="Hemidactylus frenatus")
        moths = sites.assign(species=["Feltia subterranea","Spoladea recurvalis"],lifeStage=["Adult","Larva"])
        members = sites[["key"]].assign(wcvp_accepted_name="Plant example", **{ID:"1"})
        visits = pd.DataFrame(dict(guild=["moth_family_visit"]*2,insect_name=moths.species,plant_name=["Plant example"]*2))
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            f = folder / "species_prey.csv"
            f.write_text("gecko,moth,source_stage,non_inat_trophich,relationship_rows\n"
                         "Hemidactylus frenatus,Spoladea recurvalis,adult,False,1\n"
                         "Hemidactylus mabouia,Feltia subterranea,adult,False,1\n")
            (folder/"receipt.json").write_text(json.dumps(dict(scanned=24453954,
                source_sha256="670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3",
                prey_table_sha256=hashlib.sha256(f.read_bytes()).hexdigest())))
            _, chains = linked_prey(sites,members,visits,geckos,folder,moths)
        for c in chains:
            expected = {}
            if c["gecko"]=="Hemidactylus frenatus" and c["gbif_stage"]=="any_stage":
                expected = {"island":2 if c["scale"]=="1deg" else 1}
            self.assertEqual(c["site_years_by_land"],expected)
