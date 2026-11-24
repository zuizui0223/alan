"""Checks source-unit and taxon-gate logic without mocking original data as real."""
import csv,importlib.util,io
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/"scripts/analysis/alan_sullivan2019_chironomid_source_gate_20261124.py"
sp=importlib.util.spec_from_file_location("sullivan_taxon",FILE)
mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)

def fixture():
    # Artificial *test fixture* mimicking 198 aquatic records and 963 other
    # original input records: never used to calculate empirical results.
    fields=["Station","Treatment","Group","Sample","Date","Taxa","Count"]
    rows=[]
    for group in ("1","2","3"):
        for treatment in mod.ARMS:
            for day in ("d1","d2","d3","d4"):
                rows.append(dict(Station=treatment+"-"+group,Treatment=treatment,
                    Group=group,Sample="AquaticEmergent",Date=day,Taxa="Chironomidae",Count="0"))
    for j in range(138):
        rows.append(dict(Station="DRK-1",Treatment="DRK",Group="1",
                         Sample="AquaticEmergent",Date="d1",
                         Taxa="Other aquatic family",Count="0"))
    rows[0]["Count"]="5372"
    for j in range(963):
        rows.append(dict(Station="DRK-1",Treatment="DRK",Group="1",
                         Sample="LeafFall",Date="d1",Taxa="leaf",Count=""))
    assert len(rows)==1161
    out=io.StringIO(newline="")
    w=csv.DictWriter(out,fieldnames=fields,delimiter="\t")
    w.writeheader();w.writerows(rows)
    return out.getvalue().encode("utf-8")

class TestSullivanTaxonSource(unittest.TestCase):
    def test_source_contract_identity(self):
        self.assertEqual(mod.DIGEST,"385444ab318ec4b15edd3161c740aa7f4b2d9b973fe4c2ddf6c29cc56cdec484")
        self.assertEqual(mod.URL,"https://datadryad.org/api/v2/versions/8787/download")
        self.assertEqual(len(mod.ARMS),5)

    def test_three_source_groups_not_fifteen_independent_programmes(self):
        data=mod.analyze(fixture())
        self.assertEqual(data["independent_habitat_blocks"],3)
        self.assertEqual(data["original_stations"],15)
        self.assertEqual(data["total_aquatic_catches"],5372)
        self.assertEqual(data["chironomid_source_labels"],["Chironomidae"])
        self.assertEqual(data["chironomid_count_by_treatment"]["DRK"],5372)
        self.assertEqual(data["meta_k_added"],0)

    def test_source_schema_and_event_corruption_stops(self):
        raw=fixture().decode("utf-8")
        with self.assertRaises(ValueError):
            mod.analyze(raw.replace("AquaticEmergent","BrokenLabel",1).encode("utf-8"))
        with self.assertRaises(ValueError):
            mod.analyze(raw.replace("5372","-1",1).encode("utf-8"))
        with self.assertRaises(ValueError):
            mod.analyze(raw.replace("DRK-1","DRK-4",1).encode("utf-8"))

if __name__=="__main__":
    unittest.main()
