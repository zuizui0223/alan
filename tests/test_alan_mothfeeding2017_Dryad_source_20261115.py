import importlib.util
import unittest
from pathlib import Path
SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/alan_mothfeeding2017_Dryad_source_20261115.py"
spec=importlib.util.spec_from_file_location("mothfeeding",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestMothFeeding(unittest.TestCase):
    def test_detect_source_long_schema(self):
        data=b"Day\tBlock\tTreatment\tSpecies\tSex\tunique.id\tObs\tSize\n2020-01-01\t1\tDark\tMothA\tM\tA\t1\t1\n2020-01-01\t1\tDark\tMothA\tM\tA\t0\t1\n"
        rows,metadata=m.parse_original(data)
        self.assertEqual(len(rows),2)
        self.assertEqual(metadata["source_column_names"],["Day","Block","Treatment","Species","Sex","unique.id","Obs","Size"])
    def test_individual_not_time_is_n(self):
        rows=[]
        for id,arm,values in [("a","Dark",[1,0,1,1]),("b","Green",[0,1,0,0])]:
            for obs in values:
                rows.append({"Day":"1","Block":"1","Treatment":arm,"Species":"A","Sex":"F","unique.id":id,"Obs":str(obs),"Size":"1"})
        d=m.derive(rows)
        self.assertEqual(d["total_original_longitudinal_observations"],8)
        self.assertEqual(d["unique_individual_moths_with_binary_data"],2)
        self.assertEqual(len(d["results"]),2)
        self.assertEqual(next(x for x in d["results"] if x["lamp_treatment"]=="Green")["matched_block_mean_difference_vs_dark"],-.5)
    def test_fail_closed_if_reused_id_crosses_lights(self):
        rows=[
          {"Day":"1","Block":"1","Treatment":"Dark","Species":"A","Sex":"F","unique.id":"id","Obs":"1","Size":"1"},
          {"Day":"1","Block":"1","Treatment":"Green","Species":"A","Sex":"F","unique.id":"id","Obs":"0","Size":"1"}]
        with self.assertRaises(ValueError):m.derive(rows)
    def test_reject_unknown_dark(self):
        rows=[{"Day":"1","Block":"1","Treatment":"Green","Species":"A","Sex":"F","unique.id":"id","Obs":"1","Size":"1"}]
        with self.assertRaises(ValueError):m.derive(rows)

if __name__=="__main__":
    unittest.main()
