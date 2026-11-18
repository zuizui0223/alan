"""Original plant×pollinator data: never invent missing plants' zero visits."""
import unittest
from scripts.analysis.alan_lu2026_original_sameplant_guild_replacement import (
    indexed_original,match_plants,source_results
)

def r(phase,arm,plant,block,round_,rep,moths,flies):
    return {"phase":phase,"arm":arm,"plant":plant,"block":block,
        "round":round_,"repeat":rep,
        "all":moths+flies,"Lepidoptera":moths,"Diptera":flies,
        "Coleoptera":0,"Hymenoptera":0,"Others":0}

class TestLuSameFlowerLightBlocks(unittest.TestCase):
    def test_matched_plant_does_not_assume_missing_light_arm_zero(self):
        q=[r("night","lit","A","1","t1","r1",1,3),
           r("night","dark","A","1","t1","r1",4,1),
           r("night","dark","B","1","t1","r1",3,0)]
        d=match_plants(q,["1"])
        for mode,plant in d.items():
            self.assertEqual(len(plant),1,mode)
            self.assertEqual(plant[0]["plant"],"A")
            self.assertEqual(plant[0]["moth_delta"],-3)
            self.assertEqual(plant[0]["fly_delta"],2)

    def test_round_repeat_strict_excludes_unmatched_sampling_period(self):
        q=[r("night","lit","A","1","t2","r1",1,2),
           r("night","dark","A","1","t1","r1",3,0)]
        d=match_plants(q,["1"])
        self.assertEqual(len(d["same_plant_block_any_round"]),1)
        self.assertEqual(len(d["same_plant_block_round_repeat"]),0)

    def test_source_pair_is_not_new_experimental_light_installation(self):
        q=[r("night","lit","A","1","t1","r1",1,3),
           r("night","dark","A","1","t1","r1",4,1),
           r("night","lit","A","2","t1","r1",1,2),
           r("night","dark","A","2","t1","r1",2,0)]
        out=source_results(q,["1","2"])["same_plant_block_round_repeat"]
        self.assertEqual(out["night_same_plant_moth_down_fly_up"]["source_plant_block_pairs"],2)
        self.assertEqual(out["night_source_species_with_at_least_two_independent_block_pairs"],1)
        self.assertTrue(out["all_plant_rows_not_independent_experiments"])

if __name__=="__main__":
    unittest.main()
