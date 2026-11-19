import io
import unittest
import zipfile

from scripts.analysis.alan_sanders2020_stageAB_source_census_20261119 import (
    original_source_rows,source_panel_census,REQUIRED
)
from scripts.acquisition.alan_sanders2020_source_126paper_Dryad_preflight_20261119 import parse_zip


def row(study,phylum,clazz,order,description,sub,meanL,meanD):
    d={k:"" for k in REQUIRED}
    d.update({
       "Study":study,"ID":study+"_1","Phylum":phylum,
       "Class":clazz,"Order":order,"Species":"source taxon",
       "Measure.description":description,"Category ":"populations",
       "Sub.Category":sub,"Field.Lab":"field","Control_N":"4",
       "Experimental_N":"4","Control_Mean":str(meanD),
       "Experimental_Mean":str(meanL),"Control_SD":"1",
       "Experimental_SD":"1","Night_Lux":"10",
       "Location":"test site","DOI":"10.example/test","Delete":"",
    })
    return d


class Original126ALANMetaStageScreenTests(unittest.TestCase):
    def test_synthetic_stage_taxonomy_separates_pollen_and_predator_without_meta_claim(self):
        rows=[
            row("StudyA","Arthropoda","Arachnida","Araneae","spider abundance","abundance",5,3),
            row("StudyB","Tracheophyta","Magnoliopsida","","seed output","reproductive output",3,1),
            row("StudyC","Arthropoda","Insecta","Lepidoptera","flower visitation","activity",4,2),
            row("StudyD","Chordata","Mammalia","","reproductive output","reproduction",1,1)
        ]
        out=source_panel_census(rows,{"original_CSV_SHA256":"test","original_R_code_SHA256":"test"},strict=False)
        b=out["stageAB_broad_record_screening_not_eligible_effects"]
        self.assertEqual(b["S1_arthropod_predation_abundance_candidate"]["source_unique_original_Study_ids_preliminary"],1)
        self.assertEqual(b["S2_plant_reproductive_endpoint_candidate"]["source_unique_original_Study_ids_preliminary"],1)
        self.assertEqual(b["S2_invertebrate_potential_flower_visitation_candidate"]["source_unique_original_Study_ids_preliminary"],1)
        self.assertTrue(out["source_original_effect_n_values_not_confirmed_independent_lamp_sites"])
        self.assertTrue(out["requires_publication_fulltext_taxonomy_and_habitat_before_stageA"])
        self.assertTrue(out["no_wet_plant_moisture_moderator_fitted"])
        self.assertFalse(out["source_raw_study_IDs_or_rows_exported"])
        self.assertNotIn("StudyA",str(out))
        self.assertNotIn("StudyB",str(out))

    def test_original_angiosperms_and_monocots_are_not_falsely_excluded(self):
        rows=[
            row("PlantStudyA","angiosperms","monocots","","seed output","Reproductive.output",5,3),
            row("PlantStudyB","anthophyta","dicotyledons","","fruit set","Reproductive.output",6,4),
            row("PlantStudyC","tracheophyta","liliopsida","","flower phenology","Plant.phy",8,7),
        ]
        out=source_panel_census(rows,{"original_CSV_SHA256":"synthetic","original_R_code_SHA256":"synthetic"},strict=False)
        b=out["stageAB_broad_record_screening_not_eligible_effects"]
        self.assertEqual(b["S2_plant_reproductive_endpoint_candidate"]["source_unique_original_Study_ids_preliminary"],2)
        self.assertEqual(b["S2_plant_reproductive_endpoint_candidate"]["source_positive_arms_SD_n_arithmetically_available_records"],2)
        self.assertEqual(b["S2_plant_flowering_phenology_candidate_NOT_seed_outcome"]["source_unique_original_Study_ids_preliminary"],1)
        self.assertEqual(b["S2_plant_flowering_phenology_candidate_NOT_seed_outcome"]["source_positive_arms_SD_n_arithmetically_available_records"],1)

    def test_nonpositive_means_not_lnRR_ready_even_with_large_flower_sample(self):
        rows=[row("StudyA","Tracheophyta","Magnoliopsida","","seed output","reproduction",0,3)]
        d=source_panel_census(rows,{"original_CSV_SHA256":"x","original_R_code_SHA256":"y"},strict=False)
        x=d["stageAB_broad_record_screening_not_eligible_effects"]["S2_plant_reproductive_endpoint_candidate"]
        self.assertEqual(x["source_original_effect_records_preliminary_keyword_and_taxonomy"],1)
        self.assertEqual(x["source_positive_arms_SD_n_arithmetically_available_records"],0)

    def test_sanders_real_fields_required_for_authentic_source_extraction(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:
            z.writestr("Sanders_et_al_datafile.csv","Study,Mean\nStudyA,0.1\n")
            z.writestr("Sanders_meta__rcode.r","# study")
        with self.assertRaisesRegex(ValueError,"not source verified"):
            original_source_rows(b.getvalue(),strict=False)


if __name__=="__main__":
    unittest.main()
