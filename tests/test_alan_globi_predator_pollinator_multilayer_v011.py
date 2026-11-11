"""SYNTHETIC source fixtures only: do not cite invented animal links as observations."""
from __future__ import annotations
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
import pandas as pd

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'analysis'))
import alan_globi_predator_pollinator_multilayer_v011 as m

TAXA={
 'Gekko gecko':('Gekkonidae','Squamata','Reptilia','Animalia','Animalia|Reptilia|Squamata|Gekkonidae'),
 'Phelsuma ornata':('Gekkonidae','Squamata','Reptilia','Animalia','Animalia|Reptilia|Squamata|Gekkonidae'),
 'Pteropus vampyrus':('Pteropodidae','Chiroptera','Mammalia','Animalia','Animalia|Mammalia|Chiroptera|Pteropodidae'),
 'Scolopendra subspinipes':('Scolopendridae','Scolopendromorpha','Chilopoda','Animalia','Animalia|Arthropoda|Chilopoda'),
 'Araneus diadematus':('Araneidae','Araneae','Arachnida','Animalia','Animalia|Arthropoda|Arachnida|Araneae'),
 'Actias luna':('Saturniidae','Lepidoptera','Insecta','Animalia','Animalia|Arthropoda|Insecta|Lepidoptera'),
 'Noctua pronuba':('Noctuidae','Lepidoptera','Insecta','Animalia','Animalia|Arthropoda|Insecta|Lepidoptera'),
 'Vanessa cardui':('Nymphalidae','Lepidoptera','Insecta','Animalia','Animalia|Arthropoda|Insecta|Lepidoptera'),
 'Ceiba pentandra':('Malvaceae','Malvales','Magnoliopsida','Plantae','Plantae|Embryophyta|Malvales'),
 'Hibiscus rosa':('Malvaceae','Malvales','Magnoliopsida','Plantae','Plantae|Embryophyta|Malvales')}

def row(source,target,relation,study='Field study',life_stage=None):
    v={'interactionTypeName':relation,'studyTitle':study,
       'sourceNamespace':'primary','referenceCitation':'source'}
    for side,species in (('source',source),('target',target)):
        family,order,cl,kingdom,path=TAXA[species]
        v.update({side+'TaxonName':species,side+'TaxonSpeciesName':species,
             side+'TaxonFamilyName':family,side+'TaxonOrderName':order,
             side+'TaxonClassName':cl,side+'TaxonKingdomName':kingdom,
             side+'TaxonPathNames':path})
    if life_stage:v['targetLifeStageName']=life_stage
    return v

def fixture():
    return pd.DataFrame([
      row('Gekko gecko','Actias luna','eats'),
      row('Gekko gecko','Hibiscus rosa','visitsFlowersOf'),
      row('Ceiba pentandra','Phelsuma ornata','pollinatedBy'),
      row('Phelsuma ornata','Noctua pronuba','preysOn',life_stage='adult'),
      row('Pteropus vampyrus','Ceiba pentandra','pollinates'),
      row('Pteropus vampyrus','Noctua pronuba','preysOn',study='iNaturalist archive'),
      row('Scolopendra subspinipes','Pteropus vampyrus','preysOn'),
      row('Scolopendra subspinipes','Actias luna','eats',life_stage='larva'),
      row('Vanessa cardui','Araneus diadematus','preyedUponBy'),
      row('Actias luna','Araneus diadematus','preyedUponBy'),
      row('Gekko gecko','Noctua pronuba','kills'),
    ]).fillna('')

class MultiplexTests(unittest.TestCase):
    def run_fixture(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)
        source=self.path/'synthetic.csv.gz'
        fixture().to_csv(source,index=False,compression='gzip')
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        return m.audit(source,self.path/'results',expected_sha=digest,
                       synthetic=True,chunk_size=3)

    def test_multiguild_same_species_roles(self):
        result=self.run_fixture()
        by={x['predator_group']:x for x in result['guild_summary']
            if x['scope']=='all'}
        self.assertEqual(result['all_original_interaction_rows'],11)
        self.assertEqual(by['geckos']['animal_species_both_moth_predation_and_flower_positive'],2)
        self.assertEqual(by['bats']['animal_species_both_moth_predation_and_flower_positive'],1)
        self.assertEqual(by['centipedes']['animal_species_moth_predation_recorded'],1)
        self.assertEqual(by['spiders']['animal_species_moth_predation_recorded'],1)
        self.assertEqual(result['qualified_role_counts']['moth_kill_only'],1)
        self.assertEqual(result['moth_prey_stage_counts']['adult'],1)
        self.assertEqual(result['moth_prey_stage_counts']['nonadult'],1)
        self.assertTrue(any(x['predator_group']=='centipedes' and x['prey_group']=='bats'
                            for x in result['predator_predator_group_links']))

    def test_inverse_relations_and_literature_filter(self):
        result=self.run_fixture()
        by={x['predator_group']:x for x in result['guild_summary']
            if x['scope']=='non_inaturalist'}
        self.assertEqual(by['geckos']['animal_species_both_moth_predation_and_flower_positive'],2)
        self.assertEqual(by['bats']['animal_species_moth_predation_recorded'],0)
        self.assertEqual(by['bats']['animal_species_pollinates_relation_recorded'],1)
        positive=pd.read_csv(self.path/'results'/'animal_flower_plant_positive_flags.csv.gz')
        self.assertTrue(((positive.animal_group=='bats')&
            (positive.plant_species=='Ceiba pentandra')).any())

    def test_butterfly_prey_is_not_moth(self):
        d=pd.DataFrame([row('Gekko gecko','Vanessa cardui','eats')]).astype('string')
        self.assertFalse(m.moth(d,'target').iloc[0])

    def test_unrecorded_predation_not_demonstrated_absence(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)
            src=path/'synthetic.csv.gz'
            pd.DataFrame([row('Gekko gecko','Ceiba pentandra',
                'visitsFlowersOf')]).to_csv(src,index=False,compression='gzip')
            report=m.audit(src,path/'result',
               expected_sha=hashlib.sha256(src.read_bytes()).hexdigest(),
               synthetic=True)
            self.assertEqual(report['guild_summary'][0]['animal_species_moth_predation_recorded'],0)
            self.assertTrue(report['limitations']['missing_GloBI_link_not_actual_absence'])

    def test_sha_lock_cannot_be_bypassed(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'mock.csv.gz'
            p.write_bytes(b'bad')
            with self.assertRaisesRegex(ValueError,'SHA256 changed'):
                m.audit(p,Path(t)/'output',expected_sha='0'*64,synthetic=True)
            with self.assertRaisesRegex(ValueError,'Refuse unpinned'):
                m.audit(p,Path(t)/'output',expected_sha='0'*64)

if __name__=='__main__':
    unittest.main()
