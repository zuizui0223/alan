"""Tests only parser/role logic: synthetic claims do NOT validate ecology."""
import hashlib,sys,tempfile,unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'/'analysis'))
import alan_globi_v011_predator_flower_multirole as m

BAT=('Chiroptera','Mammalia','Vespertilionidae')
GECKO=('Squamata','Reptilia','Gekkonidae')
CENTIPEDE=('Scutigeromorpha','Chilopoda','Scutigeridae')
SPIDER=('Araneae','Arachnida','Araneidae')
def row(animal,group,partner,relation,kind='plant',reverse=False,source='field'):
 order,cl,fam=group
 a={'Name':animal,'SpeciesName':animal,'OrderName':order,'ClassName':cl,
    'FamilyName':fam,'KingdomName':'Animalia','PathNames':'Animalia|'+cl+'|'+order+'|'+fam}
 if kind=='moth':p={'Name':partner,'SpeciesName':partner,'OrderName':'Lepidoptera',
    'ClassName':'Insecta','FamilyName':'Noctuidae','KingdomName':'Animalia',
    'PathNames':'Animalia|Insecta|Lepidoptera|Noctuidae'}
 elif kind=='predator':
  p={'Name':partner,'SpeciesName':partner,'OrderName':'Squamata',
    'ClassName':'Reptilia','FamilyName':'Gekkonidae','KingdomName':'Animalia',
    'PathNames':'Animalia|Reptilia|Squamata|Gekkonidae'}
 else:p={'Name':partner,'SpeciesName':partner,'OrderName':'Myrtales',
    'ClassName':'Magnoliopsida','FamilyName':'Myrtaceae','KingdomName':'Plantae',
    'PathNames':'Plantae|Tracheophyta|Magnoliopsida|Myrtales'}
 x,y=(p,a) if reverse else (a,p)
 return {'interactionTypeName':relation,'sourceNamespace':source,
   **{'sourceTaxon'+k:v for k,v in x.items()},
   **{'targetTaxon'+k:v for k,v in y.items()}}
class Tests(unittest.TestCase):
 def fixture(self):
  return pd.DataFrame([
    row('Pipistrellus kuhlii',BAT,'Agrotis ipsilon','preysOn',kind='moth'),
    row('Pipistrellus kuhlii',BAT,'Silene latifolia','visitsFlowersOf'),
    row('Phelsuma ornata',GECKO,'Trochetia blackburniana','pollinates'),
    row('Scutigera coleoptrata',CENTIPEDE,'Agrotis ipsilon','eatenBy',kind='moth',reverse=True),
    row('Scutigera coleoptrata',CENTIPEDE,'Phelsuma ornata','preysOn',kind='predator'),
    row('Scutigera coleoptrata',CENTIPEDE,'Silene latifolia','eats',source='iNaturalist'),
    row('Araneus diadematus',SPIDER,'Agrotis ipsilon','preysOn',kind='moth')])
 def test_predators_plant_visitors_and_intraguild(self):
  q=m.parse(self.fixture())
  self.assertEqual(len(q),7)
  self.assertEqual(set(q[q.animal=='Pipistrellus kuhlii'].role),{'predates_moth','visits_flowers'})
  self.assertEqual(set(q[q.animal=='Scutigera coleoptrata'].role),
     {'predates_moth','preys_on_predator_group','eats_plant_not_pollination'})
  self.assertEqual(q[q.animal=='Phelsuma ornata'].role.tolist(),['db_pollinates'])
 def test_no_flower_feeding_as_pollination(self):
  q=m.parse(self.fixture())
  self.assertFalse(q[(q.animal=='Scutigera coleoptrata')].role.isin(
       ['visits_flowers','db_pollinates']).any())
 def test_unresolved_genus_not_independent_species(self):
  q=self.fixture()
  q.loc[0,'sourceTaxonName']='Pipistrellus'
  q.loc[0,'sourceTaxonSpeciesName']=''
  self.assertEqual(len(m.parse(q)),6)
 def test_all_source_and_exclusion_tiers(self):
  with tempfile.TemporaryDirectory() as d:
   source=Path(d)/'source.csv.gz';dest=Path(d)/'out'
   self.fixture().to_csv(source,index=False,compression='gzip')
   digest=hashlib.sha256(source.read_bytes()).hexdigest()
   r=m.scan(source,dest,expected_sha=digest,fixture=True,chunk_size=4)
   b=r['grouped_positive_species_counts']['bat']
   self.assertEqual(b['same_exact_animal_name_both_predation_and_flower_use'],1)
   c=r['grouped_positive_species_counts']['centipede']
   self.assertEqual(c['db_claimed_pollination'],0)
   self.assertEqual(c['predates_other_predator_group'],1)
   self.assertEqual(c['plant_eating_not_pollination'],1)
   self.assertEqual(r['all_interaction_rows'],7)
 def test_wrong_hash_and_overwrite_fail_closed(self):
  with tempfile.TemporaryDirectory() as d:
   source=Path(d)/'s.csv.gz';self.fixture().to_csv(source,index=False,compression='gzip')
   with self.assertRaisesRegex(ValueError,'Unpinned'):m.scan(source,Path(d)/'out',expected_sha='0'*64)
   with self.assertRaisesRegex(ValueError,'SHA256'):
    m.scan(source,Path(d)/'out',expected_sha='0'*64,fixture=True)
if __name__=='__main__':unittest.main()
