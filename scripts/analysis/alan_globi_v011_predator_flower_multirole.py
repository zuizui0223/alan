#!/usr/bin/env python3
"""Full GloBI v0.11: moth predation and animal flower-use POSITIVE evidence.

All 24.45m interactions scanned with official immutable SHA. Same binomial
animal name may have multiple roles, but claim provenance may be from different
locations. A missing record is NEVER predation/pollination absence. Vegetal
feeding NEVER implies pollination; database 'pollinates' is not field efficacy.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd

RELEASE_SHA='670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3'
PRED_FWD={'preysOn','eats','kills','hasPrey'}
PRED_REV={'preyedUponBy','eatenBy','killedBy'}
FLOWER_FWD={'visitsFlowersOf','pollinates'}
FLOWER_REV={'flowersVisitedBy','pollinatedBy'}
BUTTERFLY={'Nymphalidae','Papilionidae','Pieridae','Lycaenidae','Riodinidae','Hesperiidae','Hedylidae'}
GECKOS={'Gekkonidae','Phyllodactylidae','Sphaerodactylidae','Diplodactylidae','Eublepharidae','Carphodactylidae','Pygopodidae'}
WASPS={'Vespidae','Sphecidae','Crabronidae','Pompilidae','Scoliidae'}
BINOMIAL=re.compile(r'^[A-Z][A-Za-z-]+ [a-z][a-z-]+$')
FIELDS={'interactionTypeName','studyTitle','sourceCitation','sourceNamespace','referenceCitation'}
for side in ('source','target'):
 for suffix in ('Name','SpeciesName','OrderName','ClassName','FamilyName','KingdomName','PathNames'):
  FIELDS.add(side+'Taxon'+suffix)

def require(v,msg):
 if not v:raise ValueError(msg)
def file_sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4<<20),b''):h.update(b)
 return h.hexdigest()
def col(d,side,rank):
 key=side+'Taxon'+rank
 return d[key].fillna('').astype('string').str.strip() if key in d else pd.Series('',index=d.index,dtype='string')
def binomial(d,side):
 x=col(d,side,'SpeciesName')
 x=x.where(x.ne(''),col(d,side,'Name'))
 return x.where(x.str.fullmatch(BINOMIAL).fillna(False),'')
def taxon(d,side,rank,value):
 v=col(d,side,rank+'Name')
 p=col(d,side,'PathNames')
 return v.eq(value)|p.str.contains(fr'(?:^|\|){re.escape(value)}(?:\||$)',regex=True)
def group(d,side):
 fam=col(d,side,'FamilyName')
 cl=col(d,side,'ClassName')
 o=pd.Series('',index=d.index,dtype='string')
 rules=[
  ('gecko',fam.isin(GECKOS)),
  ('bat',taxon(d,side,'Order','Chiroptera')),
  ('centipede',taxon(d,side,'Class','Chilopoda')),
  ('spider',taxon(d,side,'Order','Araneae')),
  ('mantis',taxon(d,side,'Order','Mantodea')),
  ('wasp',fam.isin(WASPS)),
  ('frog',taxon(d,side,'Order','Anura')),
  ('bird',taxon(d,side,'Class','Aves'))]
 for label,flag in rules:o=o.mask(o.eq('')&flag,label)
 return o
def plant(d,side):
 return taxon(d,side,'Kingdom','Plantae')|taxon(d,side,'Kingdom','Viridiplantae')|taxon(d,side,'Kingdom','Embryophyta')
def parse(chunk):
 rel=chunk.interactionTypeName.fillna('').astype('string')
 forward=rel.isin(PRED_FWD|FLOWER_FWD)
 reverse=rel.isin(PRED_REV|FLOWER_REV)
 d=chunk.loc[forward|reverse].copy()
 if d.empty:return pd.DataFrame(columns=('animal','guild','role','source_excluded'))
 rev=reverse.loc[d.index]
 animal=binomial(d,'source').where(~rev,binomial(d,'target'))
 category=group(d,'source').where(~rev,group(d,'target'))
 other_fam=col(d,'target','FamilyName').where(~rev,col(d,'source','FamilyName'))
 moth=taxon(d,'target','Order','Lepidoptera').where(~rev,taxon(d,'source','Order','Lepidoptera'))
 plant_target=plant(d,'target').where(~rev,plant(d,'source'))
 other_predator_group=group(d,'target').where(~rev,group(d,'source'))
 relation=rel.loc[d.index]
 ispred=relation.isin(PRED_FWD|PRED_REV)
 isflower=relation.isin(FLOWER_FWD|FLOWER_REV)
 role=pd.Series('',index=d.index,dtype='string')
 role=role.mask(ispred&moth&other_fam.ne('')&~other_fam.isin(BUTTERFLY),'predates_moth')
 role=role.mask(ispred&moth&other_fam.isin(BUTTERFLY),'predates_butterfly')
 role=role.mask(ispred&moth&other_fam.eq(''),'predates_Lepidoptera_unspecified')
 role=role.mask(ispred&plant_target,'eats_plant_not_pollination')
 role=role.mask(ispred&other_predator_group.ne('')&~moth&~plant_target,'preys_on_predator_group')
 role=role.mask(isflower&plant_target&relation.isin({'visitsFlowersOf','flowersVisitedBy'}),'visits_flowers')
 role=role.mask(isflower&plant_target&relation.isin({'pollinates','pollinatedBy'}),'db_pollinates')
 references=pd.Series('',index=d.index,dtype='string')
 for name in ('sourceNamespace','sourceCitation','studyTitle','referenceCitation'):
  if name in d:references=references+' '+d[name].fillna('').astype('string')
 exclusion=references.str.contains('inaturalist|trophich',case=False,regex=True)
 return pd.DataFrame({'animal':animal,'guild':category,'role':role,
   'source_excluded':exclusion}).loc[animal.ne('')&category.ne('')&role.ne('')]
def scan(path,out,*,expected_sha=RELEASE_SHA,fixture=False,chunk_size=160000):
 path,out=Path(path),Path(out)
 require(not out.exists(),'Refuse to overwrite')
 require(expected_sha==RELEASE_SHA or fixture,'Unpinned source cannot claim stable release')
 digest=file_sha(path);require(digest==expected_sha,'Source SHA256 mismatch')
 cols=set(pd.read_csv(path,compression='gzip',nrows=0).columns)
 require({'interactionTypeName','sourceTaxonName','targetTaxonName'}<=cols,'GloBI source columns changed')
 positive=defaultdict(set);record_counts=Counter();i=0
 for data in pd.read_csv(path,compression='gzip',usecols=lambda x:x in (FIELDS&cols),
                           dtype='string',chunksize=chunk_size):
  i+=len(data)
  q=parse(data)
  record_counts.update(q.role.to_list())
  for name,g,role,excluded in q.itertuples(index=False,name=None):
   positive[(g,name)].add(role)
   if not excluded:positive[(g,name)].add(role+'_non_inat_trophich')
  if i and i%4800000<chunk_size:print('GLOBI_V011_FULL_PREDATOR_SCAN',i,flush=True)
 groups=['gecko','bat','centipede','spider','mantis','wasp','frog','bird']
 roles=['predates_moth','predates_butterfly','predates_Lepidoptera_unspecified',
        'visits_flowers','db_pollinates','eats_plant_not_pollination','preys_on_predator_group']
 out.mkdir(parents=True)
 with gzip.open(out/'positive_species_predator_flower_roles.csv.gz','wt',newline='',encoding='utf-8') as f:
  writer=csv.writer(f);writer.writerow(['guild','animal_name']+roles+
       [r+'_non_inat_trophich' for r in roles])
  for (g,n),flags in sorted(positive.items()):
   writer.writerow([g,n]+[int(r in flags) for r in roles]+
             [int(r+'_non_inat_trophich' in flags) for r in roles])
 stats={}
 for g in groups:
  taxa=[x for (label,_),x in positive.items() if label==g]
  pred={'predates_moth','predates_Lepidoptera_unspecified'}
  visitor={'visits_flowers','db_pollinates'}
  stats[g]={'source_positive_taxon_names':len(taxa),
    'documented_predation_on_moths_or_unresolved_Lepidoptera':sum(bool(x&pred) for x in taxa),
    'documented_predation_on_moths_family_confirmed':sum('predates_moth' in x for x in taxa),
    'flower_visitation_or_db_pollination':sum(bool(x&visitor) for x in taxa),
    'same_exact_animal_name_both_predation_and_flower_use':sum(bool(x&pred) and bool(x&visitor) for x in taxa),
    'db_claimed_pollination':sum('db_pollinates' in x for x in taxa),
    'plant_eating_not_pollination':sum('eats_plant_not_pollination' in x for x in taxa),
    'predates_other_predator_group':sum('preys_on_predator_group' in x for x in taxa),
    'moth_predation_non_iNaturalist_trophiCH':sum(bool(x&{r+'_non_inat_trophich' for r in pred}) for x in taxa),
    'flower_use_non_iNaturalist_trophiCH':sum(bool(x&{r+'_non_inat_trophich' for r in visitor}) for x in taxa)}
 report={'schema':'alan_full_GloBI_v011_predator_prey_pollination_multirole_v1',
   'status':'synthetic_fixture' if fixture else 'stable_source_sha_verified',
   'source_sha256':digest,'all_interaction_rows':i,
   'candidate_record_counts':dict(record_counts),'grouped_positive_species_counts':stats,
   'source_positive_only_no_absence_inferred':True,
   'same_species_name_dual_role_not_same_site_or_behavior_event':True,
   'predation_moth_larva_vs_adult_UNRESOLVED':True,
   'db_pollinates_not_direct_verified_pollen_deposition':True,
   'flower_eating_never_upgraded_to_pollination':True,
   'known_gecko_families':sorted(GECKOS),
   'full_source_url':'https://zenodo.org/records/22691479'}
 (out/'summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 return report
if __name__=='__main__':
 a=argparse.ArgumentParser(description=__doc__)
 a.add_argument('--source',required=True,type=Path)
 a.add_argument('--out',required=True,type=Path)
 x=a.parse_args()
 result=scan(x.source,x.out)
 print(json.dumps({'all_records':result['all_interaction_rows'],
                   'by_guild':result['grouped_positive_species_counts']},sort_keys=True))
