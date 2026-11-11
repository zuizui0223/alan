#!/usr/bin/env python3
"""SHA-pinned FULL GloBI v0.11: moth predators, flower visitors and dual roles.

Species-level dual means the same ANIMAL SPECIES has both positive source
records somewhere, not that predation and effective pollination co-occurred.
Caterpillar predation is not automatically predation on flower-visiting adults.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd

RELEASE='https://zenodo.org/records/22691479'
SOURCE_SHA='670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3'
SOURCE_ROWS=24453954
SPECIES=re.compile(r'^[A-Z][A-Za-z-]+ [a-z][a-z-]+$')
BUTTERFLIES={'Papilionidae','Nymphalidae','Pieridae','Lycaenidae','Riodinidae','Hesperiidae','Hedylidae'}
GECKOS={'Gekkonidae','Eublepharidae','Sphaerodactylidae','Phyllodactylidae',
        'Diplodactylidae','Carphodactylidae','Pygopodidae'}
GROUPS=('geckos','bats','centipedes','spiders','birds','amphibians',
        'mantises','dragonflies','other_reptiles','other_mammals',
        'other_insects','other_arthropods')
FOCAL={'geckos','bats','centipedes','spiders'}
SCOPES=('all','non_inaturalist','non_inaturalist_trophich')
BITS={'moth_predation':1,'flower_visit':2,'reported_pollinates':4,'plant_dispersal':8,'moth_kill_only':16}
REL={'eats':('predation',False),'preysOn':('predation',False),
     'eatenBy':('predation',True),'preyedUponBy':('predation',True),
     'kills':('kill',False),'killedBy':('kill',True),
     'visitsFlowersOf':('flower_visit',False),'flowersVisitedBy':('flower_visit',True),
     'pollinates':('reported_pollinates',False),'pollinatedBy':('reported_pollinates',True),
     'dispersalVectorOf':('plant_dispersal',False),
     'hasDispersalVector':('plant_dispersal',True)}
FIELDS={'interactionTypeName','sourceTaxonName','targetTaxonName',
        'sourceTaxonSpeciesName','targetTaxonSpeciesName',
        'sourceTaxonFamilyName','targetTaxonFamilyName',
        'sourceTaxonOrderName','targetTaxonOrderName',
        'sourceTaxonClassName','targetTaxonClassName',
        'sourceTaxonKingdomName','targetTaxonKingdomName',
        'sourceTaxonPathNames','targetTaxonPathNames',
        'sourceLifeStage','targetLifeStage','sourceLifeStageName','targetLifeStageName',
        'sourceNamespace','studyTitle','sourceCitation','referenceCitation'}

def require(ok,msg):
    if not bool(ok):raise ValueError(msg)

def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as file:
        for b in iter(lambda:file.read(1<<22),b''):digest.update(b)
    return digest.hexdigest()

def col(d,prefix,suffix):
    field=prefix+'Taxon'+suffix
    if field not in d:return pd.Series('',index=d.index,dtype='string')
    return d[field].fillna('').str.strip()

def species_name(d,prefix):
    sp=col(d,prefix,'SpeciesName')
    v=sp.where(sp.ne(''),col(d,prefix,'Name'))
    return v.where(v.str.fullmatch(SPECIES).fillna(False),'')

def has(path,token):
    return path.str.contains(r'(?:^|\|)'+re.escape(token)+r'(?:\||$)',regex=True,na=False)

def category(d,prefix):
    family=col(d,prefix,'FamilyName')
    order=col(d,prefix,'OrderName')
    cls=col(d,prefix,'ClassName')
    path=col(d,prefix,'PathNames')
    output=pd.Series('unclassified',index=d.index,dtype='string')
    masks=[
      ('geckos',family.isin(GECKOS)|has(path,'Gekkota')),
      ('bats',order.eq('Chiroptera')|has(path,'Chiroptera')),
      ('centipedes',cls.eq('Chilopoda')|has(path,'Chilopoda')),
      ('spiders',order.eq('Araneae')|has(path,'Araneae')),
      ('birds',cls.eq('Aves')|has(path,'Aves')),
      ('amphibians',cls.eq('Amphibia')|has(path,'Amphibia')),
      ('mantises',order.eq('Mantodea')|has(path,'Mantodea')),
      ('dragonflies',order.eq('Odonata')|has(path,'Odonata')),
      ('other_reptiles',cls.eq('Reptilia')|has(path,'Reptilia')|
                        order.eq('Squamata')|has(path,'Squamata')),
      ('other_mammals',cls.eq('Mammalia')|has(path,'Mammalia')),
      ('other_insects',cls.eq('Insecta')|has(path,'Insecta')),
      ('other_arthropods',has(path,'Arthropoda'))]
    for name,flag in masks:
        output.loc[flag & output.eq('unclassified')]=name
    return output

def moth(d,prefix):
    family=col(d,prefix,'FamilyName')
    order=col(d,prefix,'OrderName')
    path=col(d,prefix,'PathNames')
    return ((order.eq('Lepidoptera')|has(path,'Lepidoptera')) &
            family.str.fullmatch('[A-Z][a-z]+idae').fillna(False) &
            ~family.isin(BUTTERFLIES))

def plant(d,prefix):
    k=col(d,prefix,'KingdomName').str.casefold()
    p=col(d,prefix,'PathNames').str.casefold()
    return (k.isin(('plantae','viridiplantae'))|
       p.str.contains(r'(?:^|\|)(?:plantae|viridiplantae|embryophyta)(?:\||$)',
                      regex=True,na=False))

def stage(d,prefix):
    v=pd.Series('',index=d.index,dtype='string')
    for field in (prefix+'LifeStageName',prefix+'LifeStage'):
        if field in d:v=v.where(v.ne(''),d[field].fillna('').str.lower())
    out=pd.Series('unknown',index=d.index,dtype='string')
    out.loc[v.str.contains('larva|caterpillar|pupa|nymph|egg',regex=True,na=False)]='nonadult'
    out.loc[v.str.contains('adult|imago',regex=True,na=False)]='adult'
    return out

def source_flags(d):
    text=pd.Series('',index=d.index,dtype='string')
    for field in ('sourceNamespace','studyTitle','sourceCitation','referenceCitation'):
        if field in d:text=text+' '+d[field].fillna('')
    text=text.str.casefold()
    ina=text.str.contains('inaturalist',regex=False,na=False)
    tro=text.str.contains('trophich',regex=False,na=False)
    return {'all':pd.Series(True,index=d.index),
            'non_inaturalist':~ina,
            'non_inaturalist_trophich':~ina & ~tro}

def annotate(d):
    """Return (group, actor species, role, partner species, partner group,
    moth life stage, trusted source scopes) for POSITIVE evidence only."""
    d=d.loc[d.interactionTypeName.isin(REL)].copy()
    if d.empty:return []
    flags=source_flags(d)
    inverted=d.interactionTypeName.map(lambda x:REL[x][1])
    result=[]
    for reverse in (False,True):
        q=d.loc[inverted.eq(reverse)].copy()
        if q.empty:continue
        actor='target' if reverse else 'source'
        other='source' if reverse else 'target'
        g=category(q,actor);a=species_name(q,actor)
        valid=g.isin(GROUPS)&a.ne('')
        if not valid.any():continue
        q=q.loc[valid];g=g.loc[valid];a=a.loc[valid]
        kind=q.interactionTypeName.map(lambda x:REL[x][0])
        is_moth=moth(q,other);is_plant=plant(q,other)
        other_name=species_name(q,other);other_group=category(q,other)
        other_stage=stage(q,other)
        accepted=((kind.isin(('predation','kill')) & is_moth) |
                  (kind.isin(('flower_visit','reported_pollinates','plant_dispersal')) &
                   is_plant & other_name.ne('')) |
                  (kind.eq('predation') & other_group.isin(GROUPS) &
                   other_name.ne('')))
        q=q.loc[accepted]
        if q.empty:continue
        tier={name:flag.reindex(q.index,fill_value=False) for name,flag in flags.items()}
        for i in q.index:
            role=kind.loc[i]
            if role in ('predation','kill'):
                if is_moth.loc[i]:
                    role='moth_predation' if role=='predation' else 'moth_kill_only'
                else:role='predator_of_predator'
            result.append((str(g.loc[i]),str(a.loc[i]),role,
                str(other_name.loc[i]),str(other_group.loc[i]),
                str(other_stage.loc[i]),
                tuple(name for name,flags in tier.items() if bool(flags.loc[i]))))
    return result

def audit(source,out,*,expected_sha=SOURCE_SHA,synthetic=False,chunk_size=120000):
    source,out=Path(source),Path(out)
    require(not out.exists(),'Refuse overwrite')
    require(synthetic or expected_sha==SOURCE_SHA,'Refuse unpinned source release')
    digest=sha256(source)
    require(digest==expected_sha,'Full GloBI SHA256 changed')
    cols=set(pd.read_csv(source,compression='gzip',nrows=0).columns)
    require({'interactionTypeName','sourceTaxonName','targetTaxonName'}<=cols,'Schema drift')
    for side in ('source','target'):
        require(side+'TaxonPathNames' in cols or side+'TaxonClassName' in cols,
                'Missing taxonomic hierarchy')
    counts=Counter()
    actor_state={scope:defaultdict(int) for scope in SCOPES}
    plant_state={scope:defaultdict(int) for scope in SCOPES}
    animal_edge_rows=Counter()
    animal_edge_focal=defaultdict(set)
    input_fields=FIELDS & cols
    for i,block in enumerate(pd.read_csv(source,compression='gzip',
        dtype='string',usecols=lambda x:x in input_fields,
        chunksize=chunk_size,low_memory=False)):
        counts['all_interaction_rows']+=len(block)
        relations=block.interactionTypeName.fillna('')
        counts.update({'source_'+key:int(relations.eq(key).sum())
                       for key in REL if relations.eq(key).any()})
        for g,animal,role,partner,other_group,moth_stage,scopes in annotate(block):
            counts['qualified_'+role]+=1
            counts['qualified_'+role+'_'+g]+=1
            if role=='moth_predation':counts['moth_stage_'+moth_stage]+=1
            for scope in scopes:
                if role=='predator_of_predator':
                    key=(scope,g,other_group)
                    animal_edge_rows[key]+=1
                    if g in FOCAL and other_group in FOCAL:
                        animal_edge_focal[key].add((animal,partner))
                else:
                    actor_state[scope][(g,animal)] |= BITS[role]
                    if role in ('flower_visit','reported_pollinates','plant_dispersal'):
                        plant_state[scope][(g,partner)] |= BITS[role]
        if (i+1)%20==0:
            print('GLOBI_MULTIPLEX_SOURCE_PROGRESS',counts['all_interaction_rows'],
                  len(actor_state['all']),flush=True)
    require(synthetic or counts['all_interaction_rows']==SOURCE_ROWS,
            'Complete GloBI source population not processed')
    out.mkdir(parents=True)
    guild_summary=[]
    for scope in SCOPES:
        for group in GROUPS:
            levels=[b for (g,_),b in actor_state[scope].items() if g==group]
            guild_summary.append({
                'scope':scope,'predator_group':group,
                'animal_species_with_any_record':len(levels),
                'animal_species_moth_predation_recorded':sum(bool(b&1) for b in levels),
                'animal_species_flower_visit_recorded':sum(bool(b&2) for b in levels),
                'animal_species_pollinates_relation_recorded':sum(bool(b&4) for b in levels),
                'animal_species_flower_positive_either':sum(bool(b&6) for b in levels),
                'animal_species_both_moth_predation_and_flower_positive':
                    sum(bool(b&1) and bool(b&6) for b in levels),
                'animal_species_both_moth_predation_and_pollinates_relation':
                    sum(bool(b&1) and bool(b&4) for b in levels),
                'animal_species_plant_dispersal_recorded':sum(bool(b&8) for b in levels),
                'plants_with_any_flower_visit_from_group':
                    sum(g==group and bool(b&6)
                        for (g,_),b in plant_state[scope].items()),
                'focused_predator_of_predator_species_edges':sum(
                    len(edges) for (s,g,p),edges in animal_edge_focal.items()
                    if s==scope and g==group)})
    with (out/'guild_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(guild_summary[0]))
        w.writeheader();w.writerows(guild_summary)
    with gzip.open(out/'animal_species_role_flags.csv.gz','wt',newline='') as f:
        w=csv.writer(f)
        w.writerow(['scope','animal_group','animal_species',
                    'moth_predation','flower_visit','reported_pollinates',
                    'plant_dispersal','moth_kill_only_not_predation'])
        for scope,v in actor_state.items():
            for (group,species),bits in sorted(v.items()):
                w.writerow([scope,group,species]+[int(bool(bits&BITS[k])) for k in
                    ('moth_predation','flower_visit','reported_pollinates',
                     'plant_dispersal','moth_kill_only')])
    with gzip.open(out/'animal_flower_plant_positive_flags.csv.gz','wt',newline='') as f:
        w=csv.writer(f)
        w.writerow(['scope','animal_group','plant_species',
                    'flower_visit','reported_pollinates','plant_dispersal'])
        for scope,v in plant_state.items():
            for (group,species),bits in sorted(v.items()):
                w.writerow([scope,group,species]+[int(bool(bits&BITS[k])) for k in
                    ('flower_visit','reported_pollinates','plant_dispersal')])
    summary={
      'schema':'alan_full_GloBI_v011_predator_pollinator_multilayer_v1',
      'status':'SYNTHETIC_FIXTURE' if synthetic else 'FULL_GLOBI_V011_VERIFIED',
      'release':RELEASE,'source_sha256':digest,
      'all_original_interaction_rows':counts['all_interaction_rows'],
      'raw_relation_counts':{k[7:]:v for k,v in counts.items() if k.startswith('source_')},
      'qualified_role_counts':{k[10:]:v for k,v in counts.items() if k.startswith('qualified_')},
      'moth_prey_stage_counts':{k[11:]:v for k,v in counts.items() if k.startswith('moth_stage_')},
      'source_columns_used':sorted(input_fields),'guild_summary':guild_summary,
      'predator_predator_group_links':[
        {'scope':s,'predator_group':g,'prey_group':p,
         'source_interaction_rows':n,
         'distinct_named_species_pairs':len(animal_edge_focal[(s,g,p)])
              if g in FOCAL and p in FOCAL else None}
        for (s,g,p),n in sorted(animal_edge_rows.items())],
      'limitations':{
        'caterpillar_predation_not_equal_to_adult_moth_predation':True,
        'missing_GloBI_link_not_actual_absence':True,
        'recorded_pollinates_not_measured_effective_pollination':True,
        'one_species_dual_role_not_same_place_or_same_time':True,
        'group_level_dual_role_not_same_species':True,
        'no_local_predator_presence_or_impact_inferred':True,
        'no_island_vs_mainland_ALAN_effect_fitted':True,
        'study_source_links_not_independent_experiments':True}}
    (out/'source_receipt.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--chunk-size',default=120000,type=int)
    a=p.parse_args()
    result=audit(a.source,a.out,chunk_size=a.chunk_size)
    print('FULL_GLOBI_PREDATOR_POLLINATOR_MULTILAYER',json.dumps(
      {'rows':result['all_original_interaction_rows'],
       'groups':result['guild_summary']},sort_keys=True),flush=True)
