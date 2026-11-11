#!/usr/bin/env python3
"""Version-frozen, complete GloBI v0.11 lepidopteran flower-interaction scan.

Reads ALL of the interpreted 2.5 GB compressed source, not the old 61,967
bounded API rows. Retains two orientations and distinct butterfly/moth guilds.
No site-specific or effective nocturnal pollination is inferred. The output
plant-name flags need separate WCVP accepted-name reconciliation before
linkage to ALAN photos; full source names are not guaranteed WCVP accepted.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd

STABLE_RELEASE='https://zenodo.org/records/22691479'
STABLE_SHA256='670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3'
BUTTERFLY_FAMILIES={'Papilionidae','Nymphalidae','Pieridae','Lycaenidae',
                    'Riodinidae','Hesperiidae','Hedylidae'}
FORWARD={'visitsFlowersOf','pollinates'}
REVERSE={'flowersVisitedBy','pollinatedBy'}
SCOPES=('all','non_inat','non_inat_trophich','non_inat_visits')
PLANT_NAME=re.compile(r'^[A-Z][A-Za-z-]+ [a-z][a-z-]+$')
FIELDS={'interactionTypeName','sourceTaxonName','targetTaxonName',
        'sourceTaxonOrderName','targetTaxonOrderName',
        'sourceTaxonFamilyName','targetTaxonFamilyName',
        'sourceTaxonKingdomName','targetTaxonKingdomName',
        'sourceTaxonPathNames','targetTaxonPathNames',
        'sourceTaxonSpeciesName','targetTaxonSpeciesName',
        'sourceNamespace','studyTitle','sourceCitation','referenceCitation','referenceDoi',
        'decimalLatitude','decimalLongitude','eventDateUnixEpoch'}

def verify(path,expected):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''): h.update(block)
    if h.hexdigest()!=expected:
        raise ValueError('Full GloBI release SHA256 mismatch; do not claim stable source')
    return h.hexdigest()

def col(chunk,prefix,suffix):
    name=prefix+'Taxon'+suffix
    if name not in chunk:return pd.Series('',index=chunk.index,dtype='string')
    return chunk[name].fillna('').astype('string').str.strip()

def species(chunk,prefix):
    specific=col(chunk,prefix,'SpeciesName')
    names=specific.where(specific.ne(''),col(chunk,prefix,'Name'))
    return names.where(names.str.fullmatch(PLANT_NAME).fillna(False),'')

def lepidoptera(chunk,prefix):
    order=col(chunk,prefix,'OrderName').eq('Lepidoptera')
    path=col(chunk,prefix,'PathNames')
    return order|path.str.contains(r'(?:^|\|)Lepidoptera(?:\||$)',regex=True)

def plant(chunk,prefix):
    kingdom=col(chunk,prefix,'KingdomName').str.casefold()
    path=col(chunk,prefix,'PathNames').str.casefold()
    return kingdom.isin(('plantae','viridiplantae'))|path.str.contains(
        r'(?:^|\|)(?:plantae|viridiplantae|embryophyta)(?:\||$)',regex=True)

OUTPUT_COLS=['plant_name','guild','family','relation','is_inat','is_trophich',
             'insect_named_species','valid_event_coordinates','event_date_available']

def annotate(chunk):
    relations=chunk.interactionTypeName.fillna('')
    fwd=relations.isin(FORWARD);rev=relations.isin(REVERSE)
    eligible=fwd|rev
    if not eligible.any():return pd.DataFrame(columns=OUTPUT_COLS)
    d=chunk.loc[eligible].copy();r=rev.loc[eligible]
    biologic=(~r&lepidoptera(d,'source')&plant(d,'target')) | (
               r&lepidoptera(d,'target')&plant(d,'source'))
    d=d.loc[biologic].copy();r=r.loc[biologic]
    if d.empty:return pd.DataFrame(columns=OUTPUT_COLS)
    names=species(d,'target').where(~r,species(d,'source'))
    family=col(d,'source','FamilyName').where(~r,col(d,'target','FamilyName'))
    insect=species(d,'source').where(~r,species(d,'target'))
    valid=names.ne('')&family.str.fullmatch(r'[A-Z][a-z]+idae').fillna(False)
    d=d.loc[valid];names=names.loc[valid];family=family.loc[valid];insect=insect.loc[valid]
    if d.empty:return pd.DataFrame(columns=OUTPUT_COLS)
    provenance=pd.Series('',index=d.index,dtype='string')
    for field in ('sourceNamespace','studyTitle','sourceCitation','referenceCitation'):
        if field in d:provenance=provenance+' '+d[field].fillna('').astype('string')
    low=provenance.str.casefold()
    lat=pd.to_numeric(d.get('decimalLatitude',pd.Series('',index=d.index)),errors='coerce')
    lon=pd.to_numeric(d.get('decimalLongitude',pd.Series('',index=d.index)),errors='coerce')
    geo=lat.between(-90,90)&lon.between(-180,180)
    stamp=pd.to_numeric(d.get('eventDateUnixEpoch',pd.Series('',index=d.index)),errors='coerce')
    dates=stamp.gt(0)&stamp.lt(2500000000)
    return pd.DataFrame({'plant_name':names.astype(str).to_numpy(),
       'guild':family.map(lambda x:'butterfly' if x in BUTTERFLY_FAMILIES else 'moth').to_numpy(),
       'family':family.astype(str).to_numpy(),
       'relation':d.interactionTypeName.astype(str).to_numpy(),
       'is_inat':low.str.contains('inaturalist',regex=False).to_numpy(),
       'is_trophich':low.str.contains('trophich',regex=False).to_numpy(),
       'insect_named_species':insect.ne('').to_numpy(),
       'valid_event_coordinates':geo.to_numpy(),
       'event_date_available':dates.to_numpy()})

def extract(source,out,*,expected_sha=STABLE_SHA256,allow_fixture=False,chunk_size=160000):
    source,out=Path(source),Path(out)
    if out.exists():raise ValueError('Refuse to overwrite source-derived result')
    if chunk_size<100:raise ValueError('Invalid chunk size')
    if expected_sha!=STABLE_SHA256 and not allow_fixture:
        raise ValueError('Cannot silently replace release SHA with synthetic data')
    digest=verify(source,expected_sha)
    header=set(pd.read_csv(source,compression='gzip',nrows=0).columns)
    required={'interactionTypeName','sourceTaxonName','targetTaxonName',
              'sourceTaxonFamilyName','targetTaxonFamilyName'}
    if not required<=header:raise ValueError('GloBI interpreted schema drift')
    for side in ('source','target'):
        if not ({side+'TaxonOrderName',side+'TaxonPathNames'}&header):
            raise ValueError('Missing Lepidoptera verification: '+side)
        if not ({side+'TaxonKingdomName',side+'TaxonPathNames'}&header):
            raise ValueError('Missing plant kingdom verification: '+side)
    cols=FIELDS&header
    counts=Counter();flags=defaultdict(set);family_counts=Counter()
    for i,chunk in enumerate(pd.read_csv(source,compression='gzip',
              usecols=lambda x:x in cols,dtype='string',
              chunksize=chunk_size,low_memory=False)):
        counts['all_interaction_rows']+=len(chunk)
        valid=annotate(chunk)
        counts['qualified_flower_interactions']+=len(valid)
        if valid.empty:continue
        family_counts.update(valid.family.tolist())
        for guild in ('moth','butterfly'):
            t=valid.loc[valid.guild.eq(guild)]
            counts[guild+'_georeferenced_rows']+=int(t.valid_event_coordinates.sum())
            counts[guild+'_dated_rows']+=int(t.event_date_available.sum())
            counts[guild+'_georef_plus_date_rows']+=int(
                (t.valid_event_coordinates&t.event_date_available).sum())
        for name,guild,inat,trophic,rel,known in valid[
            ['plant_name','guild','is_inat','is_trophich','relation',
             'insect_named_species']].itertuples(index=False,name=None):
            counts[guild+'_qualified_rows']+=1
            if known:counts['insect_species_identified_rows']+=1
            tiers=['all']
            if not inat:
                tiers.append('non_inat')
                if not trophic:tiers.append('non_inat_trophich')
                if rel in ('visitsFlowersOf','flowersVisitedBy'):
                    tiers.append('non_inat_visits')
            for scope in tiers:
                flags[name].add((scope,guild,'species_insect' if known else 'family_only'))
        if (i+1)%25==0:
            print('GLOBI_FULL_PASS',i+1,counts['all_interaction_rows'],
                  counts['qualified_flower_interactions'],flush=True)
    out.mkdir(parents=True)
    flds=['plant_name']+[f'{s}_{g}_{level}' for s in SCOPES
            for g in ('moth','butterfly') for level in ('any_family','species_insect')]
    with gzip.open(out/'GloBI_v011_positive_species_flags.csv.gz','wt',
                   encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=flds);w.writeheader()
        for name,records in sorted(flags.items()):
            row={'plant_name':name}
            for s in SCOPES:
                for g in ('moth','butterfly'):
                    row[f'{s}_{g}_any_family']=int(any((s,g,t) in records
                        for t in ('family_only','species_insect')))
                    row[f'{s}_{g}_species_insect']=int((s,g,'species_insect') in records)
            w.writerow(row)
    by_tier={f'{s}_{g}_{level}':sum(any((s,g,t) in r for t in
               (('family_only','species_insect') if level=='any_family'
                else ('species_insect',))) for r in flags.values())
             for s in SCOPES for g in ('moth','butterfly')
             for level in ('any_family','species_insect')}
    report={'status':('FULL_GLOBI_V011_SOURCE_COMPLETE' if expected_sha==STABLE_SHA256
                      else 'SYNTHETIC_FIXTURE_NOT_REAL'),
            'source_url':STABLE_RELEASE,'sha256':digest,
            'all_interaction_rows':int(counts['all_interaction_rows']),
            'qualified_flower_interaction_rows':int(counts['qualified_flower_interactions']),
            'moth_qualified_rows':int(counts['moth_qualified_rows']),
            'butterfly_qualified_rows':int(counts['butterfly_qualified_rows']),
            'plant_names_with_any_positive_claim':len(flags),
            'insect_species_identified_rows':int(counts['insect_species_identified_rows']),
            'taxon_by_evidence':by_tier,'families_top25':family_counts.most_common(25),
            'georeferenced_and_dated_rows':{g:int(counts[g+'_georef_plus_date_rows'])
                                            for g in ('moth','butterfly')},
            'butterfly_families':sorted(BUTTERFLY_FAMILIES),
            'relations':sorted(FORWARD|REVERSE),
            'source_taxon_crosswalk_to_WCVP_still_required':True,
            'missing_GloBI_not_pollinator_absence':True,
            'flower_visits_not_effective_nocturnal_pollination':True,
            'source_columns_used':sorted(cols)}
    (out/'source_receipt.json').write_text(
        json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--chunk-size',type=int,default=160000)
    a=p.parse_args()
    print(json.dumps(extract(a.source,a.out,chunk_size=a.chunk_size),ensure_ascii=False))
if __name__=='__main__':main()
