"""Exact-name EIV join for individually identified plant endpoints.

Species affinity is not experimental soil moisture or site water distance.
Community-level effects and genus-only identities are not imputed.
"""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'site_context'
src=ROOT.parents[2]/'gbif_annual/EIV_Europe.csv'
eiv=pd.read_csv(src,index_col=0,encoding='latin1')
assert hashlib.sha256(src.read_bytes()).hexdigest()=='ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba'
assert eiv.index.is_unique
moist=pd.to_numeric(eiv.MOIST,errors='coerce')
rows=[]
for table in ['hey_sources/hey2020_biomass_effects.csv','xue_sources/monoculture_effects.csv']:
    d=pd.read_csv(ROOT/table)
    for r in d.itertuples():
        rows.append(dict(effect_id=r.effect_id,study_id=r.study_id,source_species=r.species,
                         taxon_status='species_as_reported',source_table=table))
d=pd.read_csv(ROOT/'bennie_sources/final_effects.csv')
for r in d[d.endpoint=='lotus_flowerhead_count'].itertuples():
    rows.append(dict(effect_id=r.effect_id,study_id='Bennie2015',source_species='Lotus pedunculatus',
                     taxon_status='species_as_reported',source_table='bennie_sources/final_effects.csv'))
for endpoint,species,status in [('psychotria_herbivory','Psychotria asiatica','species_as_reported'),
                               ('aporosa_herbivory','Aporosa','genus_only_conflicting_source_species')]:
    rows.append(dict(effect_id='Guo2025_'+endpoint,study_id='Guo2025',source_species=species,
                     taxon_status=status,source_table='hongkong_sources/final_effects_paired.csv'))
crosswalk=pd.read_csv(ROOT/'grassland_sources/plant_taxon_crosswalk.csv').fillna('')
grass=pd.read_csv(ROOT/'grassland_sources/all_final_contrasts.csv')
assert crosswalk.material.is_unique and set(grass.material)==set(crosswalk.material)
grass=grass.merge(crosswalk,on='material',validate='many_to_one')
for r in grass[grass.taxon_status!='non_species_endpoint'].itertuples():
    rows.append(dict(effect_id=r.effect_id,study_id=r.study_id,source_species=r.source_species,
                     taxon_status=r.taxon_status,source_table='grassland_sources/all_final_contrasts.csv'))
d=pd.DataFrame(rows)
assert d.effect_id.is_unique
d['EIV_MOIST']=d.source_species.map(moist)
d['EIV_MOIST_raw']=d.source_species.map(eiv.MOIST)
d['exact_taxon_in_EIV']=d.source_species.isin(eiv.index)
d['join_method']=d.EIV_MOIST.notna().map({True:'exact_binomial',False:'no_value_no_imputation'})
d.loc[d.exact_taxon_in_EIV & d.EIV_MOIST.isna(),'join_method']='exact_taxon_non_numeric_indicator'
d.loc[d.taxon_status!='species_as_reported','join_method']='withheld_'+d.loc[d.taxon_status!='species_as_reported','taxon_status']
assert d.loc[d.taxon_status!='species_as_reported','EIV_MOIST'].isna().all()
d['species_affinity_not_site_condition']=True
d.to_csv(OUT/'plant_endpoint_moisture.csv',index=False)
taxa=d[['source_species','taxon_status','EIV_MOIST_raw','EIV_MOIST','exact_taxon_in_EIV','join_method']].drop_duplicates()
taxa.to_csv(OUT/'plant_species_moisture.csv',index=False)
catalogue=pd.read_csv(ROOT/'effect_catalogue.csv')
support=d.merge(catalogue[['effect_id','programme_id','endpoint','lnRR','sampling_variance']],
                on='effect_id',how='left',validate='one_to_one')
support['finite_effect_with_moisture']=support[['EIV_MOIST','lnRR','sampling_variance']].notna().all(axis=1)
support.to_csv(OUT/'plant_moisture_effect_support.csv',index=False)
finite=support[support.finite_effect_with_moisture]
by_programme=finite.groupby('programme_id').agg(effect_rows=('effect_id','size'),
    plant_species=('source_species','nunique'),moist_min=('EIV_MOIST','min'),moist_max=('EIV_MOIST','max'))
by_programme.to_csv(OUT/'plant_moisture_programme_support.csv')
summary=dict(source_file='external_cache/gbif_annual/EIV_Europe.csv',
    sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
    source_role='same author-distributed CSV snapshot as macro trait analysis; version identified by SHA256, not assumed latest Zenodo release',
    source_url='https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv',
    publication_doi='10.1111/jvs.13168',
    provenance_record='results/all_species_forest/decline_traits/README.md',
    endpoint_rows=len(d),unique_taxon_labels=len(taxa),
    exact_moisture_taxa=int(taxa.EIV_MOIST.notna().sum()),
    finite_moisture_effect_rows=len(finite),
    finite_moisture_species=finite.source_species.nunique(),
    finite_moisture_programmes=finite.programme_id.nunique(),
    programmes_with_multiple_moisture_species=int(by_programme.plant_species.gt(1).sum()),
    excluded_from_this_join='community biomass and mixed-flora visits have no single plant species; animal habitat not inferred from plant EIV',
    not_done=['synonym matching','community weighted moisture','site water distance','moisture meta-regression'])
(OUT/'plant_moisture_manifest.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(taxa.to_string(index=False))
print(json.dumps(summary,indent=2))
