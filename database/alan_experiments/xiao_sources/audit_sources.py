"""Audit recovered code/supplement and exact-name support; no effect extraction."""
from pathlib import Path
import hashlib
import json
import pandas as pd

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
CACHE=ROOT.parent/'gbif_annual/alan_experiment_sources/xiao2025'
z=json.loads((CACHE/'zenodo_metadata.json').read_text(encoding='utf-8'))
f=json.loads((CACHE/'figshare_metadata.json').read_text(encoding='utf-8'))
manifest=[]
for entry in z['files']:
    path=CACHE/entry['key']
    digest=hashlib.md5(path.read_bytes()).hexdigest()
    assert entry['checksum']=='md5:'+digest
    manifest.append(dict(file=path.name,bytes=path.stat().st_size,md5=digest,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),url=entry['links']['self']))
for entry in f['files']:
    path=CACHE/entry['name']
    digest=hashlib.md5(path.read_bytes()).hexdigest()
    assert entry['computed_md5']==digest
    manifest.append(dict(file=path.name,bytes=path.stat().st_size,md5=digest,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),url=entry['download_url']))
pd.DataFrame(manifest).to_csv(OUT/'recovered_sources.csv',index=False)
d=json.loads((CACHE/'files_metadata.json').read_text(encoding='utf-8'))
required=[]
for entry in d['_embedded']['stash:files']:
    path=CACHE/entry['path']
    valid=path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==entry['digest']
    required.append(dict(file=path.name,expected_bytes=entry['size'],expected_sha256=entry['digest'],
        retrieved_and_verified=valid))
pd.DataFrame(required).to_csv(OUT/'dryad_required_files.csv',index=False)

species=pd.read_csv(OUT/'species_as_reported.csv')
assert len(species)==12 and species.source_species.is_unique
eiv_path=ROOT.parent/'gbif_annual/EIV_Europe.csv'
assert hashlib.sha256(eiv_path.read_bytes()).hexdigest()=='ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba'
eiv=pd.read_csv(eiv_path,index_col=0,encoding='latin1')
assert eiv.index.is_unique
wet=pd.read_csv(ROOT/'results/all_species_forest/wet_reproduction/species.csv')
assert len(wet)==363 and int(wet.candidate.sum())==58 and wet.accepted.nunique()==361
species['exact_EIV_name']=species.source_species.isin(eiv.index)
species['EIV_MOIST_raw']=species.source_species.map(eiv.MOIST)
species['EIV_MOIST']=pd.to_numeric(species.EIV_MOIST_raw,errors='coerce')
species['wet_macro_source_name_match']=species.source_species.isin(wet.species)
species['wet_macro_accepted_name_match']=species.source_species.isin(wet.accepted)
species['wet_macro_candidate_match']=species.source_species.isin(wet.loc[wet.candidate.eq(True),'accepted'])
species['wet_macro_reference_match']=species.source_species.isin(wet.loc[wet.candidate.eq(False),'accepted'])
species['experimental_effect_available']=False
species.to_csv(OUT/'species_support.csv',index=False)
audit=dict(doi='10.1098/rspb.2024.2811',data_doi='10.5061/dryad.1rn8pk143',
    zenodo_record=14000092,figshare_article=29118899,
    recovered_files=len(manifest),raw_workbooks_verified=sum(r['retrieved_and_verified'] for r in required if r['file'].endswith('.xlsx')),
    raw_workbooks_required=2,effects_extracted=0,main_catalogue_changed=False,
    download_recheck=json.loads((CACHE/'download_recheck.json').read_text(encoding='utf-8')),
    species=12,exact_numeric_moisture_species=int(species.EIV_MOIST.notna().sum()),
    exact_nonnumeric_moisture_species=int((species.exact_EIV_name & species.EIV_MOIST.isna()).sum()),
    wet_comparison_labels=len(wet),wet_candidate_labels=int(wet.candidate.sum()),
    exact_wet_candidate_species=int(species.wet_macro_candidate_match.sum()),
    exact_wet_reference_species=int(species.wet_macro_reference_match.sum()),
    species_join='exact binomial only; no synonym imputation',
    moisture_is_species_affinity_not_greenhouse_condition=True,
    constraints=['Code and supplement are not raw observation tables',
        'Three paired light compartments for growth; pots are subsamples',
        'Separate feeding assay one compartment pair is not replicated at the lighting level',
        'Species overlap cannot identify the response direction without raw effects'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(species[['source_species','EIV_MOIST','wet_macro_candidate_match','wet_macro_reference_match']].to_string(index=False))
print(json.dumps(audit,indent=2))
