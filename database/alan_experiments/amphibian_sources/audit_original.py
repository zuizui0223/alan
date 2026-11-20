"""Check archival identity and group labels before estimating lighting effects."""
from pathlib import Path
import hashlib
import json
import pandas as pd

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/amphibian2018'
raw=(CACHE/'original.xlsx').read_bytes()
zen=json.loads((CACHE/'zenodo_metadata.json').read_text(encoding='utf-8'))['files'][0]
dry=json.loads((CACHE/'dryad_files.json').read_text())['_embedded']['stash:files'][0]
assert len(raw)==zen['size']==dry['size']==40950
assert hashlib.md5(raw).hexdigest()==dry['digest']=='4c54998fd893d7a766b34ffb9ea2af8b'
assert zen['checksum']=='md5:'+dry['digest']
book=pd.ExcelFile(CACHE/'original.xlsx')
pd.DataFrame([dict(sheet=s,rows=len(pd.read_excel(book,s)),columns=len(pd.read_excel(book,s).columns))
              for s in book.sheet_names]).to_csv(OUT/'workbook_inventory.csv',index=False)
colon=pd.read_excel(book,'AmphibianPipeSample')
toad=pd.read_excel(book,'toadmesocosm')
assert colon.Mesocosm.is_unique and toad.Mesocosm.is_unique
assert len(colon)==36 and len(toad)==18
check=toad[['Mesocosm','LightTreat','Block']].merge(
    colon[['Mesocosm','LightTreat','Block']],on='Mesocosm',validate='one_to_one',
    suffixes=('_toad','_colonization'))
mapping={'Light':'Present','No Light':'Absent'}
assert check.LightTreat_toad.map(mapping).eq(check.LightTreat_colonization).all()
check['group_labels_match']=check.Block_toad.str.upper().eq(check.Block_colonization.str.upper())
check.to_csv(OUT/'group_label_crosswalk.csv',index=False)
# A renamed block is not a contradictory partition. Verify every pair of
# shared mesocosms, not just a guessed alphabetical correspondence.
block_map=check[['Block_toad','Block_colonization']].drop_duplicates()
assert len(block_map)==9
assert block_map.Block_toad.is_unique and block_map.Block_colonization.is_unique
same_toad=check.Block_toad.to_numpy()[:,None]==check.Block_toad.to_numpy()[None,:]
same_colon=check.Block_colonization.to_numpy()[:,None]==check.Block_colonization.to_numpy()[None,:]
assert (same_toad==same_colon).all()
block_map.to_csv(OUT/'verified_block_label_map.csv',index=False)
# Export source-labelled joins only. This does not certify physical lamp layout.
join_audit=[]
for sheet in ['lengthofmeta','toadmesocosmwaterquality','InvertPipeSample']:
    frame=pd.read_excel(book,sheet)
    assert frame.Mesocosm.is_unique
    joined=frame.merge(colon[['Mesocosm','Block','LightTreat']],on='Mesocosm',
                       validate='one_to_one',suffixes=('_source','_colonization'),indicator=True)
    assert len(joined)==len(frame) and joined['_merge'].eq('both').all()
    assert joined.LightTreat_source.replace(mapping).eq(joined.LightTreat_colonization).all()
    if 'Block_source' in joined:
        mapped=joined.Block_source.str.upper().map(dict(zip(block_map.Block_toad,block_map.Block_colonization)))
        assert mapped.eq(joined.Block_colonization).all()
    join_audit.append(dict(sheet=sheet,rows=len(frame),unique_mesocosm_join=True,
                           illumination_labels_agree=True,physical_lamp_layout_verified=False))
pd.DataFrame(join_audit).to_csv(OUT/'verified_mesocosm_join_support.csv',index=False)

groups=colon.groupby(['Block','LightTreat']).size().reset_index(name='mesocosms')
assert groups.mesocosms.eq(4).all() and len(groups)==9
groups.to_csv(OUT/'colonization_sheet_groups.csv',index=False)
long=pd.read_excel(book,'toadlabgrowthreshape1')
audit=dict(source_doi='10.5061/dryad.pq54gh5',zenodo_record='4932781',
           original_filename=dry['path'],md5=dry['digest'],mirror_bytes_identical=True,
           sheets=len(book.sheet_names),outdoor_mesocosms=36,toad_mesocosms=18,
           colonization_sheet_light_groups=4,colonization_sheet_dark_groups=5,
           mismatching_toad_vs_colonization_group_labels=int((~check.group_labels_match).sum()),
           long_format_years=sorted(pd.to_datetime(long.Date).dt.year.unique().tolist()),
           block_map_bijective=True,shared_mesocosm_partition_identical=True,
           shared_mesocosm_pairs_checked=int(same_toad.size),physical_lamp_layout_verified=False,
           effects_extracted=0,
           unresolved=['Block names are a verified bijective renaming on shared mesocosms; physical lamp layout remains unverified.',
                       'Long-format calendar year differs from the primary field year; preserve rather than repair.',
                       'Indoor phase eligibility and illumination independence require Methods review.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
