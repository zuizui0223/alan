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
           effects_extracted=0,
           unresolved=['Group labels differ across response sheets; resolve physical lighting assignment from layout.',
                       'Long-format calendar year differs from the primary field year; preserve rather than repair.',
                       'Indoor phase eligibility and illumination independence require Methods review.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
