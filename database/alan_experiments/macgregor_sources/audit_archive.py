"""Verify recovered code archive and explicitly inventory missing data inputs."""
from pathlib import Path
import hashlib
import json
import re
import zipfile

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[2].parent/'gbif_annual/alan_experiment_sources/macgregor2019'
archive=CACHE/'original.zip'
meta=json.loads((CACHE/'zenodo_metadata.json').read_text())
expected=meta['files'][0]['checksum']
assert expected=='md5:'+hashlib.md5(archive.read_bytes()).hexdigest()
with zipfile.ZipFile(archive) as z:
    files=[n for n in z.namelist() if not n.endswith('/')]
    code=z.read(next(n for n in files if n.endswith('/FinalScripts.R'))).decode('utf-8-sig')
inputs=sorted(set(re.findall(r'read\.csv\("([^"\n]+)"',code)))
assert len(files)==6 and all(n.endswith('.R') for n in files)
assert inputs and not any(n.lower().endswith('.csv') for n in files)
out=dict(paper_doi='10.1002/ecs2.2550',date='2026-11-20',
         discovery='supplementary targeted search: Macgregor 2019 street lighting moth pollination seed set Dryad data',
         search_origin='outside existing indexed search; not WoS or Scopus',
         primary_pdf_url='https://eprints.whiterose.ac.uk/id/eprint/149121/1/Macgregor_et_al_2019_Ecosphere.pdf',
         primary_pdf_sha256=hashlib.sha256((CACHE/'primary.pdf').read_bytes()).hexdigest(),
         zenodo_record=1743917,archive_checksum=expected,archive_files=files,
         csv_inputs_referenced=inputs,csv_inputs_present=False,
         missing_sourced_helper='OverdispersalFunction.R',
         public_roots_checked=['https://github.com/CallumJMacgregor/BioassayPaper/tree/master',
                               'https://github.com/CallumJMacgregor/Chapter-3/tree/master'],
         public_roots_raw_csv_present=False,
         supplement_url='https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecs2.2550&file=ecs22550-sup-0001-AppendixS1.pdf',
         supplement_status=403,extracted_effects=0,main_catalogue=False,
         units=dict(physical_plots=5,rotating_runs=4,individual_female_plants=169,flowers=417,
                    successful_flowers=194,damaged_successful_flowers=46,conditional_seed_sample=148),
         blockers=['Raw plant/plot/run linkage unavailable','Repeated plants and rotating plots require dependence model',
                   'Seed outcomes are conditional; damaged capsules not zero-filled',
                   'No signed lighting coefficient and SE extracted from primary summary'])
(OUT/'audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(archive_files=len(files),missing_csv_inputs=inputs,effects=0)))
