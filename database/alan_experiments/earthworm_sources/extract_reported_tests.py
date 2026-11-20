"""Extract original supplement ANOVA rows, not signed meta-analysis effects."""
from pathlib import Path
from lxml import etree
import zipfile
import csv
import json
import hashlib

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3]/'gbif_annual/alan_experiment_sources'
source = CACHE/'earthworm2025/442_2025_5750_MOESM1_ESM.docx'
root = etree.fromstring(zipfile.ZipFile(source).read('word/document.xml'))
ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
rows=[]
for tr in root.xpath('//w:tbl/w:tr',namespaces=ns):
    cells=[''.join(c.xpath('.//w:t/text()',namespaces=ns)).strip() for c in tr.xpath('./w:tc',namespaces=ns)]
    if cells[0].startswith('Day '): day=int(cells[0].split()[1])
    elif cells[0] in {'Treatment','Treatment:Survival'}:
        rows.append(dict(day=day,term=cells[0],numerator_df=int(cells[1]),sum_sq=float(cells[2]),
                         mean_sq=float(cells[3]),F=float(cells[4]),reported_p=float(cells[5]),
                         endpoint='author_cumulative_soil_respiration',source_location='Supplement Table S1',
                         status='reported_test_not_signed_effect'))
assert len(rows)==4 and {r['day'] for r in rows}=={29,36}
with (OUT/'reported_respiration_tests.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit=dict(doi='10.1007/s00442-025-05750-z',pmcid='PMC12205020',
           supplement_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12205020/supplementaryFiles',
           supplement_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
           primary_xml_sha256=hashlib.sha256((CACHE/'fulltexts_20261120/PMC12205020.xml').read_bytes()).hexdigest(),
           table_rows=4, main_effects_added=0, primary_read_count_increment=0,
           raw_data='author_request_not_recovered', supplementary_content='8figures and TableS1; no raw mesocosm data')
(OUT/'source_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print('Four reported tests extracted; no inferred effect magnitude or sign')
