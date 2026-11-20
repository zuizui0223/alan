"""Transcribe author-model coefficients; not arm-level meta-analysis effects."""
from pathlib import Path
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET
import pandas as pd

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/middendorp2025'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with zipfile.ZipFile(CACHE / 'supplement.docx') as z:
    root = ET.fromstring(z.read('word/document.xml'))
tables = []
for table in root.findall('.//w:tbl', NS):
    tables.append([[''.join(t.text or '' for t in c.findall('.//w:t', NS))
                    for c in r.findall('w:tc', NS)] for r in table.findall('w:tr', NS)])
assert len(tables[0]) == 6 and len(tables[0][0]) == 11
rows = []
for values in tables[0][1:]:
    for column, term in [(2, 'Period (Night)'), (4, 'Group (ALAN)'), (6, 'Period x Group')]:
        estimate, dispersion = values[column].split('\u00b1')
        rows.append(dict(response=values[1], family=values[0], term=term,
            estimate=float(estimate), reported_plus_minus=float(dispersion),
            reported_p=values[column+1], source='Supplement Table S1',
            status='author_model_coefficient_coding_and_joint_covariance_unresolved'))
pd.DataFrame(rows).to_csv(OUT / 'reported_model_coefficients.csv', index=False)
# Preserve the conflicting labels exactly; do not silently swap them.
labels = pd.DataFrame(dict(date=tables[1][0][1:], supplement_period=tables[1][1][1:],
    main_methods_period=['Night', 'Night', 'Day', 'Day']))
assert labels.supplement_period.ne(labels.main_methods_period).all()
labels.to_csv(OUT / 'date_period_discrepancy.csv', index=False)
manifest=[]
for filename,url in [
 ('primary.pdf','https://static1.squarespace.com/static/67c918b80aa879635d104a6f/t/68a756aaea40e3434bbdaaf4/1755797162366/1-s2.0-S0006320725004768-main%2B%281%29.pdf'),
 ('supplement.docx','https://ars.els-cdn.com/content/image/1-s2.0-S0006320725004768-mmc1.docx')]:
    raw=(CACHE/filename).read_bytes()
    manifest.append(dict(filename=filename,url=url,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.csv',index=False)
audit=dict(doi='10.1016/j.biocon.2025.111439',reported_coefficient_rows=15,
    meta_effects_added=0,plots_total=20,plots_per_arm=10,design='alternating spatial assignment',
    cautions=['Repeated dates, bait stations and images are not independent illumination units.',
              'Table S2 dates and diel labels disagree with main Methods.',
              'Table S1 reference-level wording does not resolve contrast coding.',
              'Coefficient covariance unavailable; do not sum Group and interaction variances.',
              'Log-transformed response coefficient is not automatically arithmetic-mean lnRR.',
              'Protein-removal beta-logit coefficient is not lnRR.',
              'Dead insect bait removal is scavenging/resource removal, not live moth attack.',
              'Original plot-level data described as available on request; no author contact made.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print('15 reported coefficients; 4 date/period discrepancies; zero catalogue additions.')
