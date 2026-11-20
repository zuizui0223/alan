"""Preserve author summary cells; do not infer endpoint units or independence."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
source = ROOT / 'Mydata_light_summary.csv'
raw = source.read_bytes()
assert hashlib.md5(raw).hexdigest() == 'b116dc5eababf963e7ea03c352ee1f80'
with source.open(encoding='utf-8', newline='') as handle:
    reader = csv.DictReader(handle)
    fields = reader.fieldnames
    rows = list(reader)
endpoints = [x.removesuffix('_mean') for x in fields if x.endswith('_mean')]
assert len(rows) == 14 and len(endpoints) == 17
assert len({(r['Time'], r['Treatment']) for r in rows}) == len(rows)
assert all({r['Treatment'] for r in rows if r['Time'] == t} == {'ALAN', 'Control'}
           for t in {r['Time'] for r in rows})
out = []
for i, row in enumerate(rows, start=1):
    for endpoint in endpoints:
        mean, se = row[endpoint + '_mean'], row[endpoint + '_se']
        assert (mean == 'NA') == (se == 'NA')
        if se != 'NA':
            assert float(se) >= 0
        out.append(dict(source_data_row=i, time_label=row['Time'],
                        treatment=row['Treatment'], endpoint_code=endpoint,
                        mean_original=mean, se_original=se, independent_n='',
                        unit='', time_unit='', status='summary_only_not_pooled'))
with (ROOT / 'reported_summary.csv').open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(out[0]))
    writer.writeheader()
    writer.writerows(out)
audit = dict(article_doi='10.1111/gcb.71029', dataset_doi='10.6084/m9.figshare.31288969.v1',
             dataset_author='Min Liu', license='CC BY 4.0', retrieval_date='2026-11-21',
             api_url='https://api.figshare.com/v2/articles/31288969',
             download_url='https://ndownloader.figshare.com/files/61724725',
             source_bytes=len(raw), source_md5=hashlib.md5(raw).hexdigest(),
             source_sha256=hashlib.sha256(raw).hexdigest(), source_rows=len(rows),
             endpoint_codes=endpoints, extracted_rows=len(out),
             missing_mean_se_pairs=sum(r['mean_original'] == 'NA' for r in out),
             zero_mean_rows=sum(r['mean_original'] != 'NA' and float(r['mean_original']) == 0 for r in out),
             independent_n_confirmed=False, endpoint_units_confirmed=False,
             covariance_available=False, new_meta_effects=0,
             primary_text_recheck='PMC browser challenge; Europe PMC XML HTTP 500; no bypass')
(ROOT / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n', encoding='utf-8')
print(json.dumps(audit, indent=2))
