"""Re-extract an original experiment, not a pre-existing meta-analysis table."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from urllib.parse import quote
from urllib.request import urlopen
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('source', ROOT/'scripts/analysis/alan_true_parkinson2020_aquatic_spider_20261118.py')
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def write(name, rows):
    with (OUT/name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    cache = Path(sys.argv[1])
    cache.mkdir(parents=True, exist_ok=True)
    manifest = []
    paths = []
    for name, expected in source.SOURCE_BLOB_SHA1.items():
        url = f'https://raw.githubusercontent.com/{source.SOURCE_REPO}/{source.SOURCE_COMMIT}/{quote(name)}'
        path = cache/name
        if not path.exists():
            path.write_bytes(urlopen(url, timeout=60).read())
        if source.git_blob_sha1(path) != expected:
            raise ValueError('source hash mismatch')
        paths.append(path)
        manifest.append(dict(filename=name,url=url,git_blob=expected,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    records, audit = source.grouped_source(*paths)
    write('source_manifest.csv', manifest)
    (OUT/'source_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    endpoints = [('aquatic_pan','aquatic_origin_capture','single_overnight','primary'),
                 ('terrestrial_pan','terrestrial_origin_capture','single_overnight','primary'),
                 ('spider_final','spider_abundance','2016-08-04','primary'),
                 ('prey_final','web_prey_count','2016-08-04','primary'),
                 ('spider_mean_3dates','spider_abundance','three_date_mean','supplementary'),
                 ('prey_in_web_mean_3dates','web_prey_count','three_date_mean','supplementary')]
    for r in records:
        r['spider_final'] = r['spider_3dates'][-1]
        r['prey_final'] = r['web_prey_3dates'][-1]
    mask = np.array([r['alan'] for r in records])
    matrix = np.array([[r[k] for k,_,_,_ in endpoints] for r in records], dtype=float)
    rows = []
    for j,(key,endpoint,timepoint,role) in enumerate(endpoints):
        lit,dark = matrix[mask,j],matrix[~mask,j]
        rows.append(dict(effect_id='P2020_'+key,study_id='Parkinson2020',experiment_id='Parkinson2016_lake',
            programme_id='Parkinson2016',control_group_id='Parkinson2016_dark10',endpoint=endpoint,
            stage='insect_predator',timepoint=timepoint,analysis_role=role,design='independent_arms',
            independent_unit='mesocosm',dispersion_unit='mesocosm',mean_lit=lit.mean(),mean_dark=dark.mean(),
            dispersion_type='SD',dispersion_lit=lit.std(ddof=1),dispersion_dark=dark.std(ddof=1),n_lit=10,n_dark=10,
            ci_level='',ci_distribution='',ci_df_lit='',ci_df_dark='',source_location='raw_data: source_manifest.csv',
            digitized=0,extraction_date='2026-11-20',verification_status='independent_agent_computational_reextraction_verified; human_review_pending',
            note='same 20 units across endpoints; fish assignment unresolved; no plant outcome'))
    write('effects.csv', rows)
    # Multivariate delta covariance uses SAME units for both outcomes within each arm.
    v = sum(np.cov(matrix[m],rowvar=False,ddof=1)/m.sum()/np.outer(matrix[m].mean(axis=0),matrix[m].mean(axis=0)) for m in (mask,~mask))
    write('effect_covariance.csv',[dict(effect_id_1=rows[i]['effect_id'],effect_id_2=rows[j]['effect_id'],sampling_covariance=v[i,j],
        method='within_arm_unit_sample_covariance_delta') for i in range(len(rows)) for j in range(len(rows))])
    print('Extracted 1 experiment, 20 units, 4 primary and 2 supplementary endpoints; covariance saved')


if __name__ == '__main__':
    main()
