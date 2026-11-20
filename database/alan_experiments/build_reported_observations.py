"""Link reported statistics without converting them into pooled effect sizes."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCES = [
    ('Boom2020', 'boom2020_sources/reported_arm_summaries.csv', 'arm_summary', 'Boom2020_LON', 16),
    ('Liu2022', 'liu2022_sources/reported_tests.csv', 'likelihood_ratio_test', 'Liu2022_Konstanz', 15),
    ('Wilson2021', 'wilson2021_sources/reported_coefficients.csv', 'reported_model_coefficient', None, 7),
    ('Nunez2026', 'search/nunez2026_model_contrasts.csv', 'reported_model_odds_ratio', 'Nunez2026_MCH', 2),
]

def read(name):
    with (ROOT / name).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def write(name, rows):
    with (ROOT / name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def build():
    studies = {r['study_id']: r for r in read('studies.csv')}
    experiments = {r['experiment_id']: r for r in read('experiments.csv')}
    observations, links, hashes = [], {}, {}
    for study_id, source, kind, experiment_id, expected_count in SOURCES:
        raw = read(source)
        assert len(raw) == expected_count, source
        hashes[source] = hashlib.sha256((ROOT / source).read_bytes()).hexdigest()
        for line, row in enumerate(raw, 2):
            eid = experiment_id
            if study_id == 'Wilson2021':
                eid = 'Wilson2019_capture' if row['endpoint'].endswith('_trap_count') else 'Wilson2019_reproduction'
            study, experiment = studies[study_id], experiments[eid]
            assert experiment['study_id'] == study_id
            assert experiment['programme_id'] == study['programme_id']
            assert row.get('doi', study['doi']) == study['doi']
            payload = json.dumps(row, ensure_ascii=False, sort_keys=True)
            oid = study_id + '_' + hashlib.sha256(payload.encode()).hexdigest()[:16]
            observations.append(dict(observation_id=oid, study_id=study_id, experiment_id=eid,
                programme_id=study['programme_id'], doi=study['doi'], statistic_type=kind,
                endpoint=row['endpoint'], source_file=source, source_csv_line=line,
                source_location=row['source_location'], source_values_json=payload,
                pooling_status='not_in_effect_catalogue_not_poolable_as_registered',
                dependence_group=eid, final_eligibility='not_adjudicated'))
            assert json.loads(observations[-1]['source_values_json']) == row
            links[(study_id, eid)] = dict(study_id=study_id, experiment_id=eid,
                programme_id=study['programme_id'], doi=study['doi'])
    assert len({r['observation_id'] for r in observations}) == len(observations)
    write('reported_observations.csv', observations)
    write('reported_experiment_studies.csv', list(links.values()))
    audit = dict(observation_rows=len(observations), studies=len(SOURCES),
        experiment_links=len(links), source_sha256=hashes, added_meta_effects=0,
        limitation='Reported arms, tests and coefficients are not independent effects. Shared dependence groups are labels, not estimated covariance. Source values remain strings without scale conversion. This is partial registration, not exhaustive extraction or final adoption.')
    (ROOT / 'reported_observations_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(audit, indent=2))

if __name__ == '__main__':
    build()
