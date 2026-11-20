"""Describe the provisional catalogue; counts are not independent sample sizes."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def classify(row):
    endpoint = row['endpoint']
    if row['study_id'] in {'Giavi2021', 'Lu2026'}:
        return 'daytime_flower_visitation'
    if endpoint.endswith('_capture'):
        return 'invertebrate_capture'
    if endpoint == 'spider_abundance':
        return 'spider_abundance'
    if endpoint == 'web_prey_count':
        return 'web_prey_standing_stock'
    if endpoint.startswith('aphid_'):
        return 'aphid_abundance'
    if endpoint == 'parasitoid_mummy_abundance':
        return 'parasitoid_mummy_standing_stock'
    if endpoint.endswith('_herbivory'):
        return 'leaf_herbivory'
    if 'biomass' in endpoint:
        return 'plant_biomass_or_dead_material'
    raise ValueError(f'Unclassified endpoint: {row["effect_id"]}')


def main():
    source = ROOT / 'catalogue_evidence_view.csv'
    rows = list(csv.DictReader(source.open(encoding='utf-8-sig', newline='')))
    assert len({r['effect_id'] for r in rows}) == len(rows)
    groups = defaultdict(list)
    for row in rows:
        groups[classify(row)].append(row)
    output = []
    for stage, members in sorted(groups.items()):
        output.append(dict(
            response_group=stage, effect_rows=len(members),
            study_labels=len({r['study_id'] for r in members}),
            programme_labels=len({r['programme_id'] for r in members}),
            programme_ids=';'.join(sorted({r['programme_id'] for r in members})),
            effect_ids=';'.join(r['effect_id'] for r in members),
            final_eligibility_states=';'.join(sorted({r['final_eligibility'] for r in members}))))
    with (ROOT / 'endpoint_coverage.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    audit = dict(input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 effect_rows=len(rows), assigned_rows=sum(r['effect_rows'] for r in output),
                 study_labels=len({r['study_id'] for r in rows}),
                 programme_labels=len({r['programme_id'] for r in rows}),
                 limitation='Provisional catalogue coverage, not final eligibility or independent evidence. '
                 'Stage programme counts overlap. Unextracted and nonfinite effects are outside this table.')
    assert audit['effect_rows'] == audit['assigned_rows']
    (ROOT / 'endpoint_coverage_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    for item in output:
        print(item['response_group'], item['effect_rows'], item['programme_labels'])


if __name__ == '__main__':
    main()
