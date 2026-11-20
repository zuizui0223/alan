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
    water_source = ROOT / 'site_context/monthly_water_context.csv'
    water = list(csv.DictReader(water_source.open(encoding='utf-8-sig', newline='')))
    water_programmes = {r['programme_id'] for r in water}
    valid_water_programmes = {r['programme_id'] for r in water
                             if r['water_fraction_if_coverage_ge95'] != ''}
    assert len({r['effect_id'] for r in rows}) == len(rows)
    groups = defaultdict(list)
    for row in rows:
        groups[classify(row)].append(row)
    output = []
    for stage, members in sorted(groups.items()):
        programmes = {r['programme_id'] for r in members}
        eiv = [r for r in members if r['EIV_MOIST'] != '']
        moisture = [float(r['EIV_MOIST']) for r in eiv]
        output.append(dict(
            response_group=stage, effect_rows=len(members),
            study_labels=len({r['study_id'] for r in members}),
            programme_labels=len({r['programme_id'] for r in members}),
            programme_ids=';'.join(sorted({r['programme_id'] for r in members})),
            effect_ids=';'.join(r['effect_id'] for r in members),
            final_eligibility_states=';'.join(sorted({r['final_eligibility'] for r in members})),
            rows_with_numeric_species_EIV=len(eiv),
            species_with_numeric_EIV=len({r['source_species'] for r in eiv}),
            programmes_with_numeric_EIV=len({r['programme_id'] for r in eiv}),
            EIV_min=min(moisture) if moisture else '', EIV_max=max(moisture) if moisture else '',
            programmes_with_representative_coordinates=len({r['programme_id'] for r in members
                                                            if r['has_site_coordinate'] == 'True'}),
            programmes_with_any_monthly_water_context=len(programmes & water_programmes),
            programmes_with_at_least_one_valid_water_context=len(programmes & valid_water_programmes),
            moderator_interpretation='species_affinity_and_representative_landscape_context_not_plot_hydrology',
            moderator_pooling_status='not_established_by_coverage_counts'))
    with (ROOT / 'endpoint_coverage.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    audit = dict(input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 water_context_sha256=hashlib.sha256(water_source.read_bytes()).hexdigest(),
                 effect_rows=len(rows), assigned_rows=sum(r['effect_rows'] for r in output),
                 study_labels=len({r['study_id'] for r in rows}),
                 programme_labels=len({r['programme_id'] for r in rows}),
                 limitation='Provisional catalogue coverage, not final eligibility or independent evidence. '
                 'Stage programme counts overlap. Unextracted and nonfinite effects are outside this table. '
                 'Any valid water month is an acquisition flag, not a verified effect-specific exposure match. '
                 'Radius/month rows are not independent programmes. Zero detected open water is not dry soil.')
    assert audit['effect_rows'] == audit['assigned_rows']
    (ROOT / 'endpoint_coverage_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    for item in output:
        print(item['response_group'], item['effect_rows'], item['programme_labels'])


if __name__ == '__main__':
    main()
