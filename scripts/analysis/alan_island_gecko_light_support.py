#!/usr/bin/env python3
"""Audit island support for moth-like naturalized floras × ALAN × geckos.

Requires independently reviewed ISLAND-specific establishment and gecko data.
GRIIS country-level lists and missing occurrences cannot establish island absence.
The island table is a PRIVATE derived work; do not commit it to public repos.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from alan_moth_syndrome_validate import read_csv, read_traits, reference_profiles, require, sha256


STATUS = frozenset({'naturalized', 'native', 'introduced', 'cultivated', 'transient', 'unknown'})
GECKO = frozenset({'present', 'confirmed_absent', 'unknown'})


def unique_index(path, key):
    out = {}
    for row in read_csv(path):
        val = row[key].strip()
        require(val and val not in out, f'Blank/repeated {key}: {val}')
        out[val] = row
    return out


def load_flora(path):
    flora = defaultdict(dict)
    for row in read_csv(path):
        if row['review_status'] != 'accepted':
            continue
        island, name = row['island_id'].strip(), row['accepted_species'].strip()
        status = row['establishment_status'].strip()
        require(island and name and status in STATUS, f'Invalid island-specific status: {row}')
        require(name not in flora[island], f'Duplicate reviewed island/taxon: {island}/{name}')
        year = row.get('establishment_year', '').strip()
        if year:
            require(year.isdigit() and 1500 <= int(year) <= 2026,
                    f'Invalid documented establishment year: {island}/{name}')
        flora[island][name] = (status, int(year) if year else None)
    return flora


def load_geckos(path):
    out = {}
    for row in read_csv(path):
        if row['review_status'] != 'accepted':
            continue
        island = row['island_id'].strip()
        status = row['gecko_status'].strip()
        require(row['geographic_scope'] == 'island',
                f'Cannot infer island presence from country-only checklist: {island}')
        require(island and status in GECKO and row['source'].strip(),
                f'Invalid reviewed island gecko status: {island}')
        require(island not in out, f'Conflicting gecko evidence for island: {island}')
        out[island] = status
    return out


def load_light(path, *, first_year=2019, last_year=2025,
               min_months=6, min_coverage=0.7):
    light = defaultdict(list)
    seen = set()
    for row in read_csv(path):
        island = row['island_id'].strip()
        require(island and row['year'].isdigit(), 'Invalid light island/year')
        year = int(row['year'])
        key = (island, year)
        require(key not in seen, f'Duplicate VIIRS island/year: {key}')
        seen.add(key)
        require(row['product_version'].strip(), f'Unspecified VIIRS product: {key}')
        if not first_year <= year <= last_year:
            continue
        try:
            rad = float(row['radiance'])
            months = int(row['n_valid_months'])
            coverage = float(row['coverage_fraction'])
        except (ValueError, TypeError) as exc:
            raise ValueError(f'Invalid VIIRS radiance/coverage: {key}') from exc
        require(0 <= coverage <= 1 and 0 <= months <= 12,
                f'Invalid VIIRS quality metadata: {key}')
        if math.isfinite(rad) and rad >= 0 and months >= min_months and coverage >= min_coverage:
            light[island].append(rad)
    return light


def island_counts(islands, flora, labels, geckos, light, *, min_light_years=2):
    result = []
    for island, meta in sorted(islands.items()):
        area = float(meta['area_km2'])
        require(math.isfinite(area) and area > 0, f'Invalid island area: {island}')
        archipelago = meta['archipelago'].strip()
        require(archipelago, f'Unknown archipelago label: {island}')
        counts = Counter()
        for name, (status, established) in flora.get(island, {}).items():
            if status not in {'naturalized', 'native'}:
                counts['unqualified_status'] += 1
                continue
            kind = labels.get(name, 'unknown')
            counts[f'{status}_{kind}'] += 1
            if status == 'naturalized' and established is not None and established >= 2013:
                counts['dated_post2012_naturalized'] += 1
                if kind != 'unknown':
                    counts['dated_post2012_naturalized_known_trait'] += 1
        naturalized_known = counts['naturalized_candidate'] + counts['naturalized_noncandidate']
        records = light.get(island, [])
        result.append({
            'island_id': island, 'archipelago': archipelago, 'area_km2': area,
            'gecko_status': geckos.get(island, 'unknown'),
            'n_valid_light_years': len(records),
            'mean_annual_radiance': (sum(records) / len(records)
                                     if len(records) >= min_light_years else ''),
            'n_naturalized_mothlike': counts['naturalized_candidate'],
            'n_naturalized_other_known': counts['naturalized_noncandidate'],
            'n_naturalized_trait_unknown': counts['naturalized_unknown'],
            'n_naturalized_known': naturalized_known,
            'n_native_mothlike': counts['native_candidate'],
            'n_native_other_known': counts['native_noncandidate'],
            'n_native_trait_unknown': counts['native_unknown'],
            'n_dated_post2012_naturalized': counts['dated_post2012_naturalized'],
            'n_dated_post2012_naturalized_known_trait': counts['dated_post2012_naturalized_known_trait'],
        })
    return result


def summarize(rows, flora, geckos, light, *, min_naturalized=10, min_light_years=2):
    has_flora = {i for i, records in flora.items() if records}
    eligible = [r for r in rows if r['n_naturalized_known'] >= min_naturalized and
                r['gecko_status'] in {'present', 'confirmed_absent'} and
                r['n_valid_light_years'] >= min_light_years]
    light_threshold = statistics.median(r['mean_annual_radiance'] for r in eligible) if eligible else None
    strata = {}
    for gecko in ('present', 'confirmed_absent'):
        for light_bin in ('low_or_median', 'high'):
            subset = [r for r in eligible if r['gecko_status'] == gecko and
                      ('high' if r['mean_annual_radiance'] > light_threshold
                       else 'low_or_median') == light_bin]
            strata[gecko + '__' + light_bin] = {
                'islands': len(subset),
                'naturalized_mothlike': sum(r['n_naturalized_mothlike'] for r in subset),
                'naturalized_other_known': sum(r['n_naturalized_other_known'] for r in subset),
            }
    archipelagos = defaultdict(set)
    for row in eligible:
        archipelagos[row['archipelago']].add(row['gecko_status'])
    return {
        'schema': 'alan_island_naturalization_gecko_light_support_v1',
        'islands_in_reference_metadata': len(rows),
        'islands_with_any_reviewed_flora': len(has_flora),
        'islands_with_reviewed_gecko_status': sum(1 for r in rows if r['gecko_status'] != 'unknown'),
        'islands_with_sufficient_annual_light': sum(1 for r in rows if r['n_valid_light_years'] >= min_light_years),
        'islands_with_at_least_one_known_naturalized_trait': sum(r['n_naturalized_known'] > 0 for r in rows),
        'islands_with_all_three_sources_and_naturalized_threshold': len(eligible),
        'eligible_archipelagos': len(archipelagos),
        'archipelagos_with_both_verified_gecko_statuses': sum(len(x) == 2 for x in archipelagos.values()),
        'eligible_islands_with_dated_post2012_naturalized_records': sum(
            r['n_dated_post2012_naturalized'] > 0 for r in eligible),
        'eligible_islands_with_dated_post2012_known_trait': sum(
            r['n_dated_post2012_naturalized_known_trait'] > 0 for r in eligible),
        'diagnostic_median_light_radiance': light_threshold,
        'descriptive_2x2': strata,
        'eligibility': {'min_naturalized_known': min_naturalized,
                        'min_annual_light_years': min_light_years},
        'claim': 'source support ONLY, not an ALAN effect or naturalization causal test',
        'no_gecko_record_is_not_absence': True,
        'introduced_without_naturalization_not_counted': True,
        'native_vs_naturalized_status_is_island_specific': True,
        'dated_post2012_records_are_not_proof_of_postlight_recruitment': True,
    }


def run(profiles, traits, floras, geckos, light, islands, output,
        *, min_naturalized=10, min_light_years=2, first_year=2019, last_year=2025):
    paths = [Path(p) for p in (profiles, traits, floras, geckos, light, islands)]
    output = Path(output)
    require(not output.exists(), 'Output directory already exists: refusing overwrite')
    ref = reference_profiles(paths[0])
    labels = read_traits(paths[1], ref)
    floras_by_island = load_flora(paths[2])
    gecko_by_island = load_geckos(paths[3])
    annual_light = load_light(paths[4], first_year=first_year, last_year=last_year)
    islands_by_id = unique_index(paths[5], 'island_id')
    require(set(floras_by_island) <= set(islands_by_id), 'Flora island IDs not in pinned metadata')
    require(set(gecko_by_island) <= set(islands_by_id), 'Gecko island IDs not in pinned metadata')
    require(set(annual_light) <= set(islands_by_id), 'Light island IDs not in pinned metadata')
    rows = island_counts(islands_by_id, floras_by_island, labels, gecko_by_island,
                         annual_light, min_light_years=min_light_years)
    report = summarize(rows, floras_by_island, gecko_by_island, annual_light,
                       min_naturalized=min_naturalized, min_light_years=min_light_years)
    report['light_year_window'] = [first_year, last_year]
    report['source_sha256'] = {p.name: sha256(p) for p in paths}
    output.mkdir(parents=True)
    with (output / 'PRIVATE_island_derived_counts_DO_NOT_COMMIT.csv').open(
            'w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    (output / 'support_aggregate.json').write_text(
        json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + '\n', encoding='utf-8')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('profiles', 'traits', 'floras', 'geckos', 'light', 'islands', 'output'):
        p.add_argument('--' + key, type=Path, required=True)
    p.add_argument('--min-naturalized', type=int, default=10)
    p.add_argument('--min-light-years', type=int, default=2)
    p.add_argument('--first-year', type=int, default=2019)
    p.add_argument('--last-year', type=int, default=2025)
    a = p.parse_args()
    require(2012 <= a.first_year <= a.last_year <= 2026, 'Invalid VIIRS year window')
    require(a.min_naturalized >= 1 and a.min_light_years >= 1, 'Invalid support thresholds')
    result = run(a.profiles, a.traits, a.floras, a.geckos, a.light, a.islands,
                 a.output, min_naturalized=a.min_naturalized,
                 min_light_years=a.min_light_years, first_year=a.first_year,
                 last_year=a.last_year)
    print(json.dumps({k: result[k] for k in (
        'islands_in_reference_metadata', 'islands_with_all_three_sources_and_naturalized_threshold',
        'eligible_archipelagos', 'descriptive_2x2')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
