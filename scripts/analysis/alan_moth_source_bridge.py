#!/usr/bin/env python3
"""Source-verified ALAN/fd/island trait feasibility, aggregate output only.

Inputs: authorized GitHub Actions artifact ZIPs. Source rights prohibit release
of the underlying species-level trait/status tables. The historical island
snapshot is from September 2026, NOT a live overwrite of Chapter 1.
"""
import argparse
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pandas as pd

TRAIT_FILES = {
    'all': ('snapshot/chapter1_trait_ledger_all_analysis.csv.gz',
            '2d95136fc7df1dcce92aab738ce0364cdb1434746e1f7b794c529900192cbf4d'),
    'direct': ('snapshot/chapter1_trait_ledger_direct_only.csv.gz',
               'd7d6a0372ce7f3ebe958047a222ae0504f01564e364f93bab360273235c6661b'),
}
STATUS_FILE = ('fixed/canonical/input/chapter1_status_flora.csv.gz',
               '6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a')
GLOBI_FILE = ('GloBI_stage2_signed_extension.csv.gz',
              'a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994')
PHOTO_FILE = ('island_year_ROI_free_source_support.csv',
              'acaf9287a29ef01e2fc79fd4b3450b4c09bb5eee16737e55b4326a80359cbdf8')
RELEVANT = ('flower_primary_color', 'floral_form', 'tube_depth_class')
FORMS = frozenset({'tubular', 'salverform', 'funnel_trumpet'})


def require(value, message):
    if not value:
        raise ValueError(message)


def source(z, contract):
    path, expected = contract
    content = z.read(path)
    require(hashlib.sha256(content).hexdigest() == expected, 'Source digest mismatch: ' + path)
    buf = io.BytesIO(content)
    return pd.read_csv(gzip.GzipFile(fileobj=buf) if path.endswith('.gz') else buf, dtype=str)


def read_species_traits(frame):
    require({'accepted_species', 'trait_name', 'normalized_value'} <= set(frame),
            'Trait schema mismatch')
    d = frame[frame.trait_name.isin(RELEVANT)].copy()
    require(not d.duplicated(['accepted_species', 'trait_name']).any(),
            'Duplicated species-trait state')
    w = d.pivot(index='accepted_species', columns='trait_name', values='normalized_value')
    require(set(RELEVANT) <= set(w), 'Incomplete trait schema')
    return w


def positive_guilds(interactions):
    require({'plant_name', 'guild', 'relation'} <= set(interactions),
            'GloBI schema mismatch')
    require(interactions.guild.isin(['moth_family_visit', 'butterfly_family_visit']).all(),
            'Unrecognized visitor family')
    require(interactions.relation.isin(['visitsFlowersOf', 'pollinates']).all(),
            'Unrecognized interaction')
    links = interactions[['plant_name', 'guild']].drop_duplicates()
    flags = links.groupby('plant_name').guild.agg(lambda v: frozenset(v)).to_dict()

    def group(s):
        if s == {'moth_family_visit', 'butterfly_family_visit'}:
            return 'both'
        if s == {'moth_family_visit'}:
            return 'moth_recorded_only'
        if s == {'butterfly_family_visit'}:
            return 'butterfly_recorded_only'
        return 'no_guild_record'

    return {k: group(v) for k, v in flags.items()}


def report_scope(frame, visitors, status, photo_ids):
    traits = read_species_traits(frame)
    pair = traits.flower_primary_color.notna() & traits.floral_form.notna()
    broad = traits.flower_primary_color.eq('white') & traits.floral_form.isin(FORMS)
    depth_known = traits.tube_depth_class.notna()
    strict = broad & traits.tube_depth_class.eq('deep')
    complete = pair & depth_known

    def guild_breakdown(names):
        labels = pd.Series(names, dtype=str).map(visitors).fillna('no_guild_record')
        return {g: int(labels.eq(g).sum()) for g in
                ('moth_recorded_only', 'butterfly_recorded_only', 'both', 'no_guild_record')}

    def flora_breakdown(names):
        x = status.loc[status.accepted_species.isin(names)]
        y = x.loc[x.origin_status.eq('introduced')]
        return {
            'observed_island_taxon_pairs': int(len(x)),
            'distinct_islands': int(x.island_id.nunique()),
            'status_pairs': {s: int(x.origin_status.eq(s).sum()) for s in
                             ('introduced', 'native', 'unresolved')},
            'introduced_pairs_NOT_naturalization': int(len(y)),
            'introduced_species': int(y.accepted_species.nunique()),
            'islands_with_introduced_record': int(y.island_id.nunique()),
            'provisional_id_overlap_with_photo_islands': int(len(set(x.island_id) & photo_ids)),
            'introduced_islands_provisional_photo_id_overlap': int(len(set(y.island_id) & photo_ids)),
        }

    broad_names = traits.index[broad]
    strict_names = traits.index[strict]
    moth_positive_names = [name for name in broad_names
                           if visitors.get(name) in ('moth_recorded_only', 'both')]
    return {
        'names_with_one_or_more_relevant_traits': int(len(traits)),
        'names_with_colour_and_form': int(pair.sum()),
        'names_with_all_three_traits': int(complete.sum()),
        'white_single_form_candidate_names': int(broad.sum()),
        'white_single_form_depth_unknown': int((broad & ~depth_known).sum()),
        'white_single_form_depth_observed': int((broad & depth_known).sum()),
        'white_single_form_depth_deep': int(strict.sum()),
        'white_form_guild_positive_evidence': guild_breakdown(broad_names),
        'strict_deep_guild_positive_evidence': guild_breakdown(strict_names),
        'white_form_positive_moth_visitor_names': len(moth_positive_names),
        'white_form_moth_recorded_flora': flora_breakdown(moth_positive_names),
        'white_form_flora': flora_breakdown(broad_names),
        'strict_deep_flora': flora_breakdown(strict_names),
    }


def run(island_artifact, globi_artifact, photo_artifact, output):
    output = Path(output)
    require(not output.exists(), 'Refuse overwrite of signed aggregate')
    with zipfile.ZipFile(island_artifact) as island, zipfile.ZipFile(globi_artifact) as g, zipfile.ZipFile(photo_artifact) as p:
        status = source(island, STATUS_FILE)
        interactions = source(g, GLOBI_FILE)
        photos = source(p, PHOTO_FILE)
        require(not status.duplicated(['island_id', 'accepted_species']).any(),
                'Ambiguous island status names')
        require(status.origin_status.isin(['introduced', 'native', 'unresolved']).all(),
                'Unknown historical flora status code')
        require(status.status_conflict.eq('False').all(), 'Floristic status conflicts')
        require(not photos.duplicated(['island_id', 'year']).any(), 'Duplicate photo island year')
        require(len(interactions) == 61967, 'GloBI claim count drift')
        visitors = positive_guilds(interactions)
        photo_islands = set(photos.island_id)
        scopes = {key: report_scope(source(island, contract), visitors, status, photo_islands)
                  for key, contract in TRAIT_FILES.items()}

    result = {
        'schema': 'alan_moth_fd_island_authorized_snapshot_source_bridge_20261111_v1',
        'status': 'ACTUAL source join / feasibility, NOT gecko-naturalization-ALAN effect',
        'source_run_ids': {'island': 34232450884, 'globi': 41976552977,
                           'photo_source': 42171173575},
        'n_signed_GloBI_flower_interaction_rows': 61967,
        'n_GloBI_unique_plant_names': len(visitors),
        'n_photo_source_islands': len(photo_islands),
        'n_photo_source_island_years': len(photos),
        'n_historical_island_flora_pairs': len(status),
        'floristic_status_is_NOT_naturalization': True,
        'gecko_and_annual_VIIRS_not_tested': True,
        'photo_island_ID_overlap_NOT_polygon_equivalence': True,
        'unlisted_GloBI_not_negative_evidence': True,
        'no_restricted_identity_rows_exported': True,
        'trait_coverage_scopes': scopes,
        'verified_members_sha256': dict([*TRAIT_FILES.values(), STATUS_FILE, GLOBI_FILE, PHOTO_FILE]),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + '\n',
                      encoding='utf-8')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('island_artifact', 'globi_artifact', 'photo_artifact', 'output'):
        p.add_argument('--' + key.replace('_', '-'), type=Path, required=True)
    a = p.parse_args()
    result = run(a.island_artifact, a.globi_artifact, a.photo_artifact, a.output)
    print(json.dumps({k: result[k] for k in ('n_photo_source_islands',
           'n_GloBI_unique_plant_names', 'trait_coverage_scopes')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
