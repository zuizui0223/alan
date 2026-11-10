#!/usr/bin/env python3
"""Validate a conservative moth-like flower phenotype against positive GloBI evidence.

The source taxon-by-trait input is access-controlled. Output contains only
aggregate diagnostics and the public colour/form profile pre-filter.
GloBI flower visits are not demonstrated effective nocturnal pollination.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

MOTH_FORMS = frozenset({'tubular', 'salverform', 'funnel_trumpet'})
GUILDS = {'moth_family_visit': 'moth', 'butterfly_family_visit': 'butterfly',
          'moth': 'moth', 'butterfly': 'butterfly'}
RELATIONS = frozenset({'visitsFlowersOf', 'pollinates'})


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def read_csv(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames, f'Empty/invalid CSV: {path}')
        yield from reader


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def states(value):
    """One recorded state set, never a Cartesian product."""
    value = (value or '').strip()
    if not value or value == 'unresolved':
        return None
    parts = value.split('|')
    if any(not p or p != p.strip() or p == 'unresolved' for p in parts):
        return None
    return frozenset(parts) if len(parts) == len(set(parts)) else None


def reference_profiles(path):
    result, seen = {}, set()
    for row in read_csv(path):
        key = row['profile_id']
        pair = (states(row['flower_primary_color']), states(row['floral_form']))
        require(pair[0] and pair[1], f'Invalid public fd profile: {key}')
        require(key not in result and pair not in seen, 'Duplicate fd reference profile')
        result[key] = pair
        seen.add(pair)
    require(result, 'No fd reference profiles')
    return result


def phenotype(color, form, depth, reference):
    """Return candidate/noncandidate/unknown; missing is NOT noncandidate."""
    c, f, d = states(color), states(form), states(depth)
    if not c or not f or not d or len(c) != 1 or len(f) != 1 or len(d) != 1:
        return 'unknown'
    if (c, f) not in set(reference.values()):
        return 'unknown'
    return ('candidate' if c == {'white'} and f <= MOTH_FORMS and d == {'deep'}
            else 'noncandidate')


def read_traits(path, reference):
    result = {}
    for row in read_csv(path):
        name = row['accepted_species'].strip()
        require(name and name not in result, f'Blank/duplicate accepted taxon: {name}')
        result[name] = phenotype(row['flower_primary_color'], row['floral_form'],
                                 row['tube_depth_class'], reference)
    return result


def read_globi(path):
    positive = defaultdict(set)
    for row in read_csv(path):
        rel = (row.get('relation') or row.get('interaction_type') or '').strip()
        if rel not in RELATIONS:
            continue
        guild = GUILDS.get((row.get('guild') or row.get('visitor_guild') or '').strip())
        name = (row.get('plant_name') or '').strip()
        if guild and name:
            positive[name].add(guild)
    return positive


def summarize_validation(traits, interactions):
    bins = Counter()
    for name, phenotype_label in traits.items():
        visitor = interactions.get(name, set())
        label = ('both' if visitor == {'moth', 'butterfly'} else
                 'moth_recorded_only' if visitor == {'moth'} else
                 'butterfly_recorded_only' if visitor == {'butterfly'} else 'no_guild_record')
        bins[(label, phenotype_label)] += 1
    result = {
        'schema': 'alan_moth_shape_GloBI_positive_evidence_v1',
        'interpretation': 'Recorded moth/butterfly flower-visits, NOT true pollination absence or moth dependence',
        'total_unique_trait_taxa': len(traits),
        'counts': {k: {v: bins[(k, v)] for v in ('candidate', 'noncandidate', 'unknown')}
                   for k in ('moth_recorded_only', 'butterfly_recorded_only',
                             'both', 'no_guild_record')},
    }
    m, b = result['counts']['moth_recorded_only'], result['counts']['butterfly_recorded_only']
    a, c = m['candidate'], m['noncandidate']
    d, e = b['candidate'], b['noncandidate']
    # Descriptive comparison on disjoint, interaction-observed strata only.
    result['recorded_moth_vs_butterfly_only_candidate_odds_ratio_haldane'] = (
        (a + 0.5) * (e + 0.5) / ((c + 0.5) * (d + 0.5))
        if (a + c) and (d + e) else None
    )
    result['biological_specificity_estimable'] = False
    result['GloBI_missing_as_absence'] = False
    result['multi_state_combinations_invented'] = False
    return result


def run(profiles, traits, interactions, output):
    profiles, traits, interactions, output = map(Path, (profiles, traits, interactions, output))
    require(not output.exists(), 'Output directory exists: refuse accidental overwrite')
    ref = reference_profiles(profiles)
    labels = read_traits(traits, ref)
    positives = read_globi(interactions)
    report = summarize_validation(labels, positives)
    report['n_reference_profiles'] = len(ref)
    report['source_sha256'] = {p.name: sha256(p) for p in (profiles, traits, interactions)}
    report['sources_are_positive_evidence_only'] = True
    output.mkdir(parents=True)
    with (output / 'public_fd_colour_form_prefilter.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['profile_id', 'white_tubular_prefilter'])
        writer.writeheader()
        for key, (c, f) in sorted(ref.items()):
            writer.writerow({'profile_id': key, 'white_tubular_prefilter':
                             int(c == {'white'} and len(f) == 1 and f <= MOTH_FORMS)})
    (output / 'validation_aggregate.json').write_text(
        json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + '\n', encoding='utf-8')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ('profiles', 'traits', 'interactions', 'output'):
        p.add_argument('--' + arg, type=Path, required=True)
    args = p.parse_args()
    result = run(args.profiles, args.traits, args.interactions, args.output)
    print(json.dumps({k: result[k] for k in ('n_reference_profiles',
            'total_unique_trait_taxa', 'counts')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
