"""Source-pinned four-species PLANTATT extraction; no habitat regression."""
from pathlib import Path
import csv, hashlib, json, datetime, argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
args = parser.parse_args()
ROOT = args.out
ROOT.mkdir(parents=True, exist_ok=True)
raw = args.source.read_bytes()
blob = hashlib.sha1(f'blob {len(raw)}\0'.encode() + raw).hexdigest()
assert blob == 'e6a1da929af81f370d2a1b79177b489e3e573a06'
source_rows = list(csv.DictReader(raw.decode('latin1').splitlines()))
assert len(source_rows) == 3227
names = {
    'Agrostis tenuis': ('Agrostis capillaris', 'authoritative_Kew_POWO_synonym'),
    'Anthoxanthum odoratum': ('Anthoxanthum odoratum', 'exact_binomial'),
    'Hypochaeris radicata': ('Hypochaeris radicata', 'exact_binomial'),
    'Lotus pedunculatus': ('Lotus pedunculatus', 'exact_binomial'),
}
rows = []
for original, (matched, method) in names.items():
    candidates = [(i + 2, r) for i, r in enumerate(source_rows) if r['taxon_name_binom'] == matched]
    assert len(candidates) == 1
    source_line, r = candidates[0]
    assert r['taxonomic_status'] == 'Accepted'
    assert r['aggregate'] == ''
    rows.append({
        'source_species': original,
        'matched_Kew_species': matched,
        'join_method': method,
        'F_PLANTATT': int(r['F_PLANTATT']),
        'L_PLANTATT': int(r['L_PLANTATT']),
        'taxonomic_status': r['taxonomic_status'],
        'aggregate_flag': r['aggregate'],
        'Kew_id': r['kew_id'],
        'source_csv_line': source_line,
        'POWO_taxon_url': r['POWO_URL'],
        'synonym_evidence_url': 'https://powo.science.kew.org/taxon/386415-1' if method != 'exact_binomial' else '',
        'synonym_evidence': 'Agrostis tenuis Sibth. is a synonym of Agrostis capillaris' if method != 'exact_binomial' else '',
    })
assert [(r['F_PLANTATT'],r['L_PLANTATT']) for r in rows] == [(5,6),(6,7),(4,8),(8,7)]
manifest = {
    'purpose': 'Verified PLANTATT habitat-affinity support for four complete, survival-reconciled Penryn plant species; not a general moisture moderator estimate',
    'retrieved_and_checked_date_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'source_url': 'https://raw.githubusercontent.com/RBGKew/BIFloraExplorer/master/data-raw/BI_MAIN.csv',
    'source_browser_url': 'https://github.com/RBGKew/BIFloraExplorer/blob/master/data-raw/BI_MAIN.csv',
    'source_dataset_publication_doi': '10.1038/s41597-021-01104-5',
    'source_sha256': hashlib.sha256(raw).hexdigest(),
    'source_git_blob_sha1': blob,
    'source_encoding_used': 'latin1; exact target names and numeric fields are ASCII',
    'source_rows': len(source_rows),
    'name_source': 'Kew POWO taxon386415-1 retrieved with web search/open 2026-11-21; asserts Agrostis tenuis Sibth. synonym of Agrostis capillaris. Original biomass table reports binomial without authority.',
    'rows': rows,
    'interpretation': {
        'F_is': 'British/Irish PLANTATT ordinal habitat-moisture affinity, separate from European EIV MOIST; not measured mesocosm soil water',
        'L_is': 'British/Irish PLANTATT ordinal habitat-light affinity; not artificial-night-light sensitivity or the experiment light treatment code',
        'high_F_species': ['Lotus pedunculatus'],
        'one_high_F_taxon_confounded_with_species_identity': True,
        'reference_mean_L_PLANTATT': sum(r['L_PLANTATT'] for r in rows[:3]) / 3,
        'Lotus_L_PLANTATT': rows[3]['L_PLANTATT'],
        'same_L_pair': ['Lotus pedunculatus', 'Anthoxanthum odoratum'],
        'not_identified': ['causal effect of moisture affinity', 'response to soil moisture manipulation', 'mediation by herbivores or predators', 'general effect across wet plant species'],
    },
}
(ROOT / 'penryn_four_species_PLANTATT_verified.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
with (ROOT / 'penryn_four_species_PLANTATT_verified.tsv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t',
                       lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
    w.writeheader()
    w.writerows(rows)
print(json.dumps({'source_git_blob_sha1': blob, 'source_sha256': manifest['source_sha256'], 'rows': rows, 'interpretation':manifest['interpretation']}, indent=2))
