"""Exact-identity inventory of two frozen searches, not final PRISMA/adoption."""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parent
SEARCH = ROOT/'search'


def main():
    old_path = SEARCH/'record_identity_crosswalk.csv'
    extra_path = SEARCH/'pollination_expansion_screening.csv'
    old = pd.read_csv(old_path, keep_default_na=False)
    extra = pd.read_csv(extra_path, keep_default_na=False)
    assert len(old) == 521 and len(extra) == 51
    assert old.record_id.is_unique and extra.record_id.is_unique
    lookup = old.set_index('record_id')
    appearances = []
    for r in old.to_dict('records'):
        appearances.append(dict(search_id='europepmc_supplement_20261120',
            record_id=r['record_id'], doi=r['resolved_doi'].lower(), pmcid=r['pmcid'],
            report_key=r['report_key'], title=r['title'],
            screening_decision=r['decision'], identity_basis='original_resolved_crosswalk'))
    for r in extra.to_dict('records'):
        rid = r['record_id']
        if rid in lookup.index:
            prior = lookup.loc[rid]
            if r['doi'] and prior.resolved_doi:
                assert r['doi'].lower() == prior.resolved_doi.lower(), rid
            key, doi, basis = prior.report_key, prior.resolved_doi.lower(), 'same_record_id_preserve_original_resolution'
        else:
            doi = r['doi'].lower()
            key, basis = ('doi:'+doi if doi else 'record:'+rid), 'exact_doi_or_unresolved_record_id'
        appearances.append(dict(search_id='europepmc_pollination_expansion_20261120',
            record_id=rid, doi=doi, pmcid=r['pmcid'], report_key=key,
            title=r['title'], screening_decision=r['screening_decision'], identity_basis=basis))
    appearances = pd.DataFrame(appearances)
    assert appearances.groupby('record_id').report_key.nunique().eq(1).all()
    appearances.to_csv(SEARCH/'combined_search_appearances.csv', index=False)
    records = []
    for rid, g in appearances.groupby('record_id', sort=True):
        first = g.iloc[0]
        records.append(dict(record_id=rid, report_key=first.report_key, doi=first.doi,
            title=first.title, search_ids=';'.join(g.search_id), appearances=len(g),
            screening_history=';'.join(g.search_id+':'+g.screening_decision),
            final_eligibility='not_adjudicated_by_this_inventory'))
    records = pd.DataFrame(records)
    records.to_csv(SEARCH/'combined_search_records.csv', index=False)
    reports = records.groupby('report_key', sort=True).agg(
        record_ids=('record_id', lambda x: ';'.join(x)),
        unique_records=('record_id', 'size'),
        search_appearances=('appearances', 'sum')).reset_index()
    reports.to_csv(SEARCH/'combined_search_reports.csv', index=False)
    unresolved = pd.read_csv(SEARCH/'pollination_expansion_duplicate_candidates.csv', keep_default_na=False)
    assert set(unresolved.record1)|set(unresolved.record2) <= set(records.record_id)
    summary = dict(search_appearances=len(appearances), original_records=len(old),
        supplementary_records=len(extra), exact_record_overlap=len(set(old.record_id)&set(extra.record_id)),
        unique_records=len(records), exact_identity_reports=len(reports),
        original_reports=old.report_key.nunique(),
        new_report_keys_relative_to_original=len(set(reports.report_key)-set(old.report_key)),
        unresolved_candidate_pairs=len(unresolved), candidate_pairs_auto_merged=0,
        publication_families_reestimated=False, final_prisma=False, final_eligibility=False,
        independent_experiments_estimated=False,
        input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [old_path, extra_path]},
        limits=['Only the two named fixed EuropePMC searches; not all citation-tracked reports',
                'Distinct preprint/final DOI reports retained; no claim of family independence',
                'Missing DOI records retained by ID; title similarity not used for merging',
                'Original and supplemental screening histories preserved, not adjudicated',
                'Neither WoS nor Scopus search completed'])
    (SEARCH/'combined_search_summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
