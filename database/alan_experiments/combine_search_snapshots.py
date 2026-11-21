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
    resolution_path = SEARCH/'combined_search_identity_resolutions.csv'
    resolutions = pd.read_csv(resolution_path, keep_default_na=False).set_index('record_id')
    assert resolutions.index.is_unique
    assert set(resolutions.index) <= set(extra.record_id)-set(old.record_id)
    for row in resolutions.itertuples():
        if row.matched_record in lookup.index:
            target_doi = lookup.loc[row.matched_record, 'resolved_doi']
        else:
            target_doi = extra.set_index('record_id').loc[row.matched_record, 'doi']
        assert target_doi.lower() == row.resolved_doi
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
            if rid in resolutions.index:
                resolved = resolutions.loc[rid]
                assert not doi or doi == resolved.resolved_doi
                doi = resolved.resolved_doi
                key, basis = 'doi:'+doi, 'manual_primary_bibliography_resolution'
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
    versions = pd.concat([pd.read_csv(SEARCH/'publication_version_links.csv'),
                          pd.read_csv(SEARCH/'combined_search_additional_version_links.csv')], ignore_index=True)
    assert versions.earlier_doi.is_unique
    assert not set(versions.earlier_doi)&set(versions.final_doi)
    assert {'doi:'+x for x in list(versions.earlier_doi)+list(versions.final_doi)} <= set(reports.report_key)
    family_map = {'doi:'+r.earlier_doi:'doi:'+r.final_doi for r in versions.itertuples()}
    reports['publication_family_key'] = reports.report_key.map(lambda k: family_map.get(k,k))
    reports.to_csv(SEARCH/'combined_search_reports.csv', index=False)
    unresolved = pd.read_csv(SEARCH/'pollination_expansion_duplicate_candidates.csv', keep_default_na=False)
    assert set(unresolved.record1)|set(unresolved.record2) <= set(records.record_id)
    id_report = records.set_index('record_id').report_key.to_dict()
    candidate_audit=[]
    for row in unresolved.itertuples():
        a,b=id_report[row.record1],id_report[row.record2]
        status=('resolved_same_report_primary_metadata' if a==b else
                'resolved_version_link_primary_metadata' if family_map.get(a,a)==family_map.get(b,b) else
                'unresolved')
        candidate_audit.append(dict(record1=row.record1,record2=row.record2,status=status))
    pd.DataFrame(candidate_audit).to_csv(SEARCH/'combined_search_candidate_adjudications.csv',index=False)
    summary = dict(search_appearances=len(appearances), original_records=len(old),
        supplementary_records=len(extra), exact_record_overlap=len(set(old.record_id)&set(extra.record_id)),
        unique_records=len(records), exact_identity_reports=len(reports),
        original_reports=old.report_key.nunique(),
        new_report_keys_relative_to_original=len(set(reports.report_key)-set(old.report_key)),
        original_candidate_pairs=len(unresolved),
        unresolved_candidate_pairs=sum(r['status']=='unresolved' for r in candidate_audit),
        candidate_pairs_auto_merged=0, manual_record_identity_resolutions=len(resolutions),
        publication_families_after_verified_links=reports.publication_family_key.nunique(),
        verified_version_pairs=len(versions),
        publication_families_reestimated=True, final_prisma=False, final_eligibility=False,
        independent_experiments_estimated=False,
        input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [old_path, extra_path,
            resolution_path, SEARCH/'publication_version_links.csv', SEARCH/'combined_search_additional_version_links.csv']},
        limits=['Only the two named fixed EuropePMC searches; not all citation-tracked reports',
                'Distinct preprint/final DOI reports retained; no claim of family independence',
                'Missing DOI records retained by ID; title similarity not used for merging',
                'Original and supplemental screening histories preserved, not adjudicated',
                'Neither WoS nor Scopus search completed'])
    (SEARCH/'combined_search_summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
