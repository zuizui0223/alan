"""Reconcile search records, reviewed reports and extractions without double counting.

Coverage is not final eligibility. Earlier incremental and batch reviews can
describe the same paper; extraction status never resolves a scope hold.
"""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parent
S = ROOT / 'search'
identity = pd.read_csv(S / 'record_identity_crosswalk.csv').fillna('')
batch = pd.read_csv(S / 'fulltext_design_reviews.csv').fillna('')
incremental = pd.read_csv(S / 'fulltext_decisions_incremental.csv').fillna('')
studies = pd.read_csv(ROOT / 'studies.csv').fillna('')
catalogue = pd.read_csv(ROOT / 'effect_catalogue.csv')
assert identity.record_id.is_unique and studies.study_id.is_unique
assert set(catalogue.study_id) <= set(studies.study_id)
assert studies.doi.is_unique
evidence = {}
def add(doi, source, read=False, extracted=False, decision=''):
    key = 'doi:' + doi.strip().lower()
    x = evidence.setdefault(key, dict(sources=[], decisions=[], read=False, extracted=False))
    x['sources'].append(source)
    if decision:
        x['decisions'].append(decision)
    x['read'] |= read
    x['extracted'] |= extracted

for x in batch.itertuples():
    add(x.doi, x.review_file, x.fulltext_access_status == 'read_primary_methods_results',
        decision=x.decision)
for x in incremental.itertuples():
    add(x.doi, x.source_note, x.review_status in
        ['methods_results_and_original_data_read', 'methods_and_results_read'], decision=x.decision)
for x in studies.itertuples():
    extracted = x.study_id in set(catalogue.study_id)
    # Registry-only designs may have partial Methods reading; retain them as
    # leads unless a full review or actual catalogue extraction establishes reading.
    add(x.doi, 'studies.csv:' + x.study_id, read=extracted, extracted=extracted,
        decision=x.status)

rows = []
for key, g in identity.groupby('report_key', sort=False):
    e = evidence.get(key, dict(sources=[], decisions=[], read=False, extracted=False))
    needed = g.decision.eq('fulltext_needed').any()
    if e['read']:
        coverage = 'primary_read_eligibility_not_final'
    elif e['sources']:
        coverage = 'lead_or_access_hold_requires_full_review'
    else:
        coverage = 'no_documented_primary_review'
    rows.append(dict(report_key=key, record_ids=';'.join(g.record_id),
                     search_records=len(g), doi=g.resolved_doi.iloc[0],
                     pmcid=';'.join(sorted(set(g.pmcid) - {''})), title=g.title.iloc[0],
                     search_decisions=';'.join(sorted(set(g.decision))),
                     retained_for_fulltext=bool(needed), coverage=coverage,
                     catalogue_effects_present=e['extracted'],
                     evidence_sources=';'.join(sorted(set(e['sources']))),
                     decision_history=';'.join(sorted(set(e['decisions'])))))
out = pd.DataFrame(rows)
assert out.report_key.is_unique and out.search_records.sum() == len(identity)
out.to_csv(S / 'report_review_coverage.csv', index=False)
pending = out[out.retained_for_fulltext & out.coverage.ne('primary_read_eligibility_not_final')]
pending.to_csv(S / 'pending_primary_review.csv', index=False)
outside = [{'report_key': k, 'evidence_sources': ';'.join(v['sources']),
            'primary_read': v['read'], 'catalogue_effects_present': v['extracted']}
           for k, v in evidence.items() if k not in set(out.report_key)]
pd.DataFrame(outside, columns=['report_key','evidence_sources','primary_read',
                             'catalogue_effects_present']).to_csv(S / 'supplementary_report_sources.csv', index=False)
summary = dict(search_records=len(identity), unique_search_reports=len(out),
               retained_fulltext_reports=int(out.retained_for_fulltext.sum()),
               retained_reports_with_primary_read=int((out.retained_for_fulltext &
                   out.coverage.eq('primary_read_eligibility_not_final')).sum()),
               pending_primary_review_reports=len(pending),
               additional_source_reports_outside_search=len(outside),
               duplicate_incremental_batch_dois=sorted(set(batch.doi) & set(incremental.doi)),
               final_prisma=False, final_eligibility=False,
               note='Exact report identity only; programmes may share plots. Registry-only leads do not establish fulltext reading. No effect count is inferred from review counts.')
(S / 'report_review_coverage_summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
print(json.dumps(summary, indent=2))
