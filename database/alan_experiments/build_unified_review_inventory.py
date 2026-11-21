"""Join review evidence across search and citation routes, without final eligibility."""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parent/'search'


def main():
    names = ['combined_search_reports.csv', 'combined_search_appearances.csv',
             'report_review_coverage.csv', 'supplementary_report_sources.csv',
             'pollination_expansion_fulltext_followup.csv', 'partial_report_reviews.csv']
    reports, appearances, original, supplementary, followup, partial_reviews = [
        pd.read_csv(ROOT/n, keep_default_na=False) for n in names]
    original = original.set_index('report_key')
    supplementary = supplementary.set_index('report_key')
    report_lookup = reports.set_index('report_key')
    fulltext_statuses = {'fulltext_setting_needed', 'spectrum_comparator_fulltext_needed',
                        'visual_model_fulltext_needed'}
    rows = []
    for key in sorted(set(reports.report_key)|set(supplementary.index)):
        in_search, in_citation = key in report_lookup.index, key in supplementary.index
        read, pending_original, pending_supplement = False, False, False
        evidence, history, partial = [], [], []
        if key in original.index:
            r = original.loc[key]
            read = r.coverage == 'primary_read_eligibility_not_final'
            pending_original = str(r.retained_for_fulltext).lower() == 'true' and not read
            evidence.append(r.evidence_sources)
            history.append(r.decision_history)
        if in_citation:
            r = supplementary.loc[key]
            source_read = str(r.primary_read).lower() == 'true'
            read = read or source_read
            pending_supplement = not source_read
            evidence.append(r.evidence_sources)
            history.append(r.decision_history)
            partial.append(r.partial_review_status)
        a = appearances[appearances.report_key.eq(key)]
        for pr in partial_reviews.itertuples():
            if 'doi:'+pr.doi == key:
                partial.append(pr.review_status)
                evidence.append(pr.source_note)
                history.append(pr.decision)
        expansion_needs = a.screening_decision.isin(fulltext_statuses).any()
        for f in followup.itertuples():
            if 'doi:'+f.doi == key:
                read = read or f.primary_status == 'methods_results_read'
                evidence.append('database/alan_experiments/search/'+f.review_file)
                history.append(f.decision)
        pending = bool((pending_original or pending_supplement or expansion_needs) and not read)
        rows.append(dict(report_key=key,
            publication_family_key=report_lookup.loc[key,'publication_family_key'] if in_search else key,
            in_combined_search=in_search, in_original_search=key in original.index,
            in_original_outside_source_ledger=in_citation,
            primary_read_documented=bool(read), pending_documented_fulltext_followup=pending,
            partial_review_status=';'.join(dict.fromkeys(x for x in partial if x)),
            coverage_status='primary_read_not_final_eligibility' if read else
                'pending_documented_fulltext_followup' if pending else
                'no_fulltext_requirement_in_joined_records_not_final_exclusion',
            screening_history=';'.join(dict.fromkeys(list(a.screening_decision)+[x for x in history if x])),
            evidence_sources=';'.join(dict.fromkeys(x for x in evidence if x))))
    result = pd.DataFrame(rows)
    result.to_csv(ROOT/'unified_review_inventory.csv', index=False)
    pending = result[result.pending_documented_fulltext_followup]
    pending.to_csv(ROOT/'unified_pending_primary_review.csv', index=False)
    summary = dict(combined_search_reports=len(reports),
        supplementary_ledger_reports=len(supplementary),
        supplementary_reports_already_in_combined_search=len(set(supplementary.index)&set(report_lookup.index)),
        additional_reports_outside_both_searches=len(set(supplementary.index)-set(report_lookup.index)),
        union_reports=len(result), publication_families_after_known_links=result.publication_family_key.nunique(),
        primary_read_documented_reports=int(result.primary_read_documented.sum()),
        pending_documented_fulltext_followup_reports=len(pending),
        no_documented_current_fulltext_followup_reports=int((~result.primary_read_documented & ~result.pending_documented_fulltext_followup).sum()),
        final_eligibility=False, final_prisma=False,
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
        limits=['Original 160 read retained reports differs from all161 read original reports',
                'A primary read can support exclusion or contextual use, not necessarily experimental eligibility',
                'No fulltext requirement here is not a confirmed final exclusion',
                'Partial inspections do not become full reading; version links do not transfer reading',
                'Citation ledger is documented sources, not exhaustive citation-search hit count',
                'Programme-level independence and final adoption remain unadjudicated'])
    (ROOT/'unified_review_inventory_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
