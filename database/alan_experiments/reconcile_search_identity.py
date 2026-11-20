"""Resolve record identity from original metadata and acquired article DOI.

Same DOI identifies a report, not necessarily independent experimental data.
Preprint/journal versions with distinct DOI remain unresolved.
"""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent/'search'
d=pd.read_csv(ROOT/'screened_records.csv').fillna('')
a=pd.read_csv(ROOT/'fulltext_acquisition.csv').fillna('')
assert d.record_id.is_unique and a.record_id.is_unique
verified=a[a.status=='retrieved_body_pending_fulltext_review']
pmc_to_doi={r.pmcid:r.xml_doi for r in verified.itertuples() if r.xml_doi}
def norm(s):
    return s.strip().lower().removeprefix('https://doi.org/').removeprefix('http://dx.doi.org/')
rows=[]
for r in d.itertuples():
    original=norm(r.doi)
    fulltext=norm(pmc_to_doi.get(r.pmcid,''))
    assert not (original and fulltext and original!=fulltext),(r.record_id,original,fulltext)
    resolved=original or fulltext
    rows.append(dict(record_id=r.record_id,original_doi=original,resolved_doi=resolved,
        doi_evidence='index_metadata' if original else ('acquired_article_metadata' if fulltext else 'unresolved'),
        pmcid=r.pmcid,report_key='doi:'+resolved if resolved else ('pmc:'+r.pmcid if r.pmcid else r.record_id),
        title=r.title,decision=r.decision))
out=pd.DataFrame(rows)
out['records_sharing_report_key']=out.groupby('report_key').record_id.transform('size')
out.to_csv(ROOT/'record_identity_crosswalk.csv',index=False)
out[out.records_sharing_report_key>1].to_csv(ROOT/'confirmed_same_report_records.csv',index=False)
summary=dict(index_records=len(out),resolved_doi_records=int(out.resolved_doi.ne('').sum()),
 doi_filled_from_acquired_article=int((out.doi_evidence=='acquired_article_metadata').sum()),
 distinct_report_keys=out.report_key.nunique(),
 same_report_duplicate_records=len(out)-out.report_key.nunique(),
 independent_experiments='not established',
 warning='Exact DOI/PMCID report identity only; unmatched title duplicates, preprint-journal versions and shared experiments still require review. No screening decisions overwritten.')
(ROOT/'record_identity_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
