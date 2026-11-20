"""Join individually reviewed records without treating records as studies."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import re
import unicodedata

ROOT = Path(__file__).resolve().parent
FOLDER = ROOT/'search'


def read(name):
    with (FOLDER/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def write(name, fields, rows):
    with (FOLDER/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    records=read('europepmc_records.csv')
    screens=[]
    ranges=[('screen_000_173.csv',0,174),('screen_174_347.csv',174,348),('screen_348_520.csv',348,521)]
    for filename,start,end in ranges:
        subset=read(filename)
        assert [x['record_id'] for x in subset]==[x['record_id'] for x in records[start:end]], filename
        screens.extend(subset)
    allowed={'fulltext_needed','separate_observational','reference_mining','separate_indoor','exclude'}
    assert len(screens)==len(records)==521
    assert len(set(r['record_id'] for r in screens))==521
    assert all(r['decision'] in allowed and r['reason'].strip() and r['evidence_basis'] in {'title_abstract','title_only'} for r in screens)
    # Preserve initial judgments; explicit scope adjudications can retain records
    # that were too narrowly excluded by taxon or should be in spectrum-only.
    overrides=read('screening_adjudications.csv') if (FOLDER/'screening_adjudications.csv').exists() else []
    by_id={r['record_id']:r for r in screens}
    for change in overrides:
        assert change['record_id'] in by_id
        assert change['decision'] in allowed | {'separate_spectrum'}
        by_id[change['record_id']].update(decision=change['decision'],reason=change['reason'],reviewer=change['reviewer'])
    joined=[{**r,**s,'screening_status':s['decision'],'exclusion_reason':s['reason'] if s['decision']=='exclude' else ''}
            for r,s in zip(records,screens)]
    write('screened_records.csv',list(joined[0]),joined)
    # Metadata-only flags for manual reconciliation: do not silently merge titles.
    titles=defaultdict(list)
    for r in records:
        title=re.sub(r'[^\w]+','',unicodedata.normalize('NFKC',r['title']).casefold())
        titles[title].append(r)
    possible=[]
    for title,group in titles.items():
        if title and len(group)>1:
            for r in group:
                possible.append(dict(title_key=title,record_id=r['record_id'],doi=r['doi'],
                    title=r['title'],status='possible_duplicate_manual_resolution_pending'))
    write('possible_title_duplicates.csv',['title_key','record_id','doi','title','status'],possible)
    result=dict(raw_index_records=521,doi_deduplicated_records=521,
        title_abstract_screened_records=len(screens),decision_counts=dict(Counter(r['decision'] for r in screens)),
        evidence_counts=dict(Counter(r['evidence_basis'] for r in screens)),
        possible_title_duplicate_records=len(possible),unique_studies='not_yet_established',
        explicit_scope_adjudications=len(overrides),
        review_method='AI-assisted individual title/abstract screening; no human dual review claimed',
        fulltext_eligibility='pending',WoS='not_executed; browser endpoint returned 403',
        Scopus='not_executed; browser endpoint returned 403',
        no_missing_abstract_auto_exclusion=True,
        raw_response_sha256=[r['sha256'] for r in read('europepmc_pages.csv')])
    (FOLDER/'screening_summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
