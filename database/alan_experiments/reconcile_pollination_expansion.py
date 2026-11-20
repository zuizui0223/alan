"""Inventory a supplementary search without pretending it is fully screened."""
from pathlib import Path
import csv,json,hashlib

ROOT=Path(__file__).resolve().parent
CACHE=ROOT.parents[2]/'gbif_annual/alan_experiment_sources/pollination_search_expansion'
raw=CACHE/'response.json'
j=json.loads(raw.read_text(encoding='utf-8'))
def read(p):
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
old=read(ROOT/'search/record_identity_crosswalk.csv')
ids={r['record_id'] for r in old}
keys={r['report_key'] for r in old}
known={r['doi'].lower() for r in read(ROOT/'studies.csv')+read(ROOT/'search/fulltext_decisions_incremental.csv') if r['doi']}
rows=[]
for x in j['resultList']['result']:
    rid=x['source']+':'+x['id'];doi=x.get('doi','').lower()
    match='original_record_id' if rid in ids else 'original_report_doi' if doi and 'doi:'+doi in keys else 'supplementary_registry_doi' if doi and doi in known else 'identity_review_needed'
    rows.append(dict(record_id=rid,doi=doi,pmcid=x.get('pmcid',''),title=x['title'],identity_status=match,screening_status='not_newly_screened',primary_url='https://europepmc.org/article/'+x['source']+'/'+x['id']))
assert len(rows)==j['hitCount']==51 and len({r['record_id'] for r in rows})==51
out=ROOT/'search/pollination_expansion_records.csv'
with out.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
counts={s:sum(r['identity_status']==s for r in rows) for s in sorted({r['identity_status'] for r in rows})}
audit=dict(query=(CACHE/'query.txt').read_text(encoding='utf-8'),search_date='2026-11-20',response_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),hit_count=51,returned_records=51,identity_counts=counts,final_screening=False,original_521_records_unchanged=True,note='Exact IDs/DOIs only; unmatched records are not necessarily new papers or eligible experiments. Supplemental search removes experiment/manipulation term requirement.')
(ROOT/'search/pollination_expansion_manifest.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(counts)
