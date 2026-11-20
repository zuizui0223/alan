"""Acquire public PMC primary texts for pending records; acquisition != eligibility.

Full XML remains outside git. Preserve HTTP failures and article identity checks.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import time
import xml.etree.ElementTree as ET
import requests

ROOT=Path(__file__).resolve().parent


def normalize_doi(value):
    return (value or '').strip().lower().removeprefix('https://doi.org/').removeprefix('http://dx.doi.org/')


def acquire(row,cache):
    result=dict(record_id=row['record_id'],doi=row['doi'],pmcid=row['pmcid'],
        source_url='',status='',sha256='',xml_doi='',body_characters='',cache_filename='',error='')
    if not row['pmcid']:
        result['status']='publisher_or_repository_retrieval_needed'
        return result
    pmcid=row['pmcid']
    url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'
    result['source_url']=url
    path=cache/(pmcid+'.xml')
    try:
        if path.exists():
            raw=path.read_bytes()
        else:
            for attempt in range(3):
                response=requests.get(url,timeout=45)
                if response.status_code in (429,500,502,503,504) and attempt<2:
                    time.sleep(1+attempt)
                    continue
                response.raise_for_status()
                break
            raw=response.content
        tree=ET.fromstring(raw)
        if tree.tag!='article' or tree.find('body') is None:
            raise ValueError('response has no article body')
        ids={e.get('pub-id-type'):(e.text or '').strip() for e in tree.findall('./front/article-meta/article-id')}
        actual=ids.get('pmc') or ids.get('pmcid')
        if actual and actual.removeprefix('PMC')!=pmcid.removeprefix('PMC'):
            raise ValueError('PMC identity mismatch')
        doi=normalize_doi(ids.get('doi'))
        if doi and row['doi'] and normalize_doi(row['doi'])!=doi:
            raise ValueError('DOI identity mismatch')
        body=' '.join(tree.find('body').itertext())
        if len(body)<200:
            raise ValueError('body too short; inspect manually')
        path.write_bytes(raw)
        result.update(status='retrieved_body_pending_fulltext_review',sha256=hashlib.sha256(raw).hexdigest(),
            xml_doi=doi,body_characters=len(body),cache_filename=path.name)
    except (requests.RequestException,ValueError,ET.ParseError) as exc:
        result.update(status='retrieval_failed_not_excluded',error=str(exc)[:350])
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('cache',type=Path)
    args=parser.parse_args();args.cache.mkdir(parents=True,exist_ok=True)
    with (ROOT/'search/screened_records.csv').open(encoding='utf-8') as f:
        rows=[r for r in csv.DictReader(f) if r['decision']=='fulltext_needed']
    results={}
    with ThreadPoolExecutor(max_workers=4) as pool:
        pending={pool.submit(acquire,r,args.cache):r['record_id'] for r in rows}
        for n,future in enumerate(as_completed(pending),1):
            result=future.result();results[result['record_id']]=result
            if n%20==0:print(f'acquisition {n}/{len(rows)}',flush=True)
    ordered=[results[r['record_id']] for r in rows]
    out=ROOT/'search/fulltext_acquisition.csv'
    with out.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(ordered[0]));w.writeheader();w.writerows(ordered)
    from collections import Counter
    summary=dict(records=len(rows),status_counts=dict(Counter(r['status'] for r in ordered)),
                 interpretation='retrieval only; not fulltext eligibility decisions')
    (ROOT/'search/fulltext_acquisition_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
