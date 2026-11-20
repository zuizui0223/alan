"""Reproducible supplementary index search, not a WoS/Scopus replacement.

Raw responses, including abstracts, stay in the external cache. Published
tables contain identifiers, citation metadata and pending screening status.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys
from datetime import datetime, timezone
import requests

OUT = Path(__file__).resolve().parent / 'search'
ENDPOINT = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
QUERY = ('(TITLE_ABS:"artificial light at night" OR TITLE_ABS:"light pollution" '
         'OR TITLE_ABS:"street lighting" OR TITLE_ABS:"nocturnal illumination") '
         'AND (experiment* OR manipulat* OR mesocosm*) '
         'AND (plant* OR insect* OR arthropod* OR pollinat* OR spider* '
         'OR herbivor* OR predator* OR emergence) '
         'AND FIRST_PDATE:[1900-01-01 TO 2026-11-20]')


def write(name, rows):
    with (OUT/name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    cache = Path(sys.argv[1])
    cache.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(exist_ok=True)
    cursor = '*'
    seen_cursors = set()
    records, pages = [], []
    started = datetime.now(timezone.utc).isoformat()
    total = None
    for index in range(100):
        if cursor in seen_cursors:
            raise ValueError('repeated cursor before complete retrieval')
        seen_cursors.add(cursor)
        params = dict(query=QUERY, format='json', pageSize=1000,
                      resultType='core', cursorMark=cursor)
        response = requests.get(ENDPOINT, params=params, timeout=90)
        response.raise_for_status()
        data = response.json()
        total = data['hitCount'] if total is None else total
        if data['hitCount'] != total:
            raise ValueError('index changed during pagination; rerun snapshot')
        raw = response.content
        path = cache/f'europepmc_{index:03}.json'
        path.write_bytes(raw)
        batch = data['resultList']['result']
        pages.append(dict(page=index, request_url=response.url, filename=path.name,
                          sha256=hashlib.sha256(raw).hexdigest(), records=len(batch)))
        records.extend(batch)
        if len(records) >= total:
            break
        if not batch:
            raise ValueError('premature empty page')
        cursor = data['nextCursorMark']
    if len(records) != total:
        raise ValueError('retrieved count does not match index count')
    source_ids = [f"{r['source']}:{r['id']}" for r in records]
    if len(source_ids) != len(set(source_ids)):
        raise ValueError('duplicate source IDs in cursor retrieval')
    groups = {}
    mapping, unique = [], []
    for r in records:
        doi = r.get('doi','').strip().lower().removeprefix('https://doi.org/')
        key = 'doi:'+doi if doi else f"{r['source']}:{r['id']}"
        ident = f"{r['source']}:{r['id']}"
        representative = groups.setdefault(key, ident)
        mapping.append(dict(record_id=ident, canonical_key=key, representative_id=representative,
                            duplicate_by_doi=int(ident!=representative)))
        if representative != ident:
            continue
        unique.append(dict(record_id=ident,doi=doi,title=r.get('title',''),
            year=r.get('pubYear',''),first_publication=r.get('firstPublicationDate',''),
            journal=r.get('journalInfo',{}).get('journal',{}).get('title',''),
            pmcid=r.get('pmcid',''),is_open_access=r.get('isOpenAccess',''),
            source_url='https://europepmc.org/article/'+r['source']+'/'+r['id'],
            screening_status='pending',exclusion_reason=''))
    write('europepmc_records.csv',unique)
    write('europepmc_dedup.csv',mapping)
    write('europepmc_pages.csv',pages)
    manifest = dict(database='Europe PMC',role='supplementary_index_search',
        query=QUERY,endpoint=ENDPOINT,started_utc=started,completed_utc=datetime.now(timezone.utc).isoformat(),
        raw_hits=total,doi_deduplicated_records=len(unique),screened=0,
        scope_limit='Europe PMC is biomedical/life-science biased; not WoS or Scopus; title/abstract light terms; 1900 onward',
        cutoff='2026-11-20',abstracts_in_external_cache_only=True)
    (OUT/'europepmc_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(manifest,indent=2))


if __name__ == '__main__':
    main()
