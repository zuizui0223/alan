#!/usr/bin/env python3
"""Audit official Dryad publication files for Sullivan et al. 2019.

Do not confuse the observational 2010/2011 gradients with the 2013 lighting
field experiment; don't treat individual arthropods or visits as light-site n.
No raw author records are published by this audit. Fail closed on blocked bytes.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, os, re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote, urlparse, unquote
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

DOI = '10.5061/dryad.gb9534p'
PAPER = '10.1002/eap.1821'
API = 'https://datadryad.org/api/v2'
DATASET_URLS = (
    'https://datadryad.org/api/v2/datasets/doi:' + DOI,
    'https://datadryad.org/api/v2/datasets/doi:' + quote(DOI, safe=''),
)
EXPECTED = ('Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt',
            'README_for_Sullivan_EAP18-0555_WetlandFieldExperiment_2013.txt')
MAX_FETCH = 450_000
UA = 'ALAN-scientific-source-audit/2026-11-22 (public source; no login)'


def fetch(url):
    u = urlparse(url)
    if u.scheme != 'https' or u.hostname not in {'datadryad.org', 'www.datadryad.org'}:
        raise ValueError('Non-Dryad API URL forbidden')
    req = Request(url, headers={'User-Agent': UA, 'Accept': 'application/json, text/plain, */*'})
    with urlopen(req, timeout=20) as resp:
        final = resp.geturl()
        if (int(resp.headers.get('Content-Length') or 0) > MAX_FETCH):
            raise ValueError('source too large')
        raw = resp.read(MAX_FETCH + 1)
        if len(raw) > MAX_FETCH:
            raise ValueError('response beyond cap')
        return {'bytes': raw, 'final_url': final, 'status': resp.status,
                'mime': resp.headers.get('Content-Type','')}


def small_json(url, errors):
    try:
        response = fetch(url)
        raw = response['bytes']
        d = json.loads(raw)
        return d, {'request_url': url, 'status': 'JSON_RECEIVED', 'byte_count': len(raw),
                   'sha256': hashlib.sha256(raw).hexdigest(), 'final_url': response['final_url']}
    except (HTTPError, URLError, OSError, ValueError, json.JSONDecodeError) as e:
        issue = {'request_url': url, 'status': 'SOURCE_STOP',
                 'reason': f'{type(e).__name__}:{getattr(e,"code",str(e))}'[:180]}
        errors.append(issue)
        return None, issue


def links(d):
    """Return Dryad API urls only (bounded, no blind external navigation)."""
    uris = set()
    if isinstance(d,dict):
        for k,v in d.items():
            if k in ('href', 'url', 'link') and isinstance(v,str):
                u = urlparse(v)
                if u.scheme == 'https' and u.hostname in {'datadryad.org', 'www.datadryad.org'} and ('/api/v2/' in u.path):
                    uris.add(v)
            elif isinstance(v,(list,dict)):
                uris |= links(v)
    elif isinstance(d,list):
        for v in d: uris |= links(v)
    return uris


def file_summaries(node):
    found = []
    def visit(x):
        if isinstance(x,dict):
            name = x.get('path') or x.get('filename') or x.get('name')
            if isinstance(name,str) and ('WetlandFieldExperiment' in name or 'UrbanStream' in name):
                found.append({'name':name, 'id':x.get('id'), 'size':x.get('size'),
                              'content_type':x.get('mimeType') or x.get('mime_type'),
                              'has_download_link': bool([z for z in links(x) if '/files/' in z])})
            for v in x.values():
                if isinstance(v,(dict,list)):visit(v)
        elif isinstance(x,list):
            for v in x:visit(v)
    visit(node)
    return found


def summarize_text(content, filename):
    s=content.decode('utf-8-sig','replace')
    if filename.startswith('README_'):
        lines=[v.strip() for v in s.splitlines() if re.search(r'plot|replic|light|lamp|treat|site|habitat|control|wetland|date',v,re.I)]
        return {'format':'README_text','total_lines':len(s.splitlines()),'relevant_lines_first_15':[v[:180] for v in lines[:15]]}
    head=s[:10000]
    try:
        delim=csv.Sniffer().sniff(head,delimiters='\t,;').delimiter
    except csv.Error: delim='\t'
    dat=list(csv.reader(io.StringIO(s), delimiter=delim))
    if not dat:return {'format':'empty_table'}
    header=dat[0]
    record_rows=[r for r in dat[1:] if any(x.strip() for x in r)]
    if not 1<len(header)<200 or not header:return {'format':'unsupported_table','line_count':len(dat)}
    plausible={}
    for index,label in enumerate(header):
        if re.search(r'light|lux|treat|site|habitat|plot|block|date|time|month|type|trap',label,re.I):
            vals=[r[index].strip() for r in record_rows if len(r)>index and r[index].strip()]
            n=len(set(vals))
            plausible[label]={'observed':len(vals),'distinct':n, 'top_categories':Counter(vals).most_common(10) if n<=30 else []}
    return {'format':'delimited_table','delimiter':repr(delim), 'columns':header,
            'nonempty_rows':len(record_rows),'column_count':len(header),
            'physical_unit_candidate_fields':plausible}


def run():
    errors=[]; log=[]; root=None
    for url in DATASET_URLS:
        d,record=small_json(url,errors);log.append(record)
        if isinstance(d,dict) and ('_links' in d or 'title' in d or 'identifier' in d):
            root=d;break
    output={'doi':DOI,'publication':PAPER,'source':'original Dryad API', 'probes':log,
            'declared_expected_files':list(EXPECTED),'source_files':[],
            'new_meta_eligible_independent_programme_k':0,
            'within_site_moisture_meta_regression':False,
            'source_reanalysis':False,
            'note':'2013 wetland experiment distinct from 2010-2011 observational stream gradients; only physical independent light sites count'}
    if root is None:
        output['status']='SOURCE_STOP_DATASET_METADATA_UNAVAILABLE';return output
    output['title']=root.get('title');output['dataset_id']=root.get('id')
    queue=list(sorted(links(root)))
    # Old v2 REST often has /versions or /files links; keep only limited metadata hops.
    candidates=[];seen=set()
    while queue and len(seen)<12:
        url=queue.pop(0)
        if url in seen:continue
        seen.add(url)
        if not ('/versions/' in url or '/datasets/' in url or '/files' in url):continue
        d,rec=small_json(url,errors);log.append(rec)
        if d is None:continue
        candidates.extend(file_summaries(d))
        queue.extend(sorted(links(d)-seen-set(queue)))
    output['found_file_metadata']=candidates[:50]
    # Only read small, original expected files via official file IDs, after ID verified in manifest.
    filenames={}
    for f in candidates:
        name=f.get('name','')
        if any(name.endswith(e) for e in EXPECTED) and isinstance(f.get('id'),int):
            filenames[name.split('/')[-1]]=f['id']
    for filename,id_ in filenames.items():
        u=f'{API}/files/{id_}/download'
        try:
            r=fetch(u);raw=r['bytes']
            if raw[:1]==b'{' or b'<html' in raw[:300].lower():
                output['source_files'].append({'name':filename,'status':'SOURCE_STOP_UNEXPECTED_JSON_OR_HTML'});continue
            output['source_files'].append({'name':filename,'status':'SOURCE_BYTES_VERIFIED','byte_count':len(raw),
                      'sha256':hashlib.sha256(raw).hexdigest(), **summarize_text(raw,filename)})
        except (HTTPError, URLError, ValueError, OSError) as e:
            output['source_files'].append({'name':filename,'status':'SOURCE_STOP', 'reason':f'{type(e).__name__}:{getattr(e,"code",str(e))}'[:140]})
    output['metadata_probes_count']=len(log)
    output['status']='SOURCE_FILE_METADATA_FOUND' if candidates else 'SOURCE_STOP_NO_FILE_MANIFEST'
    if output['source_files'] and all(x['status']=='SOURCE_BYTES_VERIFIED' for x in output['source_files']):
        output['status']='SOURCE_BYTES_READ_FOR_SCHEMA_AUDIT'
    return output


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default='results/alan_sullivan2019_dryad_source_gate_20261122.json');a=p.parse_args()
    d=run();path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('SOURCE_GATE_STATUS',d['status'],'manifest_files',len(d.get('found_file_metadata',[])),'files_acquired',len(d['source_files']))
    for f in d['source_files']:print('SOURCE_FILE',f['name'],f['status'],f.get('nonempty_rows'),f.get('columns',[])[:15])
    print('RESULT',path)
if __name__=='__main__':main()
