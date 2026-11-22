#!/usr/bin/env python3
"""Bounded official Dryad audit of Burgazzi EERES-2021 DISCONNECTED 18-flume dataset.
No raw observations are re-published. No treatment effect is inferred at this stage.
"""
import csv, io, json, hashlib, os, re, zipfile, argparse
from collections import Counter
from pathlib import Path
from urllib.parse import urljoin, urlparse, quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

DOI='10.5061/dryad.rfj6q57k5'; B='https://datadryad.org/api/v2'
NAME='DISCONNECTED_alldata.xlsx'; MAX=2_000_000
HOSTS={'datadryad.org','www.datadryad.org'}

def grab(url):
    u=urlparse(url)
    if u.scheme!='https' or u.hostname not in HOSTS:raise ValueError('UNAPPROVED_HOST')
    with urlopen(Request(url,headers={'User-Agent':'ALAN-original-source-audit/2026-11-23','Accept':'application/json,*/*'}),timeout=28) as h:
        final=h.geturl();u=urlparse(final)
        if u.scheme!='https' or u.hostname not in HOSTS:raise ValueError('REDIRECTED_OUTSIDE_ALLOWLIST:'+str(u.hostname))
        if int(h.headers.get('Content-Length') or '0')>MAX:raise ValueError('PAYLOAD_TOO_LARGE')
        data=h.read(MAX+1)
        if len(data)>MAX:raise ValueError('PAYLOAD_TOO_LARGE')
        return data

def hrefs(obj):
    out=set()
    def w(o):
        if isinstance(o,dict):
            for k,v in o.items():
                if k in {'href','url','link'} and isinstance(v,str):
                    u=urljoin('https://datadryad.org',v);p=urlparse(u)
                    if p.scheme=='https' and p.hostname in HOSTS and p.path.startswith('/api/v2/'):out.add(u)
                elif isinstance(v,(dict,list)):w(v)
        elif isinstance(o,list):
            for v in o:w(v)
    w(obj);return out

def inspect_xlsx(content):
    """Source grain only: calculate sheet counts, headers and categorical keys."""
    try:
        import openpyxl
    except ImportError:return {'status':'SOURCE_STOP_MISSING_OPENPYXL'}
    book=openpyxl.load_workbook(io.BytesIO(content),read_only=True,data_only=True)
    report={}
    for name in book.sheetnames:
        ws=book[name]
        itr=ws.iter_rows(values_only=True)
        try:header=next(itr)
        except StopIteration: report[name]={'status':'EMPTY_SHEET'};continue
        cols=[str(v or '').strip() for v in header]
        field_n=len(cols);freq={k:Counter() for k in cols[:min(15,len(cols))] if re.search('ht|flume|treatment|side|section|time|date|flow|light|group|block',k,re.I)}
        n=0;sum_taxa=0;negative=0;na=0
        for row in itr:
            if not any(v is not None and str(v).strip() for v in row):continue
            n+=1
            for j,k in enumerate(cols[:min(15,len(cols))]):
                if k in freq and j<len(row):
                    v=row[j]
                    if v is not None and str(v).strip()!='':freq[k][str(v)[:65]]+=1
            if n>20000:raise ValueError('TOO_MANY_ROWS')
        report[name]={'n_nonempty_data_rows':n,'n_columns':field_n,'columns':cols[:80],
                      'physical_design_candidate_columns':{k:{'distinct':len(v),'top':v.most_common(20) if len(v)<=25 else v.most_common(6)} for k,v in freq.items()}}
    return {'status':'SOURCE_EXCEL_SHEETS_READ','sheets':report}

def find_archive(archive):
    if len(archive)>MAX:raise ValueError('ZIP_TOO_BIG')
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        entries=z.infolist()
        if len(entries)>30:raise ValueError('ZIP_ENTRIES_TOO_MANY')
        candidates=[x for x in entries if x.filename.split('/')[-1]==NAME]
        if len(candidates)!=1:raise ValueError('SOURCE_WORKBOOK_NOT_UNIQUE')
        info=candidates[0]
        if info.file_size>1_000_000:raise ValueError('XLSX_TOO_BIG')
        raw=z.read(info)
        if raw[:2]!=b'PK':raise ValueError('NOT_REAL_XLSX')
        return {'filename':info.filename,'sha256':hashlib.sha256(raw).hexdigest(),
                'size_bytes':len(raw),'metadata':inspect_xlsx(raw)}

def run():
    result={'review':'EERES2021_Disconnected','doi':DOI,'journal_publications':[
     '10.1016/j.scitotenv.2023.168836','10.1002/ecog.07758'],
     'expected':NAME,'raw_record_republication':False,'effect_adopted':False,
     'status':'SOURCE_STOP','metadata_retrieved':False,'probes':[]}
    for endpoint in [B+'/datasets/doi:'+DOI,B+'/datasets/doi:'+quote(DOI,safe='')]:
        try:
            raw=grab(endpoint);obj=json.loads(raw)
            result['dataset_metadata_sha256']=hashlib.sha256(raw).hexdigest()
            result['dataset_title']=obj.get('title')
            result['dataset_id']=obj.get('id')
            result['metadata_retrieved']=True
            result['probes'].append({'url':endpoint,'status':'JSON_RETRIEVED'})
            break
        except Exception as e:result['probes'].append({'url':endpoint,'status':'SOURCE_STOP','error':type(e).__name__+':'+str(e)[:100]})
    else:return result
    targets=sorted(hrefs(obj));result['metadata_links']=targets[:20]
    version_links=sorted(u for u in targets if re.search(r'/versions/\d+$',urlparse(u).path))
    # Verify version is advertised by official dataset metadata, not guessed.
    for version in version_links:
        try:
            content=grab(version+'/download')
            result['zip_bytes']=len(content);result['zip_sha256']=hashlib.sha256(content).hexdigest()
            result['original_workbook']=find_archive(content)
            result['status']='SOURCE_WORKBOOK_VERIFIED'
            return result
        except Exception as e:result['probes'].append({'url':version+'/download','status':'SOURCE_STOP','error':type(e).__name__+':'+str(e)[:130]})
    result['status']='SOURCE_STOP_ARCHIVE';return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='results/eeres2021_disconnected_source_audit.json');args=ap.parse_args()
    d=run();p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('EERES_STATUS',d['status'])
    if 'original_workbook' in d:
        f=d['original_workbook'];print('SOURCE_WORKBOOK',f['sha256'],f['size_bytes'])
        for nm,s in f['metadata'].get('sheets',{}).items():
            print('SHEET',nm,s.get('n_nonempty_data_rows'),s.get('n_columns'))
    print('AUDIT',p)
