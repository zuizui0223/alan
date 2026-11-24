#!/usr/bin/env python3
"""Bounded public original-data access audit for two 2026 ALAN insect studies.

Only downloads publicly released small data files from DataverseNL and Figshare.
Never treats an accessed file, a light colour, a survey night, or a species as
an independent illuminated plot. No effect-size calculations are attempted.
No raw source observations are committed; only manifest, hashes and schemas.
"""
import argparse
import csv
import hashlib
import io
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

LIMIT_BYTES = 5_000_000
MAX_FILES_PER_SOURCE = 6
USER_AGENT = 'ALAN-literature-source-provenance/2026-11-22 (+reproducible-research; no credential)'
ALLOWED_HOSTS = {'dataverse.nl', 'api.figshare.com', 'ndownloader.figshare.com', 'figshare.com'}
SOURCE = {
    'charvalakis2026': {
        'article_doi': '10.1016/j.biocon.2026.111711',
        'dataset_doi': '10.34894/SNKXWK',
        'endpoint': 'active camera light-trap catches and arrival timing',
        'design': '2023 Veluwe: four forest sites, twelve nights, eighteen traps/night, five spectra plus dark control',
        'host': 'dataversenl',
        'meta_urls': [
            'https://dataverse.nl/api/datasets/:persistentId/?persistentId=doi%3A10.34894%2FSNKXWK',
            'https://dataverse.nl/api/datasets/:persistentId/versions/:latest?persistentId=doi%3A10.34894%2FSNKXWK',
        ],
    },
    'reinhard2026': {
        'article_doi': '10.1038/s42003-026-10979-2',
        'dataset_doi': '10.6084/m9.figshare.29436527',
        'endpoint': 'lakeside fixed-roadlight HPS vs shielded-LED interception (no unlit control)',
        'design': '3 HPS + 3 shielded LED fixed roadlight fixtures, 13 nights of sampling (2023-2024)',
        'host': 'figshare',
        'meta_urls': ['https://api.figshare.com/v2/articles/29436527'],
    },
}


def check_url(url):
    u = urlparse(url)
    if u.scheme != 'https' or u.hostname not in ALLOWED_HOSTS:
        raise ValueError('Untrusted data source hostname/scheme: %s' % (u.hostname,))


def fetch_bytes(url, limit=LIMIT_BYTES):
    check_url(url)
    req = Request(url, headers={'User-Agent': USER_AGENT, 'Accept': 'application/json,text/csv,text/plain,*/*'})
    try:
        with urlopen(req, timeout=25) as resp:
            if int(resp.headers.get('Content-Length') or 0) > limit:
                raise ValueError('REMOTE_FILE_EXCEEDS_LIMIT')
            raw = resp.read(limit + 1)
            if len(raw) > limit:
                raise ValueError('REMOTE_FILE_EXCEEDS_LIMIT')
            return raw, resp.headers.get('Content-Type', '')
    except HTTPError as e:
        raise RuntimeError('HTTP_%s' % e.code) from e
    except URLError as e:
        raise RuntimeError('NETWORK_ERROR_%s' % type(e.reason).__name__) from e


def get_metadata(s):
    attempt=[]
    for url in s['meta_urls']:
        try:
            raw,typ=fetch_bytes(url, limit=2_000_000)
            obj=json.loads(raw)
            if not isinstance(obj,dict) or obj.get('status')=='ERROR':
                raise ValueError('API_BAD_STATUS_OR_SCHEMA')
            return obj, {'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'content_type':typ,'attempts':attempt}
        except (RuntimeError, ValueError, json.JSONDecodeError) as err:
            attempt.append({'url':url,'stop':str(err)[:240]})
    return None, {'attempts':attempt}


def dataverse_files(meta):
    payload=meta.get('data') or {}
    version=payload.get('latestVersion') or payload.get('datasetVersion') or payload
    files=version.get('files') if isinstance(version,dict) else None
    if files is None:files=payload.get('files',[])
    found=[]
    for v in files or []:
        df=v.get('dataFile') or {}
        fileid=df.get('id') or v.get('dataFileId')
        if not fileid:continue
        found.append({
            'label':v.get('label') or df.get('filename') or 'unknown',
            'id':fileid,
            'bytes_declared':df.get('filesize'),
            'restricted':df.get('restricted',v.get('restricted',False)),
            'mime_type':df.get('contentType'),
            'md5_from_repository':df.get('md5') or (df.get('checksum') or {}).get('value'),
            'download_url':f'https://dataverse.nl/api/access/datafile/{fileid}?format=original',
        })
    meta_license=version.get('license') or payload.get('license') or version.get('termsOfUse') or payload.get('termsOfUse')
    return found, meta_license


def figshare_files(meta):
    arr=meta.get('files',[])
    files=[]
    for f in arr:
        files.append({
            'label':f.get('name') or 'unknown',
            'id':f.get('id'),
            'bytes_declared':f.get('size'),
            'restricted':False,
            'mime_type':None,
            'md5_from_repository':f.get('supplied_md5') or f.get('computed_md5'),
            'download_url':f.get('download_url'),
        })
    return files,meta.get('license')


def tabular_summary(name, raw):
    """Inspect original tabular schema; never count records as independent lamps.

    The official Reinhard et al. source CSVs are semicolon-delimited R exports,
    which a default comma parser would silently misparse as one giant column.
    """
    low = name.lower()
    if low.endswith(('.csv', '.tsv')):
        try:
            st = raw.decode('utf-8-sig')
            first_line = st.splitlines()[0] if st else ''
            sep = '\t' if low.endswith('.tsv') else (';' if first_line.count(';') > first_line.count(',') else ',')
            reader = csv.DictReader(io.StringIO(st), delimiter=sep)
            head = reader.fieldnames
            if not head or len(head) > 500:
                return {'schema_status': 'NO_VALID_HEADER'}
            source_rows = list(reader)
            out = {
                'schema_status': 'CSV_PARSED', 'delimiter': 'tab' if sep == '\t' else sep,
                'row_count': len(source_rows), 'column_count': len(head), 'column_names': head[:70],
            }
            # Aggregation is descriptive only. The original file may contain
            # taxon × interval × repeated-trap rows including zero-filled records.
            # Do NOT infer n independent light installations from rows or dates.
            if all(c in head for c in ('count', 'treatment', 'trap', 'date')):
                groups = {}
                for rec in source_rows:
                    treatment = (rec.get('treatment') or 'MISSING').strip()
                    method = (rec.get('method') or 'MISSING').strip()
                    key = (treatment, method)
                    g = groups.setdefault(key, {'treatment': treatment, 'method': method,
                                                'source_rows': 0, 'count_numeric_rows': 0,
                                                'count_missing_or_non_numeric_rows': 0,
                                                'raw_count_column_sum': 0.0,
                                                '_traps': set(), '_dates': set(), '_intervals': set()})
                    g['source_rows'] += 1
                    if rec.get('trap'):g['_traps'].add(rec['trap'])
                    if rec.get('date'):g['_dates'].add(rec['date'])
                    if rec.get('interval'):g['_intervals'].add(rec['interval'])
                    try:
                        val = float(rec.get('count') or 'nan')
                        if not (-1e200 < val < 1e200) or val < 0:
                            raise ValueError('invalid count')
                        g['count_numeric_rows'] += 1
                        g['raw_count_column_sum'] += val
                    except ValueError:
                        g['count_missing_or_non_numeric_rows'] += 1
                summary = []
                for key, g in sorted(groups.items()):
                    g['n_trap_labels'] = len(g.pop('_traps'))
                    g['n_date_labels'] = len(g.pop('_dates'))
                    g['n_interval_labels'] = len(g.pop('_intervals'))
                    g['raw_count_column_sum'] = round(g['raw_count_column_sum'], 6)
                    summary.append(g)
                out['raw_count_by_treatment_method'] = summary[:30]
                out['group_summary_truncated'] = len(summary) > 30
                out['independence_warning'] = ('Raw column sums and date/trap labels are observational file summaries; '
                    'neither counts nor rows identify randomized lamp installations, '
                    'time-varying assignment, zero-control counterfactual or ALAN causal response.')
            return out
        except (UnicodeDecodeError, csv.Error, StopIteration):
            return {'schema_status': 'CSV_DECODE_STOP'}
    if low.endswith('.xlsx'):
        return {'schema_status': 'XLSX_NOT_PARSED'}
    return {'schema_status': 'OTHER_FORMAT_NOT_PARSED'}


def audit_one(code, spec, download=True):
    result={k:spec[k] for k in ['article_doi','dataset_doi','endpoint','design','host']}
    result['status']='SOURCE_STOP';result['meta_verified']=False
    payload,receipt=get_metadata(spec)
    result['metadata_receipt']=receipt
    if payload is None:
        result['stop_reason']='ALL_OFFICIAL_METADATA_ENDPOINTS_FAILED'
        return result
    result['meta_verified']=True
    if spec['host']=='dataversenl':
        files,license_info=dataverse_files(payload)
    else:
        files,license_info=figshare_files(payload)
    result['license_metadata']=license_info
    result['file_count']=len(files)
    result['file_registry']=[]
    confirmed=0
    for f in files:
        info={k:v for k,v in f.items() if k!='download_url'}
        size=f.get('bytes_declared')
        is_tabular=f['label'].lower().endswith(('.csv','.tsv','.xlsx'))
        info['download_status']='SKIP_NOT_TABULAR_OR_TOO_LARGE_OR_RESTRICTED'
        if download and is_tabular and not f.get('restricted') and (size is None or size<=LIMIT_BYTES) and confirmed<MAX_FILES_PER_SOURCE:
            try:
                if not f['download_url']:raise ValueError('NO_DOWNLOAD_URL')
                raw,ct=fetch_bytes(f['download_url'])
                info['download_status']='ORIGINAL_BYTES_VERIFIED'
                info['sha256_original']=hashlib.sha256(raw).hexdigest()
                info['downloaded_bytes']=len(raw)
                info['content_type']=ct
                info.update(tabular_summary(f['label'],raw))
                confirmed+=1
            except (RuntimeError,ValueError) as ex:
                info['download_status']='SOURCE_STOP_'+str(ex)[:130]
        result['file_registry'].append(info)
    result['downloaded_tabular_files']=confirmed
    result['status']='ORIGINAL_TABULAR_SOURCE_VERIFIED' if confirmed else 'METADATA_ONLY_SOURCE_GATE'
    if not files:result['stop_reason']='NO_PUBLIC_FILES_IN_RETRIEVED_MANIFEST'
    result['effect_sizes_estimated']=0
    result['independent_study_k_increase']=0
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',default='results/alan_night_insect_20261122_public_source_gate.json')
    ap.add_argument('--metadata-only',action='store_true')
    args=ap.parse_args()
    out={'date':'2026-11-22','scope':'Public original-data schema & checksum only; NOT causal meta-effect extraction', 'generated_at_utc_epoch':int(time.time()), 'sources':{}}
    for code, spec in SOURCE.items():
        out['sources'][code]=audit_one(code,spec,download=not args.metadata_only)
        print(code, out['sources'][code]['status'],'files',out['sources'][code].get('file_count',None),'verified_tabular',out['sources'][code].get('downloaded_tabular_files',0))
    p=Path(args.output);p.parent.mkdir(exist_ok=True,parents=True)
    p.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('MANIFEST',p)

if __name__=='__main__':main()
