#!/usr/bin/env python3
"""Guenat original effect data, with ALAN and daytime UV kept separate.
Metadata and bytes were independently verified in source run 43594262699.
These frozen hashes permit public direct file retrieval without live metadata.
"""
import argparse,csv,hashlib,io,json,os,re,time,urllib.request,urllib.error
from urllib.parse import quote
from collections import Counter
from pathlib import Path
DATA=(7439705,'Meta-analysis data.csv','9af3fc55b053388b5992e0e80d6c1fbd8d')
CODE=(7521548,'Pollination-meta.R','18fb1c5c6027162802879a7152bc719ea5')
PINS={7439705:(26508,'0a28da7e3687277280e98c3f438c7caeb598f159944fd4dc06596fdc6ab48d706'),
      7521548:(22175,'feb59cc139a99c41f5b2f666addabae8db73b80df07e621a02aa781242fe6368b')}
def valid_md5(raw,expected):return hashlib.md5(raw).hexdigest()==expected

def get(url):
    for attempt in range(2):
        req=urllib.request.Request(url,headers={'User-Agent':'ALAN-primary-source-meta-audit/1.0'})
        try:
            with urllib.request.urlopen(req,timeout=40) as r:raw=r.read(2000001)
            if len(raw)>2000000:raise ValueError('Unexpected large source')
            return raw
        except urllib.error.HTTPError as e:
            if e.code not in (429,500,502,503,504) or attempt==1:raise
        except (urllib.error.URLError,TimeoutError):
            if attempt==1:raise
        time.sleep(2)
    raise RuntimeError('Source retrieval exhausted')

def acquire(spec):
    record,name,digest=spec
    size,sha=PINS[record]
    cache=Path(os.environ.get('RUNNER_TEMP','/tmp'))/('alan-source-'+digest)
    receipt_path=cache.with_suffix('.json')
    def verify(raw):
        if len(raw)!=size or not valid_md5(raw,digest) or hashlib.sha256(raw).hexdigest()!=sha:
            raise ValueError('Original source identity differs from verified archive')
    if cache.exists() and receipt_path.exists():
        raw=cache.read_bytes();verify(raw);receipt=json.loads(receipt_path.read_text())
        if receipt.get('record')==record and receipt.get('filename')==name:return raw,receipt
    url=f'https://zenodo.org/records/{record}/files/{quote(name,safe="")}?download=1'
    mode='official_record_file_with_frozen_verified_metadata'
    try:raw=get(url)
    except (urllib.error.URLError,TimeoutError):
        meta=json.loads(get(f'https://zenodo.org/api/records/{record}'))
        if int(meta.get('id',-1))!=record:raise ValueError('Wrong original record')
        files=[f for f in meta['files'] if f['key']==name]
        if len(files)!=1:raise ValueError('Original file missing or duplicated')
        f=files[0]
        if f['checksum']!='md5:'+digest or f['size']!=size:raise ValueError('Original metadata changed')
        url=f['links'].get('self') or f['links']['download']
        raw=get(url);mode='official_live_metadata_file_link'
    verify(raw)
    receipt={'record':record,'filename':name,'bytes':len(raw),'md5':digest,'sha256':sha,
             'retrieval_mode':mode,'file_url':url,'initial_source_verification_run':43594262699}
    cache.write_bytes(raw);receipt_path.write_text(json.dumps(receipt))
    return raw,receipt

def parse_csv(raw):
    try:txt=raw.decode('utf-8-sig')
    except UnicodeDecodeError:txt=raw.decode('cp1252')
    candidates=[]
    for delimiter in (',',';','\t'):
        try:rows=[r for r in csv.reader(io.StringIO(txt),delimiter=delimiter,strict=True) if r]
        except csv.Error:continue
        if rows and len(rows[0])>1 and all(len(r)==len(rows[0]) for r in rows):candidates.append(rows)
    if len(candidates)!=1:raise ValueError('CSV delimiter ambiguous or inconsistent row widths')
    rows=candidates[0];fields=rows[0]
    if len(set(fields))!=len(fields):raise ValueError('Duplicated original CSV columns')
    return fields,[dict(zip(fields,r)) for r in rows[1:]]

def load_source():
    raw,receipt=acquire(DATA)
    fields,rows=parse_csv(raw)
    return fields,rows,receipt

def inventory(fields,rows):
    out={}
    for col in fields:
        vals=[str(r.get(col,'')).strip() for r in rows];counter=Counter(vals);nonnum=[]
        for v in counter:
            try:float(v)
            except ValueError:nonnum.append(v)
        out[col]={'unique':len(counter),'missing':sum(v in ('','NA','NaN','nan') for v in vals),
                  'levels':dict(counter) if len(counter)<=25 else None,'text_examples':nonnum[:12]}
    return out

def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    fields,rows,receipt=load_source();code,code_receipt=acquire(CODE)
    try:code_txt=code.decode('utf-8-sig')
    except UnicodeDecodeError:code_txt=code.decode('cp1252')
    relevant=[{'line':i,'text':line} for i,line in enumerate(code_txt.splitlines(),1)
              if re.search(r'ALAN|UV|night|light|escalc|rma\.|random|vi=|sei|Effect|filter|read.csv',line,re.I)]
    result={'source':receipt,'code_source':code_receipt,'n_effect_rows':len(rows),'columns':fields,
            'inventory':inventory(fields,rows),'first_source_effect_rows':rows[:3],
            'original_code_relevant_lines':relevant,'not_a_new_raw_plant_experiment':True,
            'new_meta_effect_not_fitted_yet':True}
    (out/'guenat_original_meta_source_inventory.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    print('ORIGINAL_META_INVENTORY',json.dumps(result,ensure_ascii=False),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);main(p.parse_args().out)
