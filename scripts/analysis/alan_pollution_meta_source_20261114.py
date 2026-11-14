#!/usr/bin/env python3
"""Read genuine published Guenat meta dataset, without treating UV as ALAN.

Original data DOI 10.5281/zenodo.7439705; code DOI 10.5281/zenodo.7521548.
Source rows are published effect estimates, not experimental individuals.
"""
import argparse,csv,hashlib,io,json,re,urllib.request
from collections import Counter
from pathlib import Path
DATA=(7439705,'Meta-analysis data.csv','9af3fc55b053388b5992e0e80d6c1fbd8d')
CODE=(7521548,'Pollination-meta.R','18fb1c5c6027162802879a7152bc719ea5')
def valid_md5(raw,expected):return hashlib.md5(raw).hexdigest()==expected

def get(url):
    req=urllib.request.Request(url,headers={'User-Agent':'ALAN-primary-source-meta-audit/1.0'})
    with urllib.request.urlopen(req,timeout=60) as r:raw=r.read(2000001)
    if len(raw)>2000000:raise ValueError('Unexpected large source')
    return raw

def acquire(spec):
    record,name,digest=spec
    meta=json.loads(get(f'https://zenodo.org/api/records/{record}'))
    if int(meta.get('id',-1))!=record:raise ValueError('Wrong original record')
    files=[f for f in meta['files'] if f['key']==name]
    if len(files)!=1:raise ValueError('Original file missing or duplicated')
    f=files[0]
    if f['checksum']!='md5:'+digest:raise ValueError('Original metadata digest changed')
    raw=get(f['links'].get('self') or f['links']['download'])
    if not valid_md5(raw,digest):raise ValueError('Original source MD5 mismatch')
    if len(raw)!=f['size']:raise ValueError('Source file truncated')
    return raw,{'record':record,'filename':name,'bytes':len(raw),'md5':digest,'sha256':hashlib.sha256(raw).hexdigest()}

def parse_csv(raw):
    try:txt=raw.decode('utf-8-sig')
    except UnicodeDecodeError:txt=raw.decode('cp1252')
    dialect=csv.Sniffer().sniff(txt[:12000],delimiters=',;\t')
    reader=csv.DictReader(io.StringIO(txt),dialect=dialect)
    return reader.fieldnames,list(reader)

def load_source():
    raw,receipt=acquire(DATA)
    fields,rows=parse_csv(raw)
    return fields,rows,receipt

def inventory(fields,rows):
    out={}
    for col in fields:
        vals=[str(r.get(col,'')).strip() for r in rows]
        counter=Counter(vals)
        nonnum=[]
        for v in counter:
            try:float(v)
            except ValueError:nonnum.append(v)
        out[col]={'unique':len(counter),'missing':sum(v in ('','NA','NaN','nan') for v in vals),
                  'levels':dict(counter) if len(counter)<=25 else None,
                  'text_examples':nonnum[:12]}
    return out

def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    fields,rows,receipt=load_source()
    code,code_receipt=acquire(CODE)
    try:code_txt=code.decode('utf-8-sig')
    except UnicodeDecodeError:code_txt=code.decode('cp1252')
    relevant=[{'line':i,'text':line} for i,line in enumerate(code_txt.splitlines(),1)
              if re.search(r'ALAN|UV|night|light|escalc|rma\.|random|vi=|sei|Effect|filter|read.csv',line,re.I)]
    result={'source':receipt,'code_source':code_receipt,'n_effect_rows':len(rows),'columns':fields,
            'inventory':inventory(fields,rows),'first_source_effect_rows':rows[:3],
            'original_code_relevant_lines':relevant,
            'not_a_new_raw_plant_experiment':True,'new_meta_effect_not_fitted_yet':True}
    (out/'guenat_original_meta_source_inventory.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    print('ORIGINAL_META_INVENTORY',json.dumps(result,ensure_ascii=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    main(p.parse_args().out)
