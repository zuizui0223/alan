#!/usr/bin/env python3
"""Primary published effect extraction; not a refit of individual observations.

Giavi 2020 original JATS coefficient tables retain uncertainty and contrast
orientation. Guenat 2023 supplemental study inventory audits light/ALAN scope.
No pooled mean is calculated across different endpoints or exposures.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

BASE='https://www.ebi.ac.uk/europepmc/webservices/rest/'

def text(node):
    return ' '.join(''.join(node.itertext()).split()) if node is not None else ''

def effect(beta, se):
    if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
        raise ValueError('Finite estimate and positive measured standard error required')
    return {'log_estimate':beta,'standard_error':se,'log_variance':se*se,
            'ratio':math.exp(beta),
            'ratio_ci95':[math.exp(beta-1.959963984540054*se),math.exp(beta+1.959963984540054*se)],
            'interval_method':'Wald normal interval reconstructed from published estimate and SE; not new model fitting'}

def common_endpoint(rows):
    names={r['endpoint'] for r in rows}
    if len(names)!=1:raise ValueError('Unlike biological endpoints must not be pooled')
    return names.pop()

def table_rows(table):
    return [[{'text':text(c),'rowspan':int(c.get('rowspan','1')),
              'colspan':int(c.get('colspan','1'))} for c in r if c.tag in ('td','th')]
            for r in table.findall('.//tr')]

def get(url, maxbytes=20000000):
    req=urllib.request.Request(url,headers={'User-Agent':'ALAN-primary-effect-source/1.0'})
    with urllib.request.urlopen(req,timeout=75) as r:raw=r.read(maxbytes+1)
    if len(raw)>maxbytes:raise ValueError('Source exceeds bounded extraction size')
    return raw

def jats(pmc):
    raw=get(BASE+pmc+'/fullTextXML')
    root=ET.fromstring(raw)
    title=text(root.find('.//article-title'))
    doi=[text(e) for e in root.findall('.//article-id') if e.get('pub-id-type')=='doi']
    tables=[{'label':text(t.find('label')),'caption':text(t.find('caption')),
             'rows':table_rows(t.find('table'))} for t in root.findall('.//table-wrap') if t.find('table') is not None]
    return root, {'pmc':pmc,'title':title,'doi':doi,'jats_sha256':hashlib.sha256(raw).hexdigest(),'tables':tables}

def docx_data(raw):
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:root=ET.fromstring(z.read('word/document.xml'))
    paragraphs=[' '.join(''.join(p.itertext()).split()) for p in root.findall('.//w:body/w:p',ns)]
    tables=[]
    for table in root.findall('.//w:tbl',ns):
        rows=[]
        for row in table.findall('w:tr',ns):
            rows.append([' '.join(''.join(c.itertext()).split()) for c in row.findall('w:tc',ns)])
        tables.append(rows)
    return {'docx_sha256':hashlib.sha256(raw).hexdigest(),'paragraphs':paragraphs,'tables':tables}

def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    _,primary=jats('PMC7368033')
    if '10.1038/s41598-020-68667-y' not in primary['doi']:
        raise ValueError('Wrong primary DOI')
    (out/'giavi2020_original_coefficient_tables.json').write_text(json.dumps(primary,indent=2))
    print('PRIMARY_TABLES',json.dumps(primary),flush=True)
    root,review=jats('PMC10111172')
    if '10.1002/ece3.9990' not in review['doi']:raise ValueError('Wrong review DOI')
    # Original text relevant to study scope, not a new outcome source.
    review['scope_paragraphs']=[text(p) for p in root.findall('.//body//p')
        if any(k in text(p).lower() for k in ('uv-b','light pollutants','light pollution','22 articles','22 studies'))]
    try:
        package=get(BASE+'PMC10111172/supplementaryFiles')
        with zipfile.ZipFile(io.BytesIO(package)) as z:
            names=z.namelist(); docs=[n for n in names if n.lower().endswith('.docx')]
            review['supplement_files']=names
            review['docx']=[{'member':n,**docx_data(z.read(n))} for n in docs]
        review['supplement_status']='ACTUAL_DOCX_READ'
    except Exception as e:
        review['supplement_status']='NOT_ACQUIRED';review['supplement_error']=str(e)
    (out/'guenat2023_scope_source_inventory.json').write_text(json.dumps(review,indent=2))
    print('REVIEW_SCOPE',json.dumps({k:v for k,v in review.items() if k!='docx'}),flush=True)
    for doc in review.get('docx',[]):
        print('SUPPLEMENT_PARAGRAPHS',json.dumps(doc['paragraphs']),flush=True)
        for i,table in enumerate(doc['tables']):print('SUPPLEMENT_TABLE',i,json.dumps(table),flush=True)
    print('PRIMARY_SOURCE_EXTRACTION_COMPLETE_NO_RAW_REFIT_OR_GLOBAL_META_MEAN',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    main(p.parse_args().out)
