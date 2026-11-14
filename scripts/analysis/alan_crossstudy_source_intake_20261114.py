#!/usr/bin/env python3
"""Bounded acquisition of original open-access review supplements.

This is source intake, not a claim of new meta-analysis. No inferred biological
observations or digitized graph estimates. Original downloaded bytes stay
transient; checksums, named tables and original publication-level summaries
are retained for independent assessment.
"""
from __future__ import annotations
import argparse, hashlib, io, json, urllib.request, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

BASE = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10111172/'
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}

def get(url: str, limit: int = 25_000_000) -> bytes:
    req=urllib.request.Request(url,headers={'User-Agent':'ALAN-open-source-scientific-intake/1.0'})
    with urllib.request.urlopen(req,timeout=60) as resp:
        raw=resp.read(limit+1)
    if len(raw)>limit: raise ValueError('Source too large: '+url)
    return raw

def text(el):
    return ' '.join(' '.join(el.itertext()).split())

def word_content(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        doc=ET.fromstring(z.read('word/document.xml'))
    body=doc.find('w:body',NS)
    paragraphs=[]; tables=[]; last=[]
    for el in body:
        tag=el.tag.rsplit('}',1)[-1]
        if tag=='p':
            t=' '.join(x.text or '' for x in el.findall('.//w:t',NS)).strip()
            if t: paragraphs.append(t); last.append(t)
        elif tag=='tbl':
            rows=[]
            for row in el.findall('w:tr',NS):
                rows.append([' '.join(x.text or '' for x in c.findall('.//w:t',NS)).strip()
                             for c in row.findall('w:tc',NS)])
            tables.append({'preceding_text':last[-3:],'rows':rows})
    return {'paragraphs':paragraphs,'tables':tables}

def main(out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    result={'schema':'ALAN_original_2023_pollution_meta_source_intake_v1',
            'article_doi':'10.1002/ece3.9990','article_pmc':'PMC10111172',
            'is_original_experiment_data':False,'new_meta_model_fitted':False}
    xml=get(BASE+'fullTextXML',5_000_000)
    root=ET.fromstring(xml)
    result['article_xml_sha256']=hashlib.sha256(xml).hexdigest()
    result['article_publication_metadata']=[text(x) for x in root.findall('.//article-title')][:1]
    result['original_results_sections']=[text(s) for s in root.findall('.//sec')
        if s.find('title') is not None and any(k in text(s.find('title')).lower()
            for k in ('impact of pollution','study characteristics','heterogeneity','data extraction'))]
    result['references']=[{'id':r.attrib.get('id'),'text':text(r)} for r in root.findall('.//ref')]
    archive=get(BASE+'supplementaryFiles')
    if not zipfile.is_zipfile(io.BytesIO(archive)): raise ValueError('Original source is not a ZIP')
    result['supplement_archive_sha256']=hashlib.sha256(archive).hexdigest()
    result['supplements']=[]
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        result['members']=z.namelist()
        for name in z.namelist():
            if not name.lower().endswith('.docx'): continue
            info=z.getinfo(name)
            if info.file_size>20_000_000: raise ValueError('DOCX exceeds safety limit')
            raw=z.read(name)
            parsed=word_content(raw)
            item={'member':name,'sha256':hashlib.sha256(raw).hexdigest(),**parsed}
            result['supplements'].append(item)
    if not result['supplements']: raise ValueError('No original Word supplement found')
    target=out/'source_tables.json'
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('SOURCE_INTAKE_RECEIPT',json.dumps({k:v for k,v in result.items()
        if k not in ('supplements','references')},ensure_ascii=False))
    for item in result['supplements']:
        print('SUPPLEMENT_PARAGRAPHS',json.dumps(item['paragraphs'],ensure_ascii=False))
        for i,t in enumerate(item['tables']):
            print('SOURCE_TABLE',i,json.dumps(t,ensure_ascii=False))
    print('SOURCE_REFERENCES',json.dumps(result['references'],ensure_ascii=False))
    print('SOURCE_INTAKE_COMPLETED_NOT_A_NEW_EFFECT_ESTIMATE')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    main(p.parse_args().out)
