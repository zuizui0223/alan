#!/usr/bin/env python3
"""Source-derived ALAN-specific evidence scope and published spillover model intake.

This does not refit individual studies or a meta-analysis. Paper membership is
not inverse-variance weight, and author-year labels are not independent trials.
"""
from __future__ import annotations
import argparse, hashlib, json, math, re, unicodedata
from pathlib import Path
from xml.etree import ElementTree as ET
from scripts.analysis.alan_crossstudy_source_intake_20261114 import get, text


def parse_references(value: str) -> set[str]:
    result=set()
    value=unicodedata.normalize('NFKC',value).strip().strip('()')
    for piece in value.split(';'):
        found=list(re.finditer(r'\b((?:19|20)\d{2}[a-z]?)\b',piece))
        if not found: raise ValueError('Unresolved source citation: '+piece)
        author=piece[:found[0].start()].strip(' ,')
        author=re.split(r'\s+et\s+al\.?|\s+and\s+',author,flags=re.I)[0]
        author=' '.join(author.casefold().split()).replace('\u2010','-').replace('\u2011','-')
        if not author: raise ValueError('Missing author')
        for match in found:
            year=match.group(1)
            # These two labels share DOI 10.1098/rsbl.2016.0874 in the source.
            # Do NOT collapse arbitrary author/year a/b distinctions.
            if author=='van langevelde' and year in {'2017a','2017b'}: year='2017'
            result.add(author+':'+year)
    return result


def classify_scope(light_cells: list[str], uv_text: str, reproduction_text: str) -> dict:
    light=set().union(*(parse_references(x) for x in light_cells))
    uv=parse_references(uv_text)
    reproduction=parse_references(reproduction_text)
    if not uv<=light or not reproduction<=light:
        raise ValueError('Scope source mismatch; missing exposure cannot be inferred')
    return {'light_papers':sorted(light),'uv_papers':sorted(uv),
            'alan_candidate_papers':sorted(light-uv),
            'reproductive_light_papers':sorted(reproduction),
            'excluded_uv_reproductive_papers':sorted(reproduction&uv),
            'alan_candidate_reproductive_papers':sorted(reproduction-uv),
            'alan_only_effect_reestimated':False,
            'citation_labels_not_independent_experiment_certification':True,
            'membership_counts_are_not_effect_size_weights':True}


def odds_ratio(beta: float,se: float) -> dict:
    if not math.isfinite(beta) or not math.isfinite(se) or se<0:
        raise ValueError('Invalid estimate or uncertainty')
    return {'log_odds':beta,'se':se,'odds_ratio':math.exp(beta),
            'ci95':[math.exp(beta-1.959963984540054*se),math.exp(beta+1.959963984540054*se)],
            'ci_method':'normal Wald from rounded published coefficient and SE',
            'new_raw_data_model':False,'not_risk_ratio_or_percent_change':True}


def require(condition, message):
    if not condition: raise ValueError(message)


def find_table(tables, appendix):
    hits=[t for t in tables if t.get('preceding_text') and
          re.search(r'Appendix\s+S'+str(appendix)+r'\s*\.',t['preceding_text'][-1])]
    require(len(hits)==1,'Ambiguous source appendix '+str(appendix))
    return hits[0]['rows']


def scope_from_original(source):
    require(source['article_doi']=='10.1002/ece3.9990','Wrong source article')
    require(len(source['supplements'])==1,'Unexpected supplement count')
    tables=source['supplements'][0]['tables']
    s4=find_table(tables,4); s7=find_table(tables,7); s13=find_table(tables,13)
    require(s4[0][2]=='Light','Original source light column moved')
    reproduction=[r[2] for r in s4[1:] if 'reproductive' in r[0].casefold()]
    require(len(reproduction)==1,'Missing reproductive source classification')
    uv=[r[2] for r in s7[1:] if r[1]=='UV-B']
    require(len(uv)==1,'Missing source UV-B classification')
    r=classify_scope([row[2] for row in s4[1:]],uv[0],reproduction[0])
    non_uv=set(); category=None
    for row in s7[1:]:
        if row[0]: category=row[0]
        if category=='Light pollution' and row[1]!='UV-B':
            non_uv.update(parse_references(row[2]))
    require(non_uv==set(r['alan_candidate_papers']),
            'Some non-UV source papers lack classified exposure; fail closed')
    r['non_uv_source_exposure_labels_verified']=True
    r['counts']={k:len(v) for k,v in r.items() if isinstance(v,list)}
    require(r['counts']=={'light_papers':22,'uv_papers':8,'alan_candidate_papers':14,
        'reproductive_light_papers':13,'excluded_uv_reproductive_papers':7,
        'alan_candidate_reproductive_papers':6},'Original 2023 source membership changed')
    require(len(s13)==2 and s13[1][0]=='(Intercept)','Source model table changed')
    model=s13[1]
    r['original_broad_light_meta_model']={
        'hedges_g':float(model[1]),'se':float(model[2]),'z':float(model[3]),
        'p':float(model[4]),'CI95':[float(x) for x in model[5].split(',')],
        'includes_UVB_and_different_outcome_levels':True,
        'not_an_ALAN_only_reproduction_estimate':True}
    r['source_supplement_sha256']=source['supplements'][0]['sha256']
    r['macgregor2019_in_source_light_membership']='macgregor:2019' in r['light_papers']
    r['cannot_infer_why_any_paper_is_unlisted']=True
    r['effect_level_yi_vi_table_available']=False
    r['scope_only_no_meta_coefficient_subtraction']=True
    return r


def acquire_giavi_tables():
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7368033/fullTextXML'
    raw=get(url,5_000_000); root=ET.fromstring(raw)
    dois=[text(a) for a in root.findall('.//article-id') if a.attrib.get('pub-id-type')=='doi']
    require('10.1038/s41598-020-68667-y' in dois,'Wrong spillover article')
    tables=[]
    for table in root.findall('.//table-wrap'):
        rows=[]
        for row in table.findall('.//tr'):
            rows.append([text(c) for c in row if c.tag in ('th','td')])
        tables.append({'id':table.attrib.get('id'),'label':text(table.find('label'))
                      if table.find('label') is not None else None,
                      'caption':text(table.find('caption')) if table.find('caption') is not None else None,
                      'rows':rows})
    sections=[]
    for section in root.findall('.//sec'):
        title=section.find('title')
        if title is not None and text(title).casefold() in {
            'data availability','statistical analysis','field experiment',
            'estimation of reproductive output and pre-dispersal seed predation'}:
            sections.append({'title':text(title),'text':text(section)})
    return {'article_doi':dois[0],'PMCID':'PMC7368033','article_xml_sha256':hashlib.sha256(raw).hexdigest(),
            'original_published_model_tables':tables,'source_method_sections':sections,
            'raw_observations_acquired':False,'new_model_fitted':False,
            'statistical_unit_warning':'97 plants are nested in 8 meadows; illuminated and adjacent share 4 meadows'}


def main(source_path,out):
    source=json.loads(Path(source_path).read_text())
    result={'schema':'ALAN_source_scope_and_spillover_published_table_v1',
            'review_scope':scope_from_original(source),
            'spillover':acquire_giavi_tables()}
    target=Path(out); target.parent.mkdir(parents=True,exist_ok=True)
    require(not target.exists(),'Refuse overwrite')
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('ALAN_ORIGINAL_SCOPE_RESULT',json.dumps(result['review_scope'],ensure_ascii=False))
    print('ALAN_GIAVI_ORIGINAL_MODEL_TABLES',json.dumps(result['spillover'],ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();main(a.source,a.out)
