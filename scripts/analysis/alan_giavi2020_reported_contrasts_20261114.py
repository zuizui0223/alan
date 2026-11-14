#!/usr/bin/env python3
"""Transform original published Giavi 2020 logit coefficients, not a data refit."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
from scripts.analysis.alan_crossstudy_scope_20261114 import odds_ratio, require


def number(value):
    return float(value.replace('\u2212','-').replace(' ',''))


def sections(rows):
    out={};current=None
    for row in rows:
        if len(row)==1 and row[0] in ('(a)','(b)','(c)'):
            current=row[0];out[current]={}
        elif current and len(row)==5:
            require(row[0] not in out[current],'Duplicate coefficient label')
            out[current][row[0]]={'beta':number(row[1]),'se':number(row[2]),'reported_p':row[4]}
    require(set(out)=={'(a)','(b)','(c)'},'Source model sections missing')
    return out


def analyze(spillover):
    tables=[t for t in spillover['original_published_model_tables'] if t['label']=='Table 1']
    require(len(tables)==1,'Original Table 1 missing')
    block=sections(tables[0]['rows'])
    require('Adjacent (intercept)' in block['(b)'],'Seed-predation reference changed')
    require('Adjacent (intercept)' in block['(c)'],'Reproductive reference changed')
    cases=[]
    for reference in ('Dark','Illuminated'):
        original=block['(b)'][reference]
        cases.append({'outcome':'capsule infestation odds conditional on developed capsules',
            'contrast':'Adjacent / '+reference,
            'source_table':'1(b)','original_reference':'Adjacent',
            'source_logit_coefficient':original['beta'],
            'source_p':original['reported_p'],
            **odds_ratio(-original['beta'],original['se'])})
    b=block['(c)']['Dark']
    result={'schema':'ALAN_Giavi2020_reported_model_spatial_spillover_contrasts_v1',
        'source_doi':spillover['article_doi'],'source_xml_sha256':spillover['article_xml_sha256'],
        'published_predation_contrasts':cases,
        'relative_reproduction_published_coefficient':{
            'contrast':'Adjacent minus Dark','beta':-b['beta'],'se':b['se'],
            'reported_p':b['reported_p'],'source_table':'1(c)',
            'raw_model_refit':False,
            'denominator_semantics_need_author_code':True,
            'issue':'Methods describe cbind(unattacked capsules, total exposed flowers), rather than explicit failures; do not silently claim a success/flower odds ratio or percent decline.'},
        'original_replication':{'meadows':8,'lit_meadows':4,'dark_control_meadows':4,
            'illuminated_and_adjacent_share_same_four_meadows':True,
            'potted_plants':97,'Silene_plants':47,'Epilobium_plants':50},
        'interpretation_limits':[
            'Published Table1 coefficient transformation only, not a refitted experiment or new independent observation.',
            'Odds ratios are not risk ratios, percentage seed loss, or measured predator density.',
            'Adjacent dark locations can experience treatment spillover; they are not untreated independent controls.',
            'Relative-reproduction result has published denominator wording ambiguity; row-level data and code are needed to resolve it.',
            'Predation is by pollinator offspring, not adult moth predators or predator-derived pollination.',
            'This study does not measure island-mainland contrasts or establish satellite dark connectivity as effective refuge.'
        ]}
    return result


def main(source,out):
    data=json.loads(Path(source).read_text())
    r=analyze(data['spillover'])
    target=Path(out)
    require(not target.exists(),'Refuse overwrite')
    target.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n')
    print('GIAVI_REPORTED_CONTRAST_TRANSFORM',json.dumps(r,ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();main(a.source,a.out)
