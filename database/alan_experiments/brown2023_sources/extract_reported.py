"""Preserve reported catch summaries and flag unresolved arithmetic differences."""
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[2].parent/'gbif_annual/alan_experiment_sources/brown2023/primary.xml'
tree=ET.parse(SOURCE)
arms=[]; checks=[]
for label in ['tbl1','tbl3','tbl4']:
    table=tree.find('.//table-wrap[@id="'+label+'"]')
    rows=[[''.join(c.itertext()).strip() for c in r] for r in table.findall('.//tbody/tr')]
    assert len(rows)==(6 if label=='tbl1' else 5)
    for row in rows:
        if label=='tbl1':
            arms.append(dict(table=label,endpoint='all_arthropods',period=row[0],arm=row[1],mean=float(row[2]),SE=float(row[3])))
        else:
            values=[]
            for period,cell in zip(['Pre','During','Post'],row[1:4]):
                mean,se=map(float,re.findall(r'\d+(?:\.\d+)?',cell))
                values.append(mean)
                arms.append(dict(table=label,endpoint=row[0],period=period,arm='Light',mean=mean,SE=se))
            for name,actual,printed in zip(['Pre-During','During-Post','Pre-Post'],
                    [values[0]-values[1],values[1]-values[2],values[0]-values[2]],row[4:]):
                estimate=float(printed.replace('−','-').split('(')[0])
                if abs(actual-estimate)>.02:
                    checks.append(dict(table=label,endpoint=row[0],comparison=name,
                        arithmetic_value=actual,reported_value=estimate,source_cell=printed,
                        status='clarify_model_adjustment_or_source_discrepancy_do_not_correct'))
assert len(arms)==36
for period in ['Pre','During','Post']:
    total=next(r['mean'] for r in arms if r['table']=='tbl1' and r['arm']=='Light' and r['period']==period)
    for table in ['tbl3','tbl4']:
        summed=sum(r['mean'] for r in arms if r['table']==table and r['period']==period)
        if (table=='tbl3' and abs(summed-total)>.04) or (table=='tbl4' and summed-total>.04):
            checks.append(dict(table=table,endpoint='sum_of_listed_groups',comparison=period,
                arithmetic_value=summed,reported_value=total,source_cell='Table1 all-arthropod mean',
                status='clarify_common_denominator_or_source_discrepancy_do_not_correct'))
for name,rows in [('reported_arm_summaries.csv',arms),('unresolved_arithmetic_checks.csv',checks)]:
    with (ROOT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit=dict(doi='10.1016/j.isci.2023.106203',source_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9982679/fullTextXML',
    sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),arm_rows=len(arms),unresolved_checks=len(checks),
    eligible_meta_effects_added=0,independent_unit='station; repeated dates not independent',
    comparison_scope='Table1 includes concurrent dark; Tables3and4 lit stations across periods only',
    SE_aggregation='not established as station-mean SE',raw_data='author_request_only_not_obtained')
(ROOT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
