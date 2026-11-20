"""Extract reported arm summaries; no unobserved paired variance is imputed."""
from pathlib import Path
from lxml import etree
import re,csv,json,hashlib

OUT=Path(__file__).resolve().parent
SOURCE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/boom2020/primary.xml'
tree=etree.parse(str(SOURCE))
endpoints=['flower_fruit_set','fruit_infestation','fertilized_ovules_uninfested_fruit','estimated_seeds_per_flower']
rows=[]
for tr in tree.xpath('//table-wrap//tbody/tr'):
    cells=[''.join(td.itertext()) for td in tr]
    for endpoint,cell in zip(endpoints,cells[1:]):
        m=re.fullmatch(r'([\d.]+)\s*±\s*([\d.]+)(?:\s*\(([\d,]+)\))?',cell)
        assert m,cell
        rows.append(dict(treatment=cells[0],endpoint=endpoint,mean=float(m[1]),SE=float(m[2]),
                         reported_organ_count=int(m[3].replace(',','')) if m[3] else '',
                         site_blocks=4,source_location='Table1',status='paired_covariance_and_SE_aggregation_hold'))
assert len(rows)==16
with (OUT/'reported_arm_summaries.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit=dict(doi='10.1038/s41598-020-75471-1',pmcid='PMC7591485',
           source_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7591485/fullTextXML',
           xml_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
           rows=16,meta_effects_added=0,independent_units='4sites x4randomized treatment transects',
           programme='Light on Nature; cross-report/site overlap unresolved',
           variance_hold='Within-site light-dark covariance absent; do not use organ counts as independent n',
           reproductive_output='Product of three proportions and average ovule count, not observed seeds per flower')
(OUT/'source_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print('16 arm summaries extracted; no effect sampling variance inferred')
