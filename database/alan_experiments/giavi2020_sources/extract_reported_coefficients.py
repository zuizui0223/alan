"""Transcribe published field-model coefficients; never convert these to lnRR."""
from pathlib import Path
import csv
import hashlib
import json
from lxml import etree

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources'
SOURCE = CACHE / 'fulltexts_20261120/PMC7368033.xml'
root = etree.parse(str(SOURCE))
table = root.xpath('//table-wrap')[0]
spec = {'(a)': ('Epilobium angustifolium', 'selected_capsule_at_least_four_seeds', 'Dark'),
        '(b)': ('Silene latifolia', 'capsule_infection', 'Adjacent'),
        '(c)': ('Silene latifolia', 'relative_reproductive_rate_denominator_hold', 'Adjacent')}
rows = []
section = None
for tr in table.xpath('.//tr'):
    values = [' '.join(x.itertext()).strip() for x in tr]
    if len(values) == 1 and values[0] in spec:
        section = values[0]
    elif section and values[0] in {'Adjacent', 'Illuminated', 'Dark'}:
        species, endpoint, reference = spec[section]
        number = lambda s: float(s.replace('−', '-').replace('\xa0','').replace(' ',''))
        rows.append(dict(species=species, endpoint=endpoint, contrast=f'{values[0]}_minus_{reference}',
                         estimate=number(values[1]), se=number(values[2]), reported_p=values[4],
                         source_location=f'Table1{section}', scale='original_binomial_model_coefficient_link_not_explicit',
                         sampling_covariance='not_recovered', pooling_status='hold_not_lnRR'))
assert len(rows) == 6
text = ' '.join(root.getroot().itertext())
assert '0.76' in text and '0.521' in text and '0.145' in text
rows.append(dict(species='Silene latifolia', endpoint='capsule_infection', contrast='Illuminated_minus_Dark',
                 estimate=.76, se=.521, reported_p='0.145', source_location='Results: seed predation paragraph',
                 scale='original_binomial_model_coefficient_link_not_explicit', sampling_covariance='not_recovered',
                 pooling_status='hold_not_lnRR'))
with (OUT/'reported_field_coefficients.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
supp = CACHE/'giavi2020/41598_2020_68667_MOESM1_ESM.pdf'
assert hashlib.md5(supp.read_bytes()).hexdigest() == '3a80b1212188fc62b3c59ff9e9cac9d1'
audit=dict(article_doi='10.1038/s41598-020-68667-y', article_xml_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
           supplement_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7368033/supplementaryFiles',
           supplement_sha256=hashlib.sha256(supp.read_bytes()).hexdigest(), supplement_content='Matlab video tracking code, not field raw observations or fitted model objects',
           coefficients=7, outcomes=3, independent_meadows=8, lit_meadows=4, dark_meadows=4,
           adjacent_sectors_share_lit_meadows=True, raw_observations_recovered=False, model_rerun=False,
           limitations=['Link not explicit in paper; no default link silently assigned',
                        'Relative reproduction denominator description ambiguous',
                        'Site random effect may refer to sectors rather than shared meadows; unresolved',
                        'Contrasts share controls and outcomes; missing covariance not zero',
                        'Relationship to other Swiss meadow programmes unresolved'])
(OUT/'source_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print('Extracted seven reported coefficients; no lnRR or pooled estimate added')
