"""Extract original reported coefficients, keeping model scales and endpoints separate."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from xml.etree import ElementTree as ET

p = argparse.ArgumentParser()
p.add_argument('primary_xml', type=Path)
args = p.parse_args()
out = Path(__file__).resolve().parent
raw = args.primary_xml.read_bytes()
assert hashlib.sha256(raw).hexdigest() == 'ffc015c8f40ed9063b13739c7148c438fc5d6776df4aabea7449d9137447f6e9'
root = ET.fromstring(raw)
rows = []
for table in root.findall('.//table-wrap'):
    response = None
    for tr in table.findall('.//tr')[1:]:
        cells = [' '.join(''.join(c.itertext()).split()) for c in tr]
        if len(cells) == 6:
            response = cells.pop(0)
        assert len(cells) == 5 and response
        term, estimate, se, chisq, p_value = cells
        rows.append(dict(source_table=table.get('id'), response=response, term=term,
                         estimate_original=estimate, se_original=se,
                         likelihood_ratio_chisq_original=chisq, p_original=p_value))
assert len(rows) == 28
def write(name, records):
    with (out / name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)
write('reported_coefficients.csv', rows)
# Only these three responses explicitly use a binomial logit model in Methods.
logit = {'Breeding site presence', 'Calling presence', 'Calling effort'}
effects = []
for r in rows:
    if r['term'] == 'Light' and r['response'] in logit:
        beta, se = float(r['estimate_original']), float(r['se_original'])
        effects.append(dict(response=r['response'], log_odds_ratio=beta, reported_se=se,
                            odds_ratio=math.exp(beta),
                            approximate_wald95_low=math.exp(beta-1.96*se),
                            approximate_wald95_high=math.exp(beta+1.96*se),
                            source_table=r['source_table'],
                            model='author_reduced_model_without_light_noise_interaction',
                            status='single_study_reported_model_not_pooled'))
assert len(effects) == 3
write('reported_light_odds_ratios.csv', effects)
(out / 'audit.json').write_text(json.dumps(dict(
    article_doi='10.1093/beheco/arac077', source_sha256=hashlib.sha256(raw).hexdigest(),
    dataset_doi='10.5061/dryad.4j0zpc8fj', dataset_version=187596,
    original_coefficients=len(rows), derived_logit_contrasts=len(effects),
    covariance_available=False, raw_data_acquired=False, new_pooled_effects=0,
    access_check='2026-11-21: Dryad metadata HTTP 200; activity file API download HTTP 401; not bypassed'
), indent=2) + '\n', encoding='utf-8')
print(json.dumps(effects, indent=2))
