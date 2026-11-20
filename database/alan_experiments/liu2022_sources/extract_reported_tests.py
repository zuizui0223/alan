"""Preserve author LR tests; these unsigned statistics are not meta effects."""
from pathlib import Path
import csv
import hashlib
import json
import re

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/liu2022'
text = (CACHE / 'primary.txt').read_text(encoding='utf-8')
table = text.split('TA B LE 1', 1)[1].split('Random effects', 1)[0]
pattern = r'^([^\n]+)\n1\n([0-9.]+)\n(\.[0-9]+)\s*$'
rows = [dict(term=term.strip(), df=1, chi_square=float(stat), reported_p=float(p),
             endpoint='sqrt_aboveground_pot_biomass', source_location='Table 1 p3679',
             status='reported_LR_test_not_signed_effect')
        for term, stat, p in re.findall(pattern, table, re.M)]
assert len(rows) == 15, len(rows)
assert any(r['chi_square'] == 0.009 and r['reported_p'] == 0.927 for r in rows)
with (OUT / 'reported_tests.csv').open('w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
audit = dict(doi='10.1111/gcb.16126', data_doi='10.5061/dryad.547d7wmb0',
             primary_sha256=hashlib.sha256((CACHE/'primary.pdf').read_bytes()).hexdigest(),
             text_sha256=hashlib.sha256((CACHE/'primary.txt').read_bytes()).hexdigest(),
             independent_cages=20, cages_per_factorial_cell=5, reported_tests=len(rows),
             effects_added=0, raw_download_status='Dryad files 1343897 and 1343899 returned HTTP401 on 2026-11-20',
             code_archive='Zenodo6007775 contains analysis Rmd, not raw CSV',
             coordinates='Printed sexagesimal minutes/seconds invalid; not geocoded',
             pot_count='Methods total660 conflicts with 13x35+7x34=693; unresolved')
(OUT/'source_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
print('15 original LR tests retained; no signed effect or sampling variance inferred')
