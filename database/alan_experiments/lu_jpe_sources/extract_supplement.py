"""Preserve publisher supplement tables; do not infer lnRR or independent n."""
import csv
import hashlib
import json
from pathlib import Path
from docx import Document

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parents[2].parent / 'gbif_annual/alan_experiment_sources/lu_jpe2026/supplement.docx'
doc = Document(SOURCE)
specs = [(2, 'S2', 18), (3, 'S3', 5), (4, 'S4a', 12), (5, 'S4b', 12)]
counts = {}
for index, label, count in specs:
    table = doc.tables[index]
    header = [c.text.strip() for c in table.rows[0].cells]
    rows = [dict(zip(header, [c.text.strip() for c in r.cells])) for r in table.rows[1:]]
    assert len(rows) == count
    for row in rows:
        row['source_table'] = label
        row['doi'] = '10.1093/jpe/rtag160'
        row['status'] = 'reported_only_main_methods_and_dependence_unverified'
    with (ROOT / f'table_{label}.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    counts[label] = len(rows)
audit = dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    source_bytes=SOURCE.stat().st_size, table_rows=counts,
    landing_page='https://academic.oup.com/jpe/advance-article/doi/10.1093/jpe/rtag160/8731923?searchresult=1',
    publisher_supplement_path='https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/jpe/PAP/10.1093_jpe_rtag160/1/rtag160_supplementary_data.docx',
    acquisition='Publisher-issued signed supplementary link returned200 on2026-11-21; unsigned path403. Main PDF not recovered. Web reader unsupported DOCX type; downloaded original signed link.',
    effect_catalogue_rows_added=0, full_primary_article_read=False)
(ROOT / 'supplement_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
print(json.dumps(audit, indent=2))
