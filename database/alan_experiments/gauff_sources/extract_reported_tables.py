"""Preserve reported comparisons without repairing inconsistent source labels."""
from pathlib import Path
import hashlib
import json
import pandas as pd
from docx import Document

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/gauff2026'
doc=Document(CACHE/'supplement.docx')
assert len(doc.tables)==14
tables=[[[cell.text.strip() for cell in row.cells] for row in table.rows] for table in doc.tables]
rows=[]
for number,endpoint in [(8,'PAM_yield_assay_stage_label_conflict'),(11,'leaf_length')]:
    table=tables[number-1]
    assert table[0]==['Comparison','Estimate','Std.Error','z.value','p.value']
    for cells in table[1:]:
        assert len(cells)==5
        rows.append(dict(source_table=number,endpoint=endpoint,comparison_as_reported=cells[0],
                         estimate=float(cells[1]),reported_standard_error=float(cells[2]),
                         reported_statistic=float(cells[3]),reported_tukey_p=float(cells[4]),
                         source_location='Archimer_supplement_docx',
                         use_for_meta=False,hold='source_labels_and_raw_data_unresolved'))
pd.DataFrame(rows).to_csv(OUT/'reported_pairwise_comparisons.csv',index=False)
assert tables[1][1:]==tables[2][1:]
manifest=[]
for filename,url in [('primary.pdf','https://archimer.ifremer.fr/doc/01074/118553/132576.pdf'),
                     ('supplement.docx','https://archimer.ifremer.fr/doc/01074/118553/132577.docx')]:
    b=(CACHE/filename).read_bytes()
    manifest.append(dict(file=filename,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
audit=dict(primary_version='institutional journal pre-proof',doi='10.1016/j.marenvres.2026.108273',
           lighting_plots=6,lit_plots=3,dark_plots=3,quadrats_per_plot=6,
           cleaned_quadrats_per_plot=3,epibiota_removal='every two weeks in field',
           allocation='alternating lighting plots; random assignment not established',
           raw_unit_data_acquired=False,reported_pairwise_rows=len(rows),
           oxygen_tables_2_and_3_identical_below_header=True,
           leaf_length_posthoc_conflict=True,PAM_stage_label_conflict=True,
           independent_studies_added_to_pooled_analysis=0,main_catalogue_effects_added=0,
           note='Signed oxygen fluxes cannot be treated as ordinary positive lnRR outcomes; no SE reconstructed from rounded t values.')
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
