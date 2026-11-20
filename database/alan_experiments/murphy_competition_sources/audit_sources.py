"""Inventory original files and subsamples; do not estimate unreplicated ALAN SE."""
from pathlib import Path
import hashlib
import json
import pandas as pd

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[2].parent / 'gbif_annual/alan_experiment_sources/murphy2026'
metadata = json.loads((CACHE / 'zenodo.json').read_text(encoding='utf-8'))
manifest = []
for item in metadata['files']:
    raw = (CACHE / item['key']).read_bytes()
    assert len(raw) == item['size']
    assert 'md5:' + hashlib.md5(raw).hexdigest() == item['checksum']
    manifest.append(dict(file=item['key'], url=item['links']['self'],
                         bytes=len(raw), md5_verified=True,
                         sha256=hashlib.sha256(raw).hexdigest()))
book = CACHE / 'Murphy et al. data 2025 ALAN competition.xlsx'
plants = pd.read_excel(book, sheet_name='Plant data')
cn = pd.read_excel(book, sheet_name='C_N data')
key = ['Species_TRT', 'Density', 'Rep', 'Trt']
for frame in (plants, cn):
    frame['Rep'] = frame['Rep'].str.lower()
assert not plants.duplicated(key).any()
chem_key = key + ['Status']
duplicates = cn[cn.duplicated(chem_key, keep=False)]
duplicates.to_csv(OUT / 'chemistry_duplicate_keys.csv', index=False)
assert set(plants.Trt) == {'L', 'NL'}
assert set(plants.Density) == {1, 2, 4}
counts = plants.groupby(['Species_TRT', 'Density', 'Trt']).size().rename('pot_rows')
assert len(counts) == 24 and counts.eq(10).all()
counts.to_csv(OUT / 'pot_cell_counts.csv')
missing = []
for label, frame, columns in [('plants', plants, ['Dry Mass Na', 'Dry Mass Ex']),
                              ('chemistry', cn, ['%N', '%C', 'C:N'])]:
    for arm, g in frame.groupby('Trt'):
        for col in columns:
            missing.append(dict(sheet=label, arm=arm, variable=col,
                                rows=len(g), missing=int(g[col].isna().sum())))
pd.DataFrame(missing).to_csv(OUT / 'missingness.csv', index=False)
merged = cn.merge(plants[key], on=key, how='left', indicator=True, validate='many_to_one')
assert merged['_merge'].eq('both').all()
audit = dict(plant_rows=len(plants), chemistry_rows=len(cn),
             unique_pot_keys=len(plants), treatment_cells=len(counts),
             chemistry_rows_matched_to_pots=len(merged),
             independent_light_rooms_per_arm=1,
             chemistry_distinct_pot_status_keys=len(cn.drop_duplicates(chem_key)),
             chemistry_rows_with_duplicated_keys=len(duplicates),
             chemistry_duplicated_key_groups=len(duplicates.drop_duplicates(chem_key)),
             missing_pot_status_keys_relative_to_two_per_pot=2*len(plants)-len(cn.drop_duplicates(chem_key)),
             missing_rows_are_not_absence=True,
             sampling_variance_for_replicated_light_effect_estimable=False,
             effects_added_to_main_catalogue=0,
             decision='retain_context_exclude_replicated_light_pooling')
(OUT / 'source_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
(OUT / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
print(json.dumps(audit, indent=2))
