"""Audit public pitfall biomass and actual treatment labels; no effect estimation."""
from pathlib import Path
import hashlib
import json
import pandas as pd
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CACHE = REPO.parent / 'gbif_annual/alan_experiment_sources/alanizer'
files = json.loads((CACHE/'zenodo_metadata.json').read_text(encoding='utf-8'))['files']
hashes = {}
for name in ['metadata.csv', 'biomass.csv']:
    raw = (CACHE/name).read_bytes()
    f = next(x for x in files if x['key'] == name)
    assert 'md5:'+hashlib.md5(raw).hexdigest() == f['checksum']
    hashes[name] = hashlib.sha256(raw).hexdigest()
m = pd.read_csv(CACHE/'metadata.csv')
b = pd.read_csv(CACHE/'biomass.csv')
keys = ['Site','StartDate']
assert not m.duplicated(keys).any() and not b.duplicated(keys).any()
j = b.merge(m, on=keys, how='left', suffixes=('_biomass','_metadata'), validate='one_to_one', indicator=True)
assert j._merge.eq('both').all()
assert j.Treatment_metadata.notna().all()
assert j.Treatment_biomass.eq(j.Treatment_metadata).all()
assert np.allclose(j.Total, j.Near+j.Far)
assert set(j.Treatment_metadata) == {'A','W','C'}
j['planned_actual_disagreement'] = j.Treatment_metadata.ne(j.PrescribedTreatment)
j.to_csv(HERE/'biomass_treatment_crosswalk.csv', index=False)
j[j.planned_actual_disagreement].to_csv(HERE/'corrected_treatment_rows.csv', index=False)
support = j.groupby(['Site','Treatment_metadata']).size().unstack(fill_value=0).reindex(columns=['A','C','W'],fill_value=0)
support.to_csv(HERE/'site_treatment_support.csv')
unmatched = m.merge(b[keys], on=keys, how='left', indicator=True)
audit = dict(
    source='https://doi.org/10.5281/zenodo.15170657', input_sha256=hashes,
    metadata_rows=len(m), metadata_sites=int(m.Site.nunique()), metadata_dates=int(m.StartDate.nunique()),
    metadata_actual_treatment_missing=int(m.Treatment.isna().sum()),
    biomass_rows=len(b), biomass_sites=int(b.Site.nunique()),
    biomass_start_date=str(b.StartDate.min()), biomass_last_start_date=str(b.StartDate.max()),
    all_biomass_rows_match_actual_treatment=True,
    planned_actual_disagreements_in_biomass=int(j.planned_actual_disagreement.sum()),
    metadata_rows_without_biomass=int(unmatched._merge.eq('left_only').sum()),
    metadata_actual_known_without_biomass=int((unmatched._merge.eq('left_only') & unmatched.Treatment.notna()).sum()),
    complete_three_arm_sites=int(support.gt(0).all(axis=1).sum()),
    raw_near_far_values=int(b[['Near','Far']].size),
    raw_total_mass_sum=float(b.Total.sum()),
    paper_reported_samples=216,
    paper_reported_mass_g=531.4,
    paper_reported_individuals=10211,
    exact_paper_sample_subset_verified=False,
    finite_effects_added=0,
    hold='Final indexed Methods/Results read; paper216samples differs from raw268trap values. Exact sample exclusions and raw mass-unit mapping unresolved; no subset chosen to force agreement')
(HERE/'paper_data_reconciliation.csv').write_text(
    'quantity,raw_file_value,paper_reported_value,status\n'
    f'near_far_sample_values,{b[["Near","Far"]].size},216,subset_not_verified\n'
    f'mass_sum,{b.Total.sum():.7f},531.4,raw_units_and_subset_not_verified\n'
    'individual_count,,10211,classification_files_not_read\n', encoding='utf-8')
(HERE/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
