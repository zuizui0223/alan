"""Original data audit; row-level descriptive effects, outside final catalogue."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/cieraad2023'
src = CACHE / 'Cieraadetal2022_JApplEcol_data.xlsx'
raw = src.read_bytes()
assert hashlib.sha1(raw).hexdigest() == '52153a23c24285223a2c798e8592f7da836b9c4c'
data = pd.concat([pd.read_excel(src, sheet_name='Cieraadetal2022_' + s)
                  for s in ['Blueberry', 'Oak']], ignore_index=True)
assert data.PlantID.is_unique and len(data) == 317
assert data.groupby('Transect')[['Site', 'Treatment']].nunique().eq(1).all().all()
assert data.groupby(['Site', 'Treatment']).Transect.nunique().eq(1).all()
assert data.Herbiv_SumFrac.between(0, 1).all()
counts = data.groupby(['Species', 'Site', 'Treatment', 'Transect', 'TransectPost']).size().rename('measured_plants')
counts.to_csv(OUT / 'source_cell_counts.csv')
# Equal posts within each illumination row; leaves and plants are subsamples.
post = data.groupby(['Species', 'Site', 'Treatment', 'Transect', 'TransectPost']).Herbiv_SumFrac.mean()
row = post.groupby(['Species', 'Site', 'Treatment', 'Transect']).mean().rename('fraction_leaf_area_lost')
row.to_csv(OUT / 'row_means.csv')
frame = row.reset_index()
blocks = frame.pivot(index='Site', columns=['Species', 'Treatment'], values='fraction_leaf_area_lost')
assert blocks.shape == (2, 8)
means = blocks.mean()
cov = blocks.cov().to_numpy() / 2
effects, jac, within = [], [], []
for species in ['Blueberry', 'Oak']:
    for colour in ['Green', 'Red', 'White']:
        a, b = (species, colour), (species, 'Dark')
        j = np.zeros(len(means))
        j[means.index.get_loc(a)] = 1 / means[a]
        j[means.index.get_loc(b)] = -1 / means[b]
        estimate = np.log(means[a] / means[b])
        variance = float(j @ cov @ j)
        half = t.ppf(.975, 1) * np.sqrt(variance)
        effects.append(dict(effect_id=f'Cieraad2023_{species}_{colour}_herbivory',
            species=species, treatment=colour, mean_lit=means[a], mean_dark=means[b],
            lnRR=estimate, sampling_variance=variance, lower=estimate-half, upper=estimate+half,
            n_independent_rows_per_arm=2, n_forest_blocks=2,
            status='descriptive_two_blocks_programme_overlap_hold'))
        jac.append(j)
        for site in blocks.index:
            within.append(dict(species=species, treatment=colour, forest=site,
                response_ratio=blocks.loc[site, a]/blocks.loc[site, b],
                independent_rows_per_arm=1, precision='not_estimable_within_forest'))
e = pd.DataFrame(effects)
v = np.asarray(jac) @ cov @ np.asarray(jac).T
assert np.allclose(np.diag(v), e.sampling_variance)
assert np.linalg.eigvalsh(v).min() > -1e-10
e.to_csv(OUT / 'descriptive_effects.csv', index=False)
pd.DataFrame(within).to_csv(OUT / 'forest_specific_descriptive_ratios.csv', index=False)
pd.DataFrame(v, index=e.effect_id, columns=e.effect_id).to_csv(OUT / 'effect_covariance.csv', index_label='effect_id')
audit = dict(article_doi='10.1111/1365-2664.14336', data_doi='10.34894/VD98DG',
    file_id=320748, sha1=hashlib.sha1(raw).hexdigest(), sha256=hashlib.sha256(raw).hexdigest(),
    measured_plants=data.Species.value_counts().to_dict(), independent_rows=8, posts=16,
    independent_rows_per_arm=2, within_forest_rows_per_arm=1, effect_rows=6,
    covariance_rank=int(np.linalg.matrix_rank(v)), catalogue_included=False,
    aggregation='plant means within posts; equal post means within rows; equal forest-block means',
    cautions=['Conditional on measured plants; not a mortality or population endpoint.',
              'Two forest blocks only; approximate t1 intervals provide very weak precision.',
              'Forest-specific uncertainty cannot be estimated from one row per arm.',
              'Shared dark rows and plant species induce covariance; six rows are not six experiments.',
              'Light on Nature programme overlap must be adjudicated before pooling.',
              'No measured predator, pollinator or reproductive response in this extraction.'])
(OUT / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
print(e[['species', 'treatment', 'lnRR', 'sampling_variance']].to_string(index=False))
