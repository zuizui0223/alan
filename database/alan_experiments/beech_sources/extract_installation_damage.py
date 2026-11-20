"""Supplementary paired-installation damage contrasts; two sites, shared controls."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
RAW = OUT.parents[3] / 'gbif_annual/alan_experiment_sources/beech2025/raw_data_ceppi.xlsx'
assert hashlib.sha256(RAW.read_bytes()).hexdigest() == 'b7c0c7a81315db06951cabcfcb070f926236a42d1f05c0b93b65333cad88e78d'
leaves = pd.read_excel(RAW, sheet_name='Leaf_area_consumed')
label = 'Treatment (color_intensity)'
doses = {'4000_30': 21, '4000_50': 35, '4000_100': 70}
selected = leaves[leaves[label].isin(['control', *doses])].copy()
assert len(selected) == 840 and selected.Perc_Herbivory.between(0, 100).all()
rows = []
for (site, installation, treatment), g in selected.groupby(['Site', 'Site_ID', label]):
    damaged = g.Perc_Herbivory.gt(0)
    assert len(g) == 105 and g.Tree_ID.nunique() == 7 and damaged.any()
    rows.append(dict(site=site, installation=installation, treatment=treatment,
                     trees=7, leaves=len(g), damaged_leaves=int(damaged.sum()),
                     damage_incidence=damaged.mean(),
                     all_leaf_mean_percent=g.Perc_Herbivory.mean(),
                     damaged_leaf_mean_percent=g.loc[damaged, 'Perc_Herbivory'].mean()))
units = pd.DataFrame(rows)
units.to_csv(OUT / 'installation_damage.csv', index=False)
endpoints = ['damage_incidence', 'all_leaf_mean_percent', 'damaged_leaf_mean_percent']
effects, scores = [], []
site_order = sorted(units.site.unique())
assert len(site_order) == 2
for endpoint in endpoints:
    wide = units.pivot(index='site', columns='treatment', values=endpoint).loc[site_order]
    assert wide.notna().all().all() and wide.gt(0).all().all()
    dark = wide['control'].to_numpy()
    for treatment, lux in doses.items():
        lit = wide[treatment].to_numpy()
        score = lit/lit.mean() - dark/dark.mean()
        var = np.var(score, ddof=1)/2
        log_rr = np.log(lit.mean()/dark.mean())
        ci = np.exp(log_rr + np.array([-1, 1])*t.ppf(.975, 1)*np.sqrt(var))
        effect_id = f'Ceppi2025_{treatment}_{endpoint}'
        effects.append(dict(effect_id=effect_id, endpoint=endpoint, treatment=treatment,
                            nominal_lux=lux, independent_site_pairs=2,
                            mean_lit=lit.mean(), mean_dark=dark.mean(),
                            lnRR=log_rr, sampling_variance=var, ratio=np.exp(log_rr),
                            lower=ci[0], upper=ci[1], df=1,
                            interval='paired_delta_t1_approximation',
                            status='supplementary_programme_overlap_pending'))
        scores.append(score)
effects = pd.DataFrame(effects)
score_matrix = np.array(scores).T
covariance = np.cov(score_matrix, rowvar=False, ddof=1)/2
assert np.allclose(np.diag(covariance), effects.sampling_variance)
assert np.linalg.eigvalsh(covariance).min() > -1e-10
effects.to_csv(OUT / 'installation_effects.csv', index=False)
pd.DataFrame(covariance, index=effects.effect_id, columns=effects.effect_id).to_csv(OUT / 'installation_covariance.csv')
pd.DataFrame(score_matrix, index=site_order, columns=effects.effect_id).to_csv(OUT / 'site_scores.csv')
audit = dict(source_leaf_rows=len(leaves), selected_leaf_rows=len(selected),
             extra_spectral_leaf_rows=len(leaves)-len(selected), installation_rows=len(units),
             sites=2, contrasts=len(effects), covariance_rank=int(np.linalg.matrix_rank(covariance)),
             control_shared=True, endpoints_share_leaves=True,
             leaf_or_tree_replication_used_for_precision=False,
             original_GLMM_reproduced=False, main_catalogue_additions=0,
             conditional_endpoint='mean leaf percent among damaged leaves within each installation',
             all_leaf_endpoint='unweighted mean leaf percent including zeros; not total-area-weighted loss',
             limitations=['Two site pairs only; t1 delta intervals are approximate.',
                          'No general dose-response threshold or causal chemical mediation estimated.',
                          'Shared Swiss programme overlap unresolved.'])
(OUT/'installation_extraction_audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
print(effects[['endpoint','nominal_lux','ratio','lower','upper']].to_string(index=False))
