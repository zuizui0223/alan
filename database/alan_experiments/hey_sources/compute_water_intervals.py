"""Uncertainty for existing water moderation contrasts; no catalogue changes.

Each arm contains five independent lamp plots. Within-plot water covariance
is retained separately for each density; incomplete pot remains flagged.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
src = ROOT / 'hey2020_biomass_units.csv'
units = pd.read_csv(src)
previous = pd.read_csv(ROOT / 'hey2020_water_moderation.csv').set_index('density')
influences = {}
means = {}
for arm in ['Lit', 'Control']:
    a = units[units.Light.eq(arm)].pivot(index='Plot', columns=['Density','Water'], values='biomass_g')
    assert a.shape == (5, 4) and a.notna().all().all()
    means[arm] = a.mean()
    influences[arm] = pd.DataFrame({density:
        a[(density, 'Water')] / a[(density, 'Water')].mean()
        - a[(density, 'NW')] / a[(density, 'NW')].mean()
        for density in ['Low','High']})
covs = {arm: x.cov()/len(x) for arm,x in influences.items()}
cov = covs['Lit'] + covs['Control']
rows = []
for density in ['Low','High']:
    estimate = sum(sign*np.log(means[arm][(density,'Water')]/means[arm][(density,'NW')])
                   for arm,sign in [('Lit',1),('Control',-1)])
    parts = np.array([covs[arm].loc[density,density] for arm in ['Lit','Control']])
    variance = parts.sum()
    df = variance**2 / np.sum(parts**2/4)
    half = t.ppf(.975, df)*np.sqrt(variance)
    assert np.isclose(estimate, previous.loc[density,'estimate'])
    assert np.isclose(variance, previous.loc[density,'variance'])
    rows.append(dict(density=density,ln_ratio_of_responses=estimate,variance=variance,
                     welch_df=df,ratio_of_responses=np.exp(estimate),
                     ci_low=np.exp(estimate-half),ci_high=np.exp(estimate+half),
                     independent_light_plots=5,independent_control_plots=5,
                     incomplete_pot_in_contrast=(density=='High'),
                     status='conditional_delta_welch_single_experiment_not_meta_pool'))
pd.DataFrame(rows).to_csv(ROOT/'water_moderation_intervals.csv',index=False)
cov.to_csv(ROOT/'water_moderation_covariance.csv',index_label='density')
assert np.linalg.eigvalsh(cov.to_numpy()).min() >= -1e-12
(ROOT/'water_interval_audit.json').write_text(json.dumps(dict(
    input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
    n_plots_by_arm={k:len(v) for k,v in influences.items()},
    matched_existing_estimates_and_variances=True,
    missing_pot_handling='existing available-individual mean retained, not zero-filled',
    interpretation='water-addition modification of light biomass response; not habitat wetness',
    model='first-order delta variance with Welch-Satterthwaite t intervals',
    final_eligibility=False),indent=2)+'\n',encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False))
