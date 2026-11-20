"""Small-sample intervals for existing zone-level effects, not new experiments.

Welch-Satterthwaite degrees of freedom on delta-method log-mean variances.
Conditional on surviving cylinders; four dependent species, no pooling.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t

HERE = Path(__file__).resolve().parent
units = pd.read_csv(HERE / 'monoculture_zone_units.csv')
effects = pd.read_csv(HERE / 'monoculture_effects.csv')
rows = []
for _, effect in effects.iterrows():
    species = effect.species.replace(' ', '_')
    groups = [units.loc[(units.Species == species) & (units.Light == arm), 'biomass_g'].to_numpy()
              for arm in ['ALAN', 'No ALAN']]
    assert all(len(g) == 3 and np.all(g > 0) for g in groups)
    terms = [g.var(ddof=1) / len(g) / g.mean()**2 for g in groups]
    variance = sum(terms)
    df = variance**2 / sum(v**2 / (len(g) - 1) for v, g in zip(terms, groups))
    estimate = np.log(groups[0].mean() / groups[1].mean())
    assert np.isclose(estimate, effect.lnRR, atol=1e-12)
    assert np.isclose(variance, effect.variance, atol=1e-12)
    half = t.ppf(.975, df) * np.sqrt(variance)
    rows.append(dict(effect_id=effect.effect_id, species=effect.species,
                     response_ratio=np.exp(estimate), ci_lower=np.exp(estimate-half),
                     ci_upper=np.exp(estimate+half), welch_df=df, lnRR=estimate,
                     variance=variance, n_light_zones=3, n_dark_zones=3,
                     conditional_on_surviving_cylinders=True,
                     interval='delta_log_mean_Welch_t_unadjusted_95_percent'))
out = pd.DataFrame(rows)
out.to_csv(HERE / 'monoculture_intervals.csv', index=False)
print(out[['species', 'response_ratio', 'ci_lower', 'ci_upper', 'welch_df']].to_string(index=False))
