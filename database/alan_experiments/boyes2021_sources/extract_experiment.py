"""Date-adjusted sweep-net counts: one experiment, two dependent contrasts."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import t
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    args = parser.parse_args()
    out = Path(__file__).resolve().parent
    raw = pd.read_csv(args.source)
    assert len(raw) == 27 and not raw.duplicated(['Site_Name', 'TransectID']).any()
    assert raw.Date.nunique() == 9 and raw.Site_Name.nunique() == 2
    assert raw.groupby('Treatment').size().to_dict() == {'HPS': 9, 'LED': 9, 'Unlit': 9}
    assert raw.Larvae.notna().all() and (raw.Larvae >= 0).all()
    raw.to_csv(out / 'source_transect_counts.csv', index=False)
    raw.groupby(['Date', 'Treatment']).size().rename('transects').reset_index().to_csv(
        out / 'date_treatment_allocation.csv', index=False)
    raw.groupby('Treatment').Larvae.agg(['size', 'sum', 'mean', 'std']).to_csv(
        out / 'arm_summaries.csv')
    data = raw.rename(columns={'Treatment': 'treatment'})
    formula = 'Larvae ~ C(treatment, Treatment(reference="Unlit")) + C(Date)'
    model = smf.glm(formula, data, family=sm.families.Poisson()).fit()
    clustered = smf.glm(formula, data, family=sm.families.Poisson()).fit(
        cov_type='cluster', cov_kwds={'groups': data.Date, 'use_correction': True})
    names = [x for x in model.params.index if x.startswith('C(treatment')]
    rows, covariance = [], []
    for label, fit, critical in [('poisson_model', model, 1.959963984540054),
                                 ('date_cluster_t8', clustered, t.ppf(.975, 8))]:
        for name in names:
            beta, se = float(fit.params[name]), float(fit.bse[name])
            rows.append(dict(model=label, treatment=name.rsplit('[T.', 1)[1][:-1],
                reference='Unlit', log_count_ratio=beta, SE=se,
                count_ratio=np.exp(beta), ci_lower=np.exp(beta-critical*se),
                ci_upper=np.exp(beta+critical*se), critical_value=critical,
                study_id='Boyes2021', endpoint='short_term_sweep_net_count',
                pooling_status='not_adopted_into_effect_catalogue'))
        for a in names:
            for b in names:
                covariance.append(dict(model=label, term1=a, term2=b,
                                       covariance=float(fit.cov_params().loc[a,b])))
    pd.DataFrame(rows).to_csv(out / 'model_contrasts.csv', index=False)
    pd.DataFrame(covariance).to_csv(out / 'model_covariance.csv', index=False)
    audit = dict(data_doi='10.5285/7e7ce7c1-aa16-484f-86b1-90c5b088a64c',
        article_doi='10.1126/sciadv.abi8322', source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        transects=27, sites=2, sampling_nights=9, formula=formula,
        raw_years=sorted(pd.to_datetime(raw.Date, dayfirst=True).dt.year.unique().tolist()),
        article_methods_year=2020, year_discrepancy_unresolved=True,
        pearson_dispersion=float(model.pearson_chi2/model.df_resid),
        fitting='date fixed effects; not an exact replay of original random-intercept GLMM',
        cluster_interval='small-sample-corrected date sandwich SE with t8 interval; nine clusters only',
        statsmodels=statsmodels.__version__, scipy=scipy.__version__,
        independent_experiments=1, meta_effects_added=0,
        limitation='Counts measure short-term sweep-net availability/activity, not population mortality, consumption, or pollination; same rig per light type may limit fixture generalization.')
    (out/'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
