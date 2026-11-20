"""Exploratory two-programme meta-analysis; not a final systematic-review result."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
CACHE = REPO.parent / 'gbif_annual/alan_experiment_sources'
OUT = ROOT / 'capture_origin_synthesis'
OUT.mkdir(exist_ok=True)

raw = CACHE / 'parkinson2020/Aquatic and Terrestrial Invertebrate Inputs.csv'
blob = raw.read_bytes()
assert hashlib.sha1(b'blob ' + str(len(blob)).encode() + b'\0' + blob).hexdigest() == '4e2b7065d67e0604ca121ea26c96b5db6d64de94'
p = pd.read_csv(raw)
u20 = pd.DataFrame(dict(
    programme_id='Parkinson2016', unit_id=p['Mesocosm ID'].astype(str),
    habitat='freshwater_littoral', colour=np.where(p.Treatment.eq('ALAN'), 'original_treatment', 'dark'),
    aquatic=p['# individual aquatic invertberates'],
    terrestrial=p['# of individual terrestrial invertebrates'],
    ephydridae=p['# Ephydridae (I)'], conflicting_origin='aquatic'))
u20['family_aquatic'] = p[[c for c in p if re.search(r'\([A-Z]\)', c)]].sum(axis=1)
u20['family_terrestrial'] = p[[c for c in p if re.search(r'\([a-z]\)', c)]].sum(axis=1)
u23 = pd.read_csv(ROOT / 'parkinson2023_sources/trap_units.csv')
taxa = pd.read_csv(ROOT / 'parkinson2023_sources/taxon_counts.csv')
eph = taxa[taxa.family_source_spelling.eq('ephydridae')].groupby('unit_id')['count'].sum()
u23 = u23.rename(columns={'aquatic_origin_capture':'aquatic', 'terrestrial_origin_capture':'terrestrial'})
# An absent taxon row means zero in this exhaustive trap-count table, not an
# absence inferred from opportunistic occurrence records.
u23['ephydridae'] = u23.unit_id.map(eph).fillna(0)
u23['conflicting_origin'] = 'terrestrial'
u23['family_aquatic'] = u23.aquatic
u23['family_terrestrial'] = u23.terrestrial
units = pd.concat([u20, u23[u20.columns]], ignore_index=True)
assert len(units) == 50
units.to_csv(OUT / 'source_units.csv', index=False)
units[(units.aquatic != units.family_aquatic) | (units.terrestrial != units.family_terrestrial)].to_csv(OUT / 'summary_family_discrepancies.csv', index=False)
units.groupby(['programme_id','colour'], sort=False).ephydridae.sum().to_csv(OUT / 'conflicting_family_counts.csv')
original = pd.read_csv(ROOT / 'capture_origin_contrasts.csv')
original_v = pd.read_csv(ROOT / 'capture_origin_contrast_covariance.csv', index_col=0)
programme_rows, pooled_rows, weights = [], [], []
audits = {}

# First verify the old nine contrasts from unit-level values, independently of
# the previous covariance assembly. This does not endorse their taxonomy.
estimates, scores = [], []
for r in original.itertuples():
    score = np.zeros(len(units)); estimate = 0.
    for colour, sign in [(r.colour, 1), ('dark', -1)]:
        g = units[(units.programme_id == r.programme_id) & (units.habitat == r.habitat) & (units.colour == colour)]
        a, b = g.aquatic.to_numpy(), g.terrestrial.to_numpy()
        estimate += sign * np.log(a.mean()/b.mean())
        score[g.index] = sign * ((a-a.mean())/a.mean() - (b-b.mean())/b.mean()) / np.sqrt(len(g)*(len(g)-1))
    estimates.append(estimate); scores.append(score)
assert np.allclose(estimates, original.lnRR_aquatic_minus_terrestrial, atol=1e-12)
assert np.allclose(np.array(scores) @ np.array(scores).T,
                   original_v.loc[original.contrast_id, original.contrast_id], atol=1e-12)
pd.DataFrame(np.array(scores) @ np.array(scores).T,
             index=original.contrast_id, columns=original.contrast_id).to_csv(
                 OUT / 'author_defined_contrast_covariance.csv')

for variant in ['source_defined_menu_mean_counts', 'family_counts_excluding_Ephydridae']:
    u = units.copy()
    if variant != 'source_defined_menu_mean_counts':
        for endpoint in ['aquatic', 'terrestrial']:
            u[endpoint] = u['family_'+endpoint]
            mask = u.conflicting_origin.eq(endpoint)
            u.loc[mask, endpoint] -= u.loc[mask, 'ephydridae']
    assert (u[['aquatic','terrestrial']] >= 0).all().all()
    estimates, scores, arm_means = [], [], []
    programmes = list(original.programme_id.unique())
    for programme in programmes:
        score = np.zeros(len(u)); estimate = 0.; means = {}
        up = u[u.programme_id.eq(programme)]
        for is_lit, sign in [(True, 1), (False, -1)]:
            arm = up[up.colour.ne('dark') == is_lit]
            cells = list(arm.groupby(['habitat', 'colour'], sort=False))
            count_means = np.mean([g[['aquatic','terrestrial']].mean().to_numpy() for _,g in cells], axis=0)
            ma, mb = count_means
            assert ma > 0 and mb > 0
            estimate += sign*np.log(ma/mb)
            label = 'lit' if is_lit else 'dark'
            means.update({label+'_aquatic':ma, label+'_terrestrial':mb})
            for (habitat, colour), g in cells:
                n = len(g)
                aq, te = g.aquatic.to_numpy(), g.terrestrial.to_numpy()
                # Fixed, equal treatment-cell weights; do not count variation
                # between deliberately different colours as sampling replication.
                score[g.index] = sign/len(cells)*((aq-aq.mean())/ma-(te-te.mean())/mb)/np.sqrt(n*(n-1))
                weights.append(dict(variant=variant, programme_id=programme, habitat=habitat,
                    colour=colour, arm=label, cell_weight=1/len(cells), units=n))
        estimates.append(estimate); scores.append(score); arm_means.append(means)
    py = np.array(estimates)
    pv = np.array(scores) @ np.array(scores).T
    assert len(py) == 2 and abs(pv[0,1]) < 1e-12
    # k=2 REML depends only on the independent difference: Var(y1-y2)=v1+v2+2*tau2.
    vs = np.diag(pv)
    tau2 = max(0., ((py[0]-py[1])**2 - vs.sum())/2)
    def nll(x):
        vv = vs + x
        ww = 1/vv
        mu = np.sum(ww*py)/ww.sum()
        return .5*(np.log(vv).sum() + np.log(ww.sum()) + np.sum(ww*(py-mu)**2))
    fit = minimize_scalar(nll, bounds=(0, max(10., tau2*4)), method='bounded', options={'xatol':1e-12})
    numeric_tau2 = 0. if nll(0.) <= fit.fun else float(fit.x)
    assert np.isclose(tau2, numeric_tau2, atol=1e-7)
    w = 1/(vs+tau2)
    mu = float(np.sum(w*py)/w.sum())
    q = float(np.sum(w*(py-mu)**2))  # df=k-1=1
    variance = float(max(1., q)/w.sum())
    radius = float(t.ppf(.975, 1)*np.sqrt(variance))
    pooled_rows.append(dict(variant=variant, programmes=2, source_contrasts=9, log_ratio=mu,
        ratio=float(np.exp(mu)), lower=float(np.exp(mu-radius)), upper=float(np.exp(mu+radius)),
        tau2_REML=tau2, q=q, mKH_variance=variance, df=1,
        status='exploratory_not_final_adjudication', estimand='origin_response_of_fixed_menu_arithmetic_mean_captures'))
    for j, programme in enumerate(programmes):
        programme_rows.append(dict(variant=variant, programme_id=programme,
            log_ratio=py[j], sampling_variance=pv[j,j], ratio=np.exp(py[j]),
            meta_weight_fraction=w[j]/w.sum(), **arm_means[j]))
    audits[variant] = dict(covariance_min_eigenvalue=float(np.linalg.eigvalsh(pv).min()),
        analytic_tau2=tau2, numerical_tau2=numeric_tau2, source_units=len(u),
        cross_programme_covariance=float(pv[0,1]))

for filename, rows in [('programme_estimates.csv',programme_rows), ('pooled_estimates.csv',pooled_rows),
                       ('aggregation_weights.csv',weights)]:
    pd.DataFrame(rows).to_csv(OUT / filename, index=False)
audit = dict(variants=audits, source_defined_estimates_and_covariance_reproduced=True,
    input_sha256={str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path):
                  hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in [raw, ROOT/'parkinson2023_sources/trap_units.csv',
                               ROOT/'parkinson2023_sources/taxon_counts.csv',
                               ROOT/'capture_origin_contrasts.csv',
                               ROOT/'capture_origin_contrast_covariance.csv']},
    covariance_method='independent-arm centered delta scores, shared controls retained',
    interval_method='modified Knapp-Hartung max(1,q), t1',
    master_catalogue_changed=False, final_eligibility=False,
    limitations=['Two programmes in Michigan; no broad geographic inference',
        'Treatment menus differ; not a common spectrum or irradiance effect',
        'Ephydridae exclusion resolves one documented conflict, not a universal taxonomic validation',
        'Arithmetic capture-menu mean is not the geometric mean of contrast ratios',
        'Fish allocation identity unavailable in Parkinson2016',
        'Three traps per treatment cell and possible attraction interference in SevenPonds2022',
        'Capture is not emergence, predation, pollination or plant decline'])
(OUT / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
print(pd.DataFrame(pooled_rows).to_string(index=False))
