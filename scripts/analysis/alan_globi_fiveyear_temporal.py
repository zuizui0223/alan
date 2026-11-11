#!/usr/bin/env python3
"""GloBI-positive plant photographs: five observed-year island/mainland ALAN dynamics.

Uses PRIVATE source-verified original photo × accepted taxon × same-year 4-month
VIIRS site-year aggregates. No raw site/species records are exported.
This is not an annual plant population series or documented local pollination.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsmr
from scipy.stats import t as student_t, f as fdist

YEARS = (2019, 2021, 2023, 2024, 2025)
NUISANCE = ('logpop', 'logphotos', 'logspecies', 'bio5', 'lograin', 'elevation_km')

def load(path):
    d = pd.read_csv(path, low_memory=False)
    required = {'target_id', 'land_class', 'region5', 'year', 'n_species',
        'photo_attempts', 'viirs_median', 'pop_density', 'bio5', 'rain',
        'elevation', 'all_moth', 'all_butterfly', 'non_inat_moth',
        'non_inat_butterfly'}
    if required - set(d):
        raise ValueError(f'Missing original-source columns: {required - set(d)}')
    if d.target_id.duplicated().any() or set(d.year) != set(YEARS):
        raise ValueError('Repeated photographic source or changed five-year window')
    if set(d.land_class) != {'island', 'mainland'}:
        raise ValueError('Unexpected land classification')
    for col in ('all_moth', 'all_butterfly', 'non_inat_moth', 'non_inat_butterfly'):
        if not (d[col].ge(0) & d[col].le(d.n_species)).all():
            raise ValueError('Positive GloBI taxon count exceeds photographed accepted species')
    d['island'] = d.land_class.eq('island').astype(float)
    d['region_year'] = d.year.astype(str) + '|' + d.region5.astype(str)
    d['site'] = d.target_id.str.split('|', regex=False).map(
        lambda x: '|'.join(x[1:]) if len(x) == 3 else '')
    if d.site.eq('').any() or d.duplicated(['year', 'site']).any():
        raise ValueError('Invalid original year|cell|GSHHG landmass identity')
    d['loglight'] = np.log1p(d.viirs_median.clip(lower=0))
    d['logpop'] = np.log1p(d.pop_density)
    d['logphotos'] = np.log1p(d.photo_attempts)
    d['logspecies'] = np.log1p(d.n_species)
    d['lograin'] = np.log1p(d.rain.clip(lower=0))
    d['elevation_km'] = d.elevation / 1000
    d['light_island'] = d.island * d.loglight
    d['share_all_moth'] = 100 * d.all_moth / d.n_species
    d['share_all_butterfly'] = 100 * d.all_butterfly / d.n_species
    d['share_non_inat_moth'] = 100 * d.non_inat_moth / d.n_species
    d['share_non_inat_butterfly'] = 100 * d.non_inat_butterfly / d.n_species
    d['share_moth_minus_butterfly'] = d.share_all_moth - d.share_all_butterfly
    d['share_non_inat_moth_minus_butterfly'] = (
        d.share_non_inat_moth - d.share_non_inat_butterfly)
    return d

def select(d, min_photos=3, min_species=2, light_common=False):
    z = d.loc[d.photo_attempts.ge(min_photos) & d.n_species.ge(min_species)].copy()
    z = z.replace([np.inf, -np.inf], np.nan).dropna(
        subset=['loglight', 'logpop', 'logphotos', 'logspecies', 'bio5',
                'lograin', 'elevation_km'])
    n = z.groupby('region_year').land_class.nunique()
    z = z.loc[z.region_year.isin(n.index[n.eq(2)])].copy()
    if light_common:
        bounds = z.groupby(['region_year', 'land_class']).loglight.agg(['min', 'max'])
        lower = bounds['min'].unstack().max(axis=1)
        upper = bounds['max'].unstack().min(axis=1)
        valid = lower.index[upper.gt(lower)]
        z = z.loc[z.region_year.isin(valid)].copy()
        z = z.loc[z.loglight.ge(z.region_year.map(lower)) &
                  z.loglight.le(z.region_year.map(upper))].copy()
        n = z.groupby('region_year').land_class.nunique()
        z = z.loc[z.region_year.isin(n.index[n.eq(2)])].copy()
    return z

def balanced_weights(d):
    n = d.groupby(['region_year', 'land_class']).size()
    index = pd.MultiIndex.from_frame(d[['region_year', 'land_class']])
    w = n.groupby(level=0).min().reindex(d.region_year).to_numpy(float) / (
        n.reindex(index).to_numpy(float))
    if not np.isfinite(w).all() or not (w > 0).all():
        raise ValueError('Cannot balance region-year site counts')
    return w

def fit(d, terms, response, *, fe=('region_year',), site_fe_df=False,
        min_clusters=12):
    if len(d) < 50:
        return {'estimable': False, 'reason': 'too_few_source_photographs'}
    a = d.copy()
    for key, val in terms.items():
        a[key] = a.eval(val) if isinstance(val, str) else val(a)
    columns = list(terms)
    X0 = a[columns].to_numpy(float)
    Y0 = a[response].to_numpy(float)
    if not np.isfinite(X0).all() or not np.isfinite(Y0).all():
        raise ValueError('Nonfinite photographic response or source predictor')
    w = balanced_weights(a)
    sw = np.sqrt(w)
    V = np.column_stack([X0, Y0]) * sw[:, None]
    geo, names = pd.factorize(a.region5, sort=True)
    G = len(names)
    if G < min_clusters:
        return {'estimable': False, 'reason': 'too_few_geographic_blocks', 'n_regions': G}
    fe_counts = {}
    if len(fe) == 1:
        group, ids = pd.factorize(a[fe[0]], sort=True)
        if (group < 0).any():
            raise ValueError('Missing geographic fixed effect')
        totals = np.bincount(group, weights=w, minlength=len(ids))
        fe_counts[fe[0]] = len(ids)
        for j in range(V.shape[1]):
            means = np.bincount(group, weights=w*(V[:, j]/sw), minlength=len(ids)) / totals
            V[:, j] -= sw * means[group]
    else:
        blocks = []
        for col in fe:
            idx, levels = pd.factorize(a[col], sort=True)
            if (idx < 0).any():
                raise ValueError('Missing site/year fixed effect')
            fe_counts[col] = len(levels)
            blocks.append(sparse.coo_matrix((np.ones(len(a)),
                (np.arange(len(a)), idx)), shape=(len(a), len(levels))).tocsr())
        Z = sparse.hstack(blocks, format='csr').multiply(sw[:, None]).tocsr()
        for j in range(V.shape[1]):
            coef, status, iters, residual, *_ = lsmr(
                Z, V[:, j], atol=1e-11, btol=1e-11, maxiter=5000)
            if status not in (0, 1, 2, 4, 5) or not np.isfinite(residual):
                raise ValueError('Site/year fixed-effect projection failed')
            V[:, j] -= Z @ coef
    X, Y = V[:, :-1], V[:, -1]
    coefficient, _, rank, _ = np.linalg.lstsq(X, Y, rcond=1e-12)
    if rank != len(columns):
        return {'estimable': False, 'reason': 'rank_deficiency', 'rank': int(rank)}
    residual = Y - X @ coefficient
    scores = X * residual[:, None]
    bygeo = np.column_stack([
        np.bincount(geo, weights=scores[:, j], minlength=G)
        for j in range(len(columns))])
    nfe = sum(fe_counts.values()) - max(0, len(fe_counts)-1) if site_fe_df else 0
    df = len(a) - len(columns) - nfe
    if df <= 0:
        return {'estimable': False, 'reason': 'no_residual_df'}
    bread = np.linalg.inv(X.T @ X)
    covariance = bread @ bygeo.T @ bygeo @ bread * (G/(G-1))*((len(a)-1)/df)
    covariance = (covariance + covariance.T)/2
    tcrit = float(student_t.ppf(.975, G-1))
    effects = {}
    for i, col in enumerate(columns):
        se = math.sqrt(float(covariance[i, i])) if covariance[i, i] > 0 else None
        effects[col] = {
            'beta': float(coefficient[i]),
            'se': se,
            'ci95': [float(coefficient[i]-tcrit*se),
                     float(coefficient[i]+tcrit*se)] if se is not None else None,
            'p': float(2*student_t.sf(abs(coefficient[i]/se), G-1))
                  if se is not None else None}
    return {
        'estimable': True, 'n_site_years': len(a),
        'n_island': int(a.island.sum()),
        'n_mainland': int((1-a.island).sum()),
        'n_sites': int(a.site.nunique()),
        'n_regions': G, 'n_year_region_groups': int(a.region_year.nunique()),
        'fixed_effects': fe_counts, 'coefficients': effects,
        'terms': columns, '_coef': coefficient, '_cov': covariance}

def basic_terms():
    return {'island': 'island', 'loglight': 'loglight',
            'light_island': 'light_island', **{x:x for x in NUISANCE}}

def dynamic_terms():
    terms = {}
    for year in YEARS:
        terms[f'island_{year}'] = lambda d, year=year: d.island * d.year.eq(year)
    for year in YEARS:
        terms[f'light_{year}'] = lambda d, year=year: d.loglight * d.year.eq(year)
    for year in YEARS:
        terms[f'island_light_{year}'] = lambda d, year=year: d.light_island * d.year.eq(year)
    return {**terms, **{x:x for x in NUISANCE}}

def trend_terms():
    return {'island': 'island', 'loglight': 'loglight',
        'light_island': 'light_island',
        'year_island': lambda d: (d.year-2019)*d.island,
        'year_light': lambda d: (d.year-2019)*d.loglight,
        'year_island_light': lambda d: (d.year-2019)*d.light_island,
        **{x:x for x in NUISANCE}}

def year_heterogeneity(result):
    if not result.get('estimable'):
        return {'estimable': False}
    cols = result['terms']
    R = np.zeros((4, len(cols)))
    for i, yr in enumerate(YEARS[1:]):
        R[i, cols.index(f'island_light_{yr}')] = 1
        R[i, cols.index(f'island_light_{YEARS[0]}')] = -1
    delta = R @ result['_coef']
    cov = R @ result['_cov'] @ R.T
    if np.linalg.matrix_rank(cov) < 4:
        return {'estimable': False, 'reason': 'singular_year_difference_covariance'}
    F = float(delta @ np.linalg.solve(cov, delta) / 4)
    return {'estimable': True, 'F': F, 'df_num': 4,
            'df_den': result['n_regions']-1,
            'p_exploratory': float(fdist.sf(F, 4, result['n_regions']-1))}

def safe(result):
    return {k:v for k,v in result.items() if not k.startswith('_')}

def repeated(d):
    n = d.groupby('site').year.nunique()
    z = d.loc[d.site.isin(n.index[n.ge(2)])].copy()
    support = z.groupby('region_year').land_class.nunique()
    z = z.loc[z.region_year.isin(support.index[support.eq(2)])].copy()
    n = z.groupby('site').year.nunique()
    return z.loc[z.site.isin(n.index[n.ge(2)])].copy()

def run(source, photo_season=None):
    d = load(source)
    if photo_season is not None:
        months = pd.read_csv(photo_season)
        fields = {'target_id', 'valid_dates', 'sin_month','cos_month','south_sin','south_cos'}
        if fields - set(months) or months.target_id.duplicated().any():
            raise ValueError('Photo month source missing fields or duplicated')
        if set(months.target_id) != set(d.target_id):
            raise ValueError('Photo month source not aligned to original photo cell-year population')
        d = d.merge(months[list(fields)], on='target_id', how='inner',
                    validate='one_to_one')
        if d[list(fields-{'target_id'})].isna().any().any():
            raise ValueError('No missing photograph season imputation allowed')
    result = {
        'schema': 'alan_globi_all_source_photo_temporal_v1',
        'years': list(YEARS), 'missing_calendar_years': [2020, 2022],
        'light_source': 'SAME-YEAR 4-month VIIRS, NOT annual 12-month light',
        'original_PRIVATE_source_sha256': hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        'original_PRIVATE_season_sha256': hashlib.sha256(Path(photo_season).read_bytes()).hexdigest()
               if photo_season else None,
        'GloBI_plant_record_somewhere_NOT_documented_local_pollination': True,
        'no_population_loss_or_causal_ALAN_claim': True,
        'cohorts': {}}
    outcomes = {
        'all_moth': 'share_all_moth',
        'non_inat_moth': 'share_non_inat_moth',
        'butterfly': 'share_all_butterfly',
        'moth_minus_butterfly': 'share_moth_minus_butterfly',
        'non_inat_moth_minus_butterfly': 'share_non_inat_moth_minus_butterfly'}
    for name, minphoto, minspecies in (('broad', 3, 2), ('strict', 5, 3)):
        matched = select(d, minphoto, minspecies)
        source_n = matched.groupby(['year','land_class']).size().unstack(fill_value=0)
        c = {'n':len(matched),'geographic_clusters':matched.region5.nunique(),
             'observed_photo_site_years': {str(int(yr)): {'island':int(row.island),
                 'mainland':int(row.mainland)} for yr,row in source_n.iterrows()}}
        for label, outcome in outcomes.items():
            yearly = {}
            for yr in YEARS:
                fit_year = fit(matched.loc[matched.year.eq(yr)],basic_terms(),outcome)
                yearly[str(yr)] = safe(fit_year)
            joint = fit(matched,dynamic_terms(),outcome)
            trend = fit(matched,trend_terms(),outcome)
            pooled = fit(matched,basic_terms(),outcome)
            omitted = {str(yr):safe(fit(matched.loc[matched.year.ne(yr)],
                basic_terms(),outcome)) for yr in YEARS}
            c[label]={'per_year':yearly,'joint_year_slope_model':safe(joint),
                      'year_heterogeneity':year_heterogeneity(joint),
                      'pooled':safe(pooled),
                      'pooled_leave_year_out':omitted,
                      'linear_calendar_trend':safe(trend)}
            if photo_season:
                term = trend_terms()
                term.update({x:x for x in ('sin_month','cos_month',
                                           'south_sin','south_cos')})
                c[label]['season_adjusted_calendar_trend']=safe(
                    fit(matched,term,outcome))
        sr = repeated(matched)
        c['repeated_site_support']={'site_years':len(sr),'sites':sr.site.nunique(),
             'island_site_years':int(sr.island.sum()),
             'mainland_site_years':int((1-sr.island).sum()),
             'geographic_blocks':sr.region5.nunique()}
        site_terms={'loglight':'loglight','light_island':'light_island',
            'logpop':'logpop','logphotos':'logphotos','logspecies':'logspecies'}
        site_time = {**site_terms,
            'year_island':lambda d:(d.year-2019)*d.island,
            'year_light':lambda d:(d.year-2019)*d.loglight,
            'year_island_light':lambda d:(d.year-2019)*d.light_island}
        for label in ('all_moth','non_inat_moth','moth_minus_butterfly',
                      'non_inat_moth_minus_butterfly'):
            c[label]['site_FE']=safe(fit(sr,site_terms,outcomes[label],
                                fe=('site','region_year'),site_fe_df=True))
            c[label]['site_FE_time_trend']=safe(fit(sr,site_time,outcomes[label],
                                fe=('site','region_year'),site_fe_df=True))
        eligible_by_year=[yr for yr,v in c['all_moth']['per_year'].items()
                     if not v.get('estimable') or v.get('n_island',0)<30]
        c['year_specific_support_warning']={
            'years_with_too_few_photo_island_sites_or_clusters':eligible_by_year,
            'joint_year_heterogeneity_not_reliable_if_nonempty':bool(eligible_by_year)}
        result['cohorts'][name]=c
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--private-source',type=Path,required=True)
    p.add_argument('--private-season',type=Path,default=None)
    p.add_argument('--aggregate-output',type=Path,required=True)
    args=p.parse_args()
    if args.aggregate_output.exists():
        raise ValueError('Refuse to overwrite temporal research result')
    results=run(args.private_source,args.private_season)
    args.aggregate_output.parent.mkdir(parents=True,exist_ok=True)
    args.aggregate_output.write_text(json.dumps(results,ensure_ascii=False,
        indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'years':results['years'],
         'broad_source_n':results['cohorts']['broad']['n']},ensure_ascii=False))
if __name__=='__main__':main()
