"""Named-species ALAN contrasts with covariance from shared physical mesocosms.

This is one post hoc, complete-species subset of Penryn2012_2016, not a
wetness meta-regression. Run from any directory with --out DIRECTORY.
"""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import t

BASE = Path(__file__).resolve().parent
FOCAL = 'Lotus pedunculatus'
COMPARATORS = ['Agrostis tenuis', 'Anthoxanthum odoratum', 'Hypochaeris radicata']
EXPECTED_SPECIES = set(COMPARATORS + [FOCAL])
LIGHTS = ['L2', 'L3']
TROPHIC = ['T1', 'T2', 'T3']


def select_units(data, survival, community):
    """Require complete positive rows, consistent allocation and S1 agreement."""
    if data.duplicated(['species', 'mesocosm']).any():
        raise ValueError('Duplicate species-by-mesocosm observations')
    if community.mesocosm.duplicated().any() or len(community) != 54:
        raise ValueError('Expected 54 unique physical mesocosms')
    if survival.duplicated(['species', 'light', 'trophic']).any():
        raise ValueError('Duplicate survival-table cells')
    cells = community.groupby(['light', 'trophic']).size()
    if len(cells) != 9 or not cells.eq(6).all():
        raise ValueError('Expected nine treatment cells with six mesocosms each')
    joined = data.merge(community, on='mesocosm', validate='many_to_one',
                        suffixes=('', '_allocation'), how='left', indicator=True)
    if not (joined['_merge'].eq('both').all()
            and joined.light.eq(joined.light_allocation).all()
            and joined.trophic.eq(joined.trophic_allocation).all()):
        raise ValueError('Biomass and community treatment allocations disagree')
    selection = []
    for species, group in data.groupby('species', sort=True):
        s = survival[survival.species.eq(species)]
        complete = (len(group) == 54 and set(group.mesocosm) == set(community.mesocosm)
                    and np.isfinite(group.total_g).all() and group.total_g.gt(0).all())
        agreement = (len(s) == 9 and s.status.eq('agrees').all()
                     and s.surviving_mesocosms_S1.eq(6).all())
        components = not group.component_inconsistency.any()
        selection.append(dict(species=species, recorded_units=len(group),
            complete_positive_54=complete, survival_all_nine_cells_agree=agreement,
            component_values_consistent=components,
            included=bool(complete and agreement and components)))
    selection = pd.DataFrame(selection)
    eligible = set(selection.loc[selection.included, 'species'])
    if eligible != EXPECTED_SPECIES:
        raise ValueError(f'Eligible source subset changed: {sorted(eligible)}; review design before rerunning')
    return data[data.species.isin(eligible)].copy(), selection


def species_effects(units):
    """Return arithmetic-mean lnRRs and their mesocosm-aligned Gram covariance.

    Each cell contributes signed (x / mean(x) - 1) / sqrt(n*(n-1)).
    This exactly reproduces the first-order delta-method covariance and
    preserves species covariance and shared dark controls.
    """
    if units.duplicated(['species', 'mesocosm']).any():
        raise ValueError('Duplicate species-by-unit rows')
    if not np.isfinite(units.total_g).all() or not units.total_g.gt(0).all():
        raise ValueError('Biomass must be finite and positive; no pseudocounts')
    if not units.groupby('mesocosm')[['light', 'trophic']].nunique().eq(1).all().all():
        raise ValueError('Treatment assignment changes within a physical mesocosm')
    plot_ids = sorted(units.mesocosm.unique())
    position = {p: i for i, p in enumerate(plot_ids)}
    rows, columns = [], []
    for (species, trophic), group in units.groupby(['species', 'trophic'], sort=True):
        dark = group[group.light.eq('L1')].set_index('mesocosm').total_g
        for light in LIGHTS:
            lit = group[group.light.eq(light)].set_index('mesocosm').total_g
            if len(dark) < 2 or len(lit) < 2:
                raise ValueError('At least two independent units per arm required')
            score = np.zeros(len(plot_ids))
            for values, sign in [(lit, 1), (dark, -1)]:
                n = len(values)
                centered = sign*(values/values.mean()-1)/np.sqrt(n*(n-1))
                for plot, value in centered.items():
                    score[position[plot]] = value
            rows.append(dict(effect_id=f'{species}|{trophic}|{light}', species=species,
                trophic=trophic, light=light, n_lit=len(lit), n_dark=len(dark),
                mean_lit_g=lit.mean(), mean_dark_g=dark.mean(),
                lnRR=np.log(lit.mean()/dark.mean())))
            columns.append(score)
    effects = pd.DataFrame(rows)
    scores = np.column_stack(columns)
    covariance = scores.T@scores
    effects['sampling_variance'] = np.diag(covariance)
    effects['ratio'] = np.exp(effects.lnRR)
    return effects, covariance, scores


def contrast_matrix(effects, comparators=COMPARATORS):
    """Six cell, two marginal, four treatment-modification contrasts."""
    if not comparators or FOCAL in comparators:
        raise ValueError('Nonempty comparator set must exclude the focal species')
    lookup = {(r.species, r.trophic, r.light): i for i, r in enumerate(effects.itertuples())}
    cell, rows, vectors = {}, [], []
    for trophic in TROPHIC:
        for light in LIGHTS:
            vector = np.zeros(len(effects))
            vector[lookup[FOCAL, trophic, light]] = 1
            for species in comparators:
                vector[lookup[species, trophic, light]] = -1/len(comparators)
            cell[trophic, light] = vector
            rows.append(dict(contrast_id=f'relative_{trophic}_{light}', kind='cell',
                             trophic=trophic, light=light, family_size=6))
            vectors.append(vector)
    for light in LIGHTS:
        rows.append(dict(contrast_id=f'relative_average_{light}', kind='marginal',
                         trophic='equal_weight_T1_T2_T3', light=light, family_size=2))
        vectors.append(np.mean([cell[t, light] for t in TROPHIC], axis=0))
    for high, low, label in [('T2', 'T1', 'herbivore_package'), ('T3', 'T2', 'predator_addition')]:
        for light in LIGHTS:
            rows.append(dict(contrast_id=f'{label}_{light}', kind='modifier',
                             trophic=f'{high}_minus_{low}', light=light, family_size=4))
            vectors.append(cell[high, light]-cell[low, light])
    return pd.DataFrame(rows), np.vstack(vectors)


def estimate_contrasts(meta, matrix, theta, covariance):
    out = meta.copy()
    out['estimate'] = matrix@np.asarray(theta)
    v = matrix@covariance@matrix.T
    if np.linalg.eigvalsh(v).min() < -1e-10:
        raise ValueError('Contrast covariance is not positive semidefinite')
    out['sampling_variance'] = np.maximum(np.diag(v), 0)
    out['se'] = np.sqrt(out.sampling_variance)
    out['ratio'] = np.exp(out.estimate)
    out['df_conservative'] = 5
    half = t.ppf(.975, 5)*out.se
    out['lower_95_t5'] = np.exp(out.estimate-half)
    out['upper_95_t5'] = np.exp(out.estimate+half)
    family_half = t.ppf(1-.05/(2*out.family_size), 5)*out.se
    out['lower_95_family_t5'] = np.exp(out.estimate-family_half)
    out['upper_95_family_t5'] = np.exp(out.estimate+family_half)
    return out, v


def mean_log_sensitivity(units):
    """Change estimand: mean mesocosm log-ratio instead of log of arm means."""
    wide = units.pivot(index=['mesocosm', 'light', 'trophic'], columns='species', values='total_g')
    z = np.log(wide[FOCAL])-np.log(wide[COMPARATORS]).mean(axis=1)
    d = z.rename('index').reset_index()
    rows = []
    for trophic in TROPHIC:
        g = d[d.trophic.eq(trophic)]
        for light in LIGHTS:
            estimate = g.loc[g.light.eq(light), 'index'].mean()-g.loc[g.light.eq('L1'), 'index'].mean()
            rows.append(dict(trophic=trophic, light=light, estimate=estimate, ratio=np.exp(estimate)))
    cell = pd.DataFrame(rows)
    for light in LIGHTS:
        estimate = cell.loc[cell.light.eq(light), 'estimate'].mean()
        rows.append(dict(trophic='equal_weight_T1_T2_T3', light=light,
                         estimate=estimate, ratio=np.exp(estimate)))
    return pd.DataFrame(rows), d


def save_cov(path, ids, matrix):
    pd.DataFrame([dict(effect_id_i=a, effect_id_j=b, covariance=matrix[i, j])
                  for i, a in enumerate(ids) for j, b in enumerate(ids)]).to_csv(path, index=False)


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    inputs = ['species_biomass_source_rows.csv', 'species_survival_reconciliation.csv',
              'anic2022_community_units.csv']
    hashes = {name: hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in inputs}
    data = pd.read_csv(BASE/inputs[0], dtype={name: str for name in
                       ['total_source', 'root_source', 'shoot_source']})
    units, selection = select_units(data, pd.read_csv(BASE/inputs[1]), pd.read_csv(BASE/inputs[2]))
    effects, covariance, scores = species_effects(units)
    meta, matrix = contrast_matrix(effects)
    estimates, contrast_cov = estimate_contrasts(meta, matrix, effects.lnRR, covariance)
    units.to_csv(out/'selected_source_rows.csv', index=False)
    selection.to_csv(out/'species_selection.csv', index=False)
    effects.to_csv(out/'species_effects.csv', index=False)
    save_cov(out/'species_covariance.csv', effects.effect_id.tolist(), covariance)
    pd.DataFrame(matrix, index=meta.contrast_id, columns=effects.effect_id).to_csv(out/'contrast_weights.csv')
    estimates.to_csv(out/'composition_and_trophic_contrasts.csv', index=False)
    save_cov(out/'contrast_covariance.csv', meta.contrast_id.tolist(), contrast_cov)
    mean_log, unit_index = mean_log_sensitivity(units)
    mean_log.to_csv(out/'mean_log_estimand_sensitivity.csv', index=False)
    unit_index.to_csv(out/'mesocosm_log_indices.csv', index=False)
    sensitivity = []
    for species in COMPARATORS:
        for mode, comp in [('single_comparator', [species]),
                           ('leave_one_comparator_out', [s for s in COMPARATORS if s != species])]:
            m, a = contrast_matrix(effects, comp)
            e, _ = estimate_contrasts(m, a, effects.lnRR, covariance)
            e['mode'] = mode; e['named_comparator'] = species
            sensitivity.append(e)
    sensitivity = pd.concat(sensitivity, ignore_index=True)
    sensitivity.to_csv(out/'comparator_sensitivity.csv', index=False)
    marginal = []
    for species in sorted(EXPECTED_SPECIES):
        for light in LIGHTS:
            vector = (effects.species.eq(species)&effects.light.eq(light)).to_numpy(dtype=float)/3
            estimate = vector@effects.lnRR.to_numpy()
            se = np.sqrt(vector@covariance@vector)
            marginal.append(dict(species=species, light=light, estimate=estimate,
                ratio=np.exp(estimate), se=se, lower_95_t5=np.exp(estimate-t.ppf(.975,5)*se),
                upper_95_t5=np.exp(estimate+t.ppf(.975,5)*se)))
    pd.DataFrame(marginal).to_csv(out/'species_marginal_effects.csv', index=False)
    audit = dict(schema='Penryn_complete_species_trophic_v1',
        base_commit='db82976f834cc048e76e69e6f621ba11e3ef33e5', input_sha256=hashes,
        original_pdf_sha256='baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0',
        source_doi='10.6084/m9.figshare.20304039', paper_doi='10.1002/ece3.9157',
        source_rows=len(units), physical_mesocosms=units.mesocosm.nunique(),
        independent_programmes=1, selected_species=sorted(EXPECTED_SPECIES),
        species_effects=len(effects), covariance_elements=int(covariance.size),
        derived_contrasts=len(estimates), covariance_min_eigenvalue=float(np.linalg.eigvalsh(covariance).min()),
        covariance_rank=int(np.linalg.matrix_rank(covariance)),
        imputed_values=0, pooled_estimates_added=0, main_catalogue_changed=False,
        marginal_results=estimates[estimates.kind.eq('marginal')].to_dict('records'),
        trophic_modifiers=estimates[estimates.kind.eq('modifier')].to_dict('records'),
        sensitivity_all_marginal_point_ratios_below_one=bool(
            sensitivity.loc[sensitivity.kind.eq('marginal'), 'ratio'].lt(1).all()),
        limitations=[
            'Post hoc complete-survivor species subset; focal species chosen after prior evidence',
            'One species above moisture F7; wetness, species identity and legume status inseparable',
            'Daytime light-preference scores are not sensitivity to artificial nighttime light',
            'T1 receives insecticide/molluscicide; T2 minus T1 is a treatment package',
            'First-order delta covariance and approximate t5 intervals with six units per cell',
            'No unmeasured competition mechanism, pollination, mediation, or macro decline effect',
            'Same Penryn programme as existing aphid, flowerhead and community-biomass effects',
            'AI calculation checks do not constitute human duplicate extraction'])
    (out/'analysis_summary.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(estimates[['contrast_id','ratio','lower_95_t5','upper_95_t5']].to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    run(parser.parse_args().out)
