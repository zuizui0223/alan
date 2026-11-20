"""Extract aphid-excluded sentinel biomass; audit distinct insect cohort."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source_dir', type=Path)
    source = parser.parse_args().source_dir
    out = Path(__file__).resolve().parent
    biomass = pd.read_csv(source/'Spectra_biomass.csv')
    insects = pd.read_csv(source/'Spectra_fieldexperiment.csv')
    assert len(biomass) == 56 and biomass.Cage.is_unique
    assert not biomass.duplicated(['Block', 'Treat']).any()
    assert biomass.groupby('Treat').size().eq(7).all()
    assert biomass[['ab.plant', 'root']].notna().all().all()
    assert (biomass[['ab.plant', 'root']] > 0).all().all()
    assert not insects.duplicated(['Cage', 'Week']).any()
    allocation = insects[['Cage', 'Block', 'Treatment']].drop_duplicates()
    assert allocation.Cage.is_unique
    joined = biomass.merge(allocation, on='Cage', how='outer',
                           suffixes=('_biomass', '_insects'), indicator=True)
    both = joined._merge.eq('both')
    assert joined.loc[both, 'Block_biomass'].eq(joined.loc[both, 'Block_insects']).all()
    assert joined.loc[both, 'Treat'].eq(joined.loc[both, 'Treatment']).all()
    joined.to_csv(out/'cross_endpoint_allocation.csv', index=False)
    biomass.to_csv(out/'source_biomass.csv', index=False)
    treatments = sorted(set(biomass.Treat)-{'C'})
    rows, influences = [], []
    for endpoint in ['ab.plant', 'root']:
        wide = biomass.pivot(index='Block', columns='Treat', values=endpoint).sort_index()
        assert wide.shape == (7, 8) and wide.notna().all().all()
        dark = wide['C'].to_numpy()
        for treatment in treatments:
            lit = wide[treatment].to_numpy()
            effect = np.log(lit.mean()/dark.mean())
            influence = (lit-lit.mean())/lit.mean() - (dark-dark.mean())/dark.mean()
            variance = np.var(influence, ddof=1)/7
            radius = t.ppf(.975, 6)*np.sqrt(variance)
            effect_id = f'Sanders2017_sentinel_{endpoint}_{treatment}'
            rows.append(dict(effect_id=effect_id, endpoint=endpoint, treatment=treatment,
                control='C', n_blocks=7, mean_lit=lit.mean(), mean_dark=dark.mean(),
                lnRR=effect, sampling_variance=variance, RR=np.exp(effect),
                ci_lower=np.exp(effect-radius), ci_upper=np.exp(effect+radius),
                status='original_dataset_extracted_final_article_pending',
                interpretation='aphid_excluded_sentinel_biomass_week13_not_trophic_mediation'))
            influences.append(influence)
    effects = pd.DataFrame(rows)
    cov = np.cov(np.stack(influences), ddof=1)/7
    assert np.allclose(np.diag(cov), effects.sampling_variance)
    assert np.linalg.eigvalsh(cov).min() > -1e-12
    effects.to_csv(out/'biomass_effects.csv', index=False)
    pd.DataFrame(cov, index=effects.effect_id, columns=effects.effect_id).to_csv(out/'biomass_covariance.csv')
    retained = allocation.groupby('Treatment').size().to_dict()
    audit = dict(dataset_doi='10.5285/dfe9a1ed-fc98-4abd-8040-6ee9f01cfcd0',
        candidate_article_doi='10.1016/j.scitotenv.2022.154893', final_article_read=False,
        biomass_cages=56, biomass_blocks=7, biomass_missing=0,
        insect_cages=len(allocation), insect_cage_weeks=len(insects),
        observed_insect_cages_by_treatment=retained,
        description_expected_cages_by_treatment={'Bl':7,'C':7,'G':6,'O':7,'P':7,'R':6,'Ub':7,'Uv':5},
        biomass_only_cages=joined.loc[joined._merge.eq('left_only'),'Cage'].tolist(),
        allocation_conflicts_in_shared_cages=0, extracted_biomass_contrasts=len(effects),
        covariance_rank=int(np.linalg.matrix_rank(cov)), covariance_method='paired_block_delta_lnRR; t6 intervals',
        independent_programmes=1, added_catalogue_effects=0,
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in source.glob('*') if f.suffix in ['.csv','.docx','.zip']},
        limitation='Sentinel plants exclude aphids; not final17week yield or a predator-removal contrast. Insect exclusions unresolved. UV photon flux differs from other colours.')
    (out/'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    print(effects[['endpoint','treatment','RR','ci_lower','ci_upper']].to_string(index=False))
    print('Biomass-only cages:', audit['biomass_only_cages'])


if __name__ == '__main__':
    main()
