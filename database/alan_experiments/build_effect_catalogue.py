"""Catalogue extracted effects without claiming approval or statistical pooling."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
rows=[]
def append(table, covariance, programme, habitat, role, mapping=None, select=None):
    d=pd.read_csv(ROOT/table)
    if select is not None:
        d=d.loc[select(d)].copy()
    if mapping:
        d=d.rename(columns=mapping)
    for _,r in d.iterrows():
        e={k:r.get(k,'') for k in ['effect_id','study_id','endpoint','mean_lit','mean_dark','n_lit','n_dark','lnRR','sampling_variance']}
        if not e['effect_id']:
            e['effect_id']='Lu2026_day_'+r.endpoint
        if not e['study_id']:
            e['study_id']='Lu2026'
        e.update(programme_id=programme,habitat=habitat,
                 analysis_role=r.get('analysis_role',role),source_table=table,
                 covariance_table=covariance,pooled=False,
                 approval_status='extracted_not_fully_adjudicated',
                 paired_units=r.get('paired_units',r.get('n_pairs','')))
        rows.append(e)
    assert (ROOT/covariance).exists()

append('effects_computed.csv','effect_covariance.csv','Parkinson2016','freshwater_littoral','primary')
append('plant_sources/anic2022_effects.csv','plant_sources/anic2022_covariance.csv',
       'Penryn2012_2016','artificial_sand_mesocosm_moisture_unclassified','factorial_stratum')
append('giavi_sources/giavi2021_effects.csv','giavi_sources/giavi2021_covariance.csv',
       'KnopGiavi_overlap_unresolved','meadow_moisture_unclassified','overlapping_daytime_endpoints')
append('hey_sources/hey2020_biomass_effects.csv','hey_sources/hey2020_biomass_covariance.csv',
       'Hey2017_Blandy','outdoor_pots_water_addition_not_habitat_class','factorial_stratum',
       dict(mean_light='mean_lit',mean_control='mean_dark',n_light='n_lit',n_control='n_dark',variance='sampling_variance'))
append('lu_sources/paired_block_descriptive_effects.csv','lu_sources/daytime_effect_covariance.csv',
       'Lu2024_Gannan_companion_overlap_unresolved','alpine_meadow_moisture_unclassified','overlapping_daytime_endpoints',
       dict(mean_unlit='mean_dark',paired_delta_variance='sampling_variance'),
       lambda d:(d.diel=='Diurnal') & (d.metric=='source_row_mean'))
d=pd.DataFrame(rows)
mask=d.study_id=='Lu2026'
d.loc[mask,'n_lit']=5
d.loc[mask,'n_dark']=5
assert d.effect_id.is_unique
assert np.isfinite(d.lnRR.astype(float)).all()
assert (d.sampling_variance.astype(float)>0).all()
assert np.allclose(np.log(d.mean_lit.astype(float)/d.mean_dark.astype(float)),d.lnRR.astype(float))
d.to_csv(ROOT/'effect_catalogue.csv',index=False)
summary=dict(extracted_effect_rows=len(d),study_labels=d.study_id.nunique(),
             programme_labels=d.programme_id.nunique(),
             counts_by_study=d.groupby('study_id').size().to_dict(),
             pooled_estimates=0,
             warning='Rows include overlapping endpoints, shared controls and supplementary time summaries; not independent studies. Eligibility and moisture moderators incomplete.')
(ROOT/'effect_catalogue_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
