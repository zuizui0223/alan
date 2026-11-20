"""Observed stem-diameter response before transfers; supplementary time window."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[2].parent/'gbif_annual/alan_experiment_sources/lockett2022'
metadata=json.loads((CACHE/'files.json').read_text(encoding='utf-8'))
entry=next(f for f in metadata if f['filename'].endswith('Data (trees).xlsx'))
path=CACHE/entry['filename']
assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['content_details']['sha256_hash']
g=pd.read_excel(path,'Growth')
g=g[g.week.isin([1,18])].copy()
assert len(g)==288 and not g.duplicated(['tree','week']).any()
assert g.diameter.notna().all() and g.diameter.gt(0).all()
assert g.groupby('plot1').light_treatment1.nunique().eq(1).all()
assert g.groupby('tree').plot1.nunique().eq(1).all()
assert g.groupby(['plot1','week']).size().eq(6).all()
g[['tree','plot1','light_treatment1','week','diameter','diamweek1','diam_change',
   'use4weeks1-18','dead_stunted1','dead_stunted11','dead_stunted18']].to_csv(OUT/'pretransfer_tree_rows.csv',index=False)
baseline=g[g.week.eq(1)].set_index('tree').diameter
end=g[g.week.eq(18)].set_index('tree')
assert np.allclose(end.diamweek1,baseline.loc[end.index])
assert np.allclose(end.diam_change,end.diameter/end.diamweek1)
assert end['use4weeks1-18'].eq('YES').sum()==138
# Preserve the full recorded initial cohort. Poor-condition flags mix death,
# damage, disease and stunting; no conversion to death or zero biomass.
plots=g.groupby(['plot1','light_treatment1','week']).diameter.mean().unstack('week').reset_index()
plots=plots.rename(columns={1:'diameter_week1',18:'diameter_week18'})
assert len(plots)==24 and plots.groupby('light_treatment1').size().eq(12).all()
plots.to_csv(OUT/'pretransfer_plot_means.csv',index=False)
value=0.;variance=0.;score_variance=0.;arms=[]
for treatment,sign in [('light',1),('dark',-1)]:
    x=plots[plots.light_treatment1.eq(treatment)][['diameter_week1','diameter_week18']].to_numpy()
    n=len(x);m=x.mean(axis=0);c=np.cov(x,rowvar=False,ddof=1)/n
    gradient=np.array([-1/m[0],1/m[1]])
    value+=sign*np.log(m[1]/m[0])
    variance+=float(gradient @ c @ gradient)
    score=(x-m) @ gradient/np.sqrt(n*(n-1))
    score_variance+=float(score@score)
    arms.append(dict(light=treatment,plots=n,trees=n*6,mean_week1=m[0],mean_week18=m[1],
        ratio_week18_to_week1=m[1]/m[0],variance_mean_week1=c[0,0],variance_mean_week18=c[1,1],
        covariance_means=c[0,1],units='archived_diameter_units_no_conversion'))
assert np.isclose(variance,score_variance,atol=1e-12)
half=float(t.ppf(.975,11)*np.sqrt(variance))
row=dict(effect_id='Lockett2022_pretransfer_recorded_diameter_growth_ratio',
    programme_id='Kongwak2018_2019_tree_psyllid',endpoint='recorded_stem_diameter_week18_to_week1_ratio',
    analysis_role='supplementary_pretransfer_not_final_harvest',log_ratio=float(value),
    sampling_variance=variance,ratio=float(np.exp(value)),approx_t11_lower=float(np.exp(value-half)),
    approx_t11_upper=float(np.exp(value+half)),n_lit_plots=12,n_dark_plots=12,
    trees=144,main_catalogue=False,interpretation='recorded_diameter_not_living_biomass_or_population_growth')
pd.DataFrame([row]).to_csv(OUT/'pretransfer_growth_effect.csv',index=False)
pd.DataFrame(arms).to_csv(OUT/'pretransfer_arm_summary.csv',index=False)
audit=dict(raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    original_trees=144,plots=24,pretransfer_week=18,baseline_week=1,
    recorded_poor_condition_flags_not_excluded=6,rows_imputed=0,
    author_flagged_138_analysis_reproduced=False,within_plot_time_covariance_retained=True,
    scalar_score_variance_agreement=True,main_catalogue_changed=False,
    final_harvest_label_conflict_resolved=False,
    limitations=['Supplementary pre-transfer window, not terminal endpoint replacement',
        'Full recorded cohort, not living biomass or randomized intent-to-treat estimate',
        'Systematic checkerboard lighting; no random lighting allocation or invented plot pairing',
        'Independent plot-error assumption, spatial dependence not estimated',
        'Approximate delta-method interval, t11; no predation or pollination mediation'])
(OUT/'pretransfer_growth_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(row,indent=2))
