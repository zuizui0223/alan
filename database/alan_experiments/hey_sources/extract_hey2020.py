"""Original Springer workbook; final biomass, plot-level delta-method covariance.

Run from repository root. Missing individual biomass is not zero-imputed.
The available-individual pot average is explicitly flagged for audit.
"""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
RAW = OUT.parents[2].parent / 'gbif_annual/alan_experiment_sources/hey2020'
src = RAW / 'supplement1.xlsx'
assert hashlib.sha256(src.read_bytes()).hexdigest() == '4593975a9df91fa79831d696bc084c2d11c6e3c06b643cfc3501a078372f74b7'
d = pd.read_excel(src, sheet_name='allocation&defense')
assert len(d) == 80 and d.ID.nunique() == 80
assert d.groupby('Plot').Light.nunique().eq(1).all()
d['biomass_g'] = d.DryWeightS + d.DryWeightR
units = d.groupby(['Plot', 'Pot', 'Light', 'Water', 'Density']).agg(
    biomass_g=('biomass_g', 'mean'), n_individuals=('ID', 'size'),
    n_observed=('biomass_g', 'count')).reset_index()
assert len(units) == 40
units['incomplete_pot'] = units.n_observed < units.n_individuals
assert units.incomplete_pot.sum() == 1
units.to_csv(OUT / 'hey2020_biomass_units.csv', index=False)
strata = [('NW', 'Low'), ('NW', 'High'), ('Water', 'Low'), ('Water', 'High')]
frames = {}
for light in ['Lit', 'Control']:
    x = units[units.Light.eq(light)].pivot(index='Plot', columns=['Water','Density'], values='biomass_g').reindex(columns=pd.MultiIndex.from_tuples(strata))
    assert x.shape == (5,4) and x.notna().all().all()
    frames[light] = x
means = {k:v.mean().to_numpy() for k,v in frames.items()}
cov = sum(x.cov().to_numpy()/len(x)/np.outer(means[k],means[k]) for k,x in frames.items())
effects = []
for i,(water,density) in enumerate(strata):
    a,b = frames['Lit'].iloc[:,i],frames['Control'].iloc[:,i]
    effects.append(dict(study_id='hey2020',experiment_id='hey2020_blandy',effect_id=f'hey2020_biomass_{water}_{density}',species='Asclepias syriaca',endpoint='final_total_dry_biomass_per_plant',time_weeks=4,unit='g',water_treatment=water,density=density,mean_light=a.mean(),sd_light=a.std(),n_light=5,mean_control=b.mean(),sd_control=b.std(),n_control=5,lnRR=np.log(a.mean()/b.mean()),variance=cov[i,i],independent_unit='lamp_plot',control_group_id=f'hey2020_control_{water}_{density}',dependency_block='hey2020_shared_plots',incomplete_individual_pot=bool(water=='NW' and density=='High'),extraction_status='original_raw_available_individual_pot_mean_requires_second_check',source_location='MOESM1 allocation&defense DryWeightS+DryWeightR',digitized=0))
e=pd.DataFrame(effects); e.to_csv(OUT/'hey2020_biomass_effects.csv',index=False)
pd.DataFrame([dict(effect_id_i=effects[i]['effect_id'],effect_id_j=effects[j]['effect_id'],covariance=cov[i,j],method='plot_matched_empirical_delta_method') for i in range(4) for j in range(4)]).to_csv(OUT/'hey2020_biomass_covariance.csv',index=False)
contrasts=[]
for i,j in [(2,0),(3,1)]:
    v=cov[i,i]+cov[j,j]-2*cov[i,j]
    contrasts.append(dict(density=strata[i][1],contrast='watered_minus_ambient_rain_lnRR',estimate=effects[i]['lnRR']-effects[j]['lnRR'],variance=v,se=np.sqrt(v),interpretation='within_study_water_addition_moderation_not_wet_vs_dry_habitat'))
pd.DataFrame(contrasts).to_csv(OUT/'hey2020_water_moderation.csv',index=False)
assert np.allclose(cov,cov.T) and np.linalg.eigvalsh(cov).min()>-1e-12
print(e[['water_treatment','density','lnRR','variance']].to_string(index=False))
print(pd.DataFrame(contrasts).to_string(index=False))
