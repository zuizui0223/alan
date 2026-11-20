"""All final-year material endpoints, preserving shared controls and zero means."""
from pathlib import Path
import hashlib
import zipfile
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/grassland2018'
archive=(CACHE/'archive_data.zip').read_bytes()
assert hashlib.md5(archive).hexdigest()=='6332d8b00cf937916248dc68c7cdd230'
with zipfile.ZipFile(CACHE/'archive_data.zip') as z:
    d=pd.read_excel(z.open('archive data/biomass_v3.xlsx'))
d.columns=d.columns.str.strip()
assert len(d)==120 and not d.duplicated(['Plot','Year']).any()
assert d.groupby('Plot').Light.nunique().eq(1).all() and d.notna().all().all()
material=[c for c in d if c not in ['Plot','Light','Year','total']]
discrepancy=d.total-d[material].sum(axis=1)
bad=d.loc[~np.isclose(discrepancy,0),['Plot','Light','Year']].copy()
bad['total_minus_components']=discrepancy[~np.isclose(discrepancy,0)]
bad.to_csv(OUT/'historical_total_discrepancies.csv',index=False)
last=d[d.Year==d.Year.max()].copy()
assert last.Year.eq(2015).all() and len(last)==30 and last.groupby('Light').size().eq(6).all()
assert np.allclose(last.total,last[material].sum(axis=1))
last['live_total']=last[[c for c in material if c!='Dead']].sum(axis=1)
last.to_csv(OUT/'final_plot_materials.csv',index=False)
endpoints=material+['total','live_total']
order=['CON','A','HIW','LIW','LIWT'];n=len(endpoints)
arms={a:last[last.Light==a].set_index('Plot')[endpoints] for a in order}
means=np.concatenate([arms[a].mean().to_numpy() for a in order]);cov=np.zeros((5*n,5*n))
for i,a in enumerate(order):cov[i*n:(i+1)*n,i*n:(i+1)*n]=arms[a].cov().to_numpy()/6
dose={'A':18.2,'HIW':29.6,'LIW':14.6,'LIWT':14.4}
effects=[];jac=[];ids=[]
for k,endpoint in enumerate(endpoints):
    for arm in order[1:]:
        i=order.index(arm)*n+k;j=k;lit,dark=means[i],means[j]
        eid=f'Bennie2018_{endpoint.replace(" ","_")}_{arm}'
        result=dict(effect_id=eid,study_id='Bennie2018',endpoint='grassland_biomass_'+endpoint,
            material=endpoint,treatment_code=arm,light_lux=dose[arm],
            mean_lit=lit,mean_dark=dark,n_lit=6,n_dark=6,
            lnRR=np.nan,sampling_variance=np.nan,lower=np.nan,upper=np.nan,df=np.nan,
            status='zero_arm_mean_no_lnRR',mean_units='source_workbook_mass_units',
            experiment_id='Lizard2012_2015_30plots',control_group_id='Lizard_dark6_2015',
            analysis_role='final_year_shared_dark_overlapping_material_endpoints',
            interpretation='2015 dry material; live_total excludes dead; total includes dead; absolute area scaling unresolved.')
        if lit>0 and dark>0:
            g=np.zeros(5*n);g[i]=1/lit;g[j]=-1/dark
            variance=float(g@cov@g);estimate=np.log(lit/dark)
            va=cov[i,i]/lit**2;vc=cov[j,j]/dark**2
            df=variance**2/(va**2/5+vc**2/5)
            half=t.ppf(.975,df)*np.sqrt(variance)
            result.update(lnRR=estimate,sampling_variance=variance,lower=estimate-half,upper=estimate+half,
                df=df,status='finite_final_year_lnRR')
            jac.append(g);ids.append(eid)
        effects.append(result)
e=pd.DataFrame(effects);e.to_csv(OUT/'all_final_contrasts.csv',index=False)
finite=e[e.status.eq('finite_final_year_lnRR')];finite.to_csv(OUT/'final_effects.csv',index=False)
v=np.asarray(jac)@cov@np.asarray(jac).T
assert np.allclose(np.diag(v),finite.sampling_variance)
assert np.linalg.eigvalsh(v).min()>-1e-10
pd.DataFrame([dict(effect_id_1=a,effect_id_2=b,sampling_covariance=v[i,j])
 for i,a in enumerate(ids) for j,b in enumerate(ids)]).to_csv(OUT/'final_covariance.csv',index=False)
audit=dict(raw_plot_years=120,final_year=2015,independent_plots=30,plots_per_arm=6,
    material_columns=len(material),all_contrasts=len(e),finite_effects=len(finite),
    zero_arm_contrasts=int(e.lnRR.isna().sum()),covariance_rank=int(np.linalg.matrix_rank(v)),
    historical_total_discrepancies=len(bad),final_totals_verified=True,
    cautions=['No pairing inferred from plot names; arms comprise randomized independent plots.',
              'Live biomass, dead matter and total including dead are different, overlapping outcomes.',
              'Common multiplicative mass/area conversion cancels from lnRR and delta covariance.',
              'Absolute source mass units not relabelled as g per square metre.',
              'All zero-arm endpoints retained; no pseudocount or significance selection.',
              'Final2015 differs from author longitudinal models and earlier years.'])
(OUT/'biomass_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
