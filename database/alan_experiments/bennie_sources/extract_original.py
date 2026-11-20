"""Bennie2015 author workbook. Run from repo root; no meta-table mutation."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
O=Path(__file__).resolve().parent
S=Path('../gbif_annual/alan_experiment_sources/bennie2015')
file=S/'rstb20140131supp1.xlsx'
sha=hashlib.sha256(file.read_bytes()).hexdigest()
assert sha=='009ea860f7493cb3ccd8eb0fa12503702f30eed1c38b2b570c4e0e6eb695f412'
allrows=[]
for sheet,endpoint in [('Acyrthosiphon pisum','aphid_count_3min'),('Lotus pedunculatus','lotus_flowerhead_count')]:
 d=pd.read_excel(file,sheet_name=sheet,header=1)
 d.columns=['plot','herbivores','predators','light']+list(d.columns[4:])
 assert len(d)==54 and d['plot'].nunique()==54
 d['trophic']=np.where(d.herbivores=='no','T1',np.where(d.predators=='no','T2','T3'))
 d['source_excel_row']=np.arange(3,57)
 dates=[c for c in d.columns if isinstance(c,pd.Timestamp) or hasattr(c,'year')]
 long=d.melt(id_vars=['plot','herbivores','predators','light','trophic','source_excel_row'],value_vars=dates,var_name='date',value_name='value')
 long['endpoint']=endpoint;long['source_sheet']=sheet;long['date']=pd.to_datetime(long.date)
 allrows.append(long)
raw=pd.concat(allrows,ignore_index=True)
assert raw.query("endpoint=='aphid_count_3min' and trophic!='T1'").value.notna().all()
raw.to_csv(O/'repeated_original_units.csv',index=False)
final=raw[raw.date==raw.groupby('endpoint').date.transform('max')].copy()
final.to_csv(O/'final_original_units.csv',index=False)
summ=final.groupby(['endpoint','date','trophic','light']).value.agg(['count','mean','std']).reset_index()
summ.to_csv(O/'final_arm_summaries.csv',index=False)
effects=[]; components={}
for (endpoint,trophic),g in final.groupby(['endpoint','trophic']):
 if endpoint=='aphid_count_3min' and trophic=='T1':continue
 dark=g[g.light=='Dark'].set_index('plot').value
 assert len(dark)==6 and dark.notna().all()
 for light in ['Amber','White']:
  lit=g[g.light==light].set_index('plot').value
  assert len(lit)==6 and lit.notna().all()
  mx,my=lit.mean(),dark.mean();valid=mx>0 and my>0
  ident=f'Bennie2015_{endpoint}_{trophic}_{light}'
  effects.append(dict(effect_id=ident,programme_id='Penryn2012_2016',experiment_id='Penryn2012_2016',endpoint=endpoint,timepoint=str(g.date.iloc[0].date()),trophic=trophic,light=light,control_group_id=f'Penryn_{trophic}_dark_2013_{endpoint}',mean_lit=mx,sd_lit=lit.std(ddof=1),n_lit=6,mean_dark=my,sd_dark=dark.std(ddof=1),n_dark=6,lnRR=np.log(mx/my) if valid else np.nan,sampling_variance=lit.var(ddof=1)/6/mx**2+dark.var(ddof=1)/6/my**2 if valid else np.nan,mean_difference=mx-my,mean_difference_variance=(lit.var(ddof=1)+dark.var(ddof=1))/6,status='extracted_pending_independent_verification' if valid else 'undefined_lnRR_zero_mean_no_pseudocount'))
  if valid:components[ident]={(trophic,light):lit/mx,(trophic,'Dark'):-dark/my}
pd.DataFrame(effects).to_csv(O/'final_effects.csv',index=False)
# Same mesocosm IDs independently verify this is the same programme as Anic2022.
units=final.query("endpoint=='lotus_flowerhead_count'")[['plot','trophic','light']]
anic=pd.read_csv('database/alan_experiments/plant_sources/anic2022_community_units.csv')
cross=units.merge(anic,left_on='plot',right_on='mesocosm',suffixes=('_2013','_2016'),validate='one_to_one')
assert len(cross)==54 and (cross.trophic_2013==cross.trophic_2016).all()
assert (cross.light_2013.map({'Dark':'L1','Amber':'L2','White':'L3'})==cross.light_2016).all()
cross.to_csv(O/'anic2022_same_mesocosm_crosswalk.csv',index=False)
for trophic,g in cross.groupby('trophic_2013'):
 for name,code in [('Amber','L2'),('White','L3')]:
  lit=g[g.light_2013==name].set_index('plot').total_biomass_g
  dark=g[g.light_2013=='Dark'].set_index('plot').total_biomass_g
  components[f'Anic2022_community_{trophic}_{code}']={(trophic,name):lit/lit.mean(),(trophic,'Dark'):-dark/dark.mean()}
def cov(a,b):
 v=0.
 for key in set(a)&set(b):
  x,y=a[key].align(b[key],join='inner');assert len(x)==6
  v+=x.cov(y)/6
 return v
covs=[dict(effect_id_1=a,effect_id_2=b,sampling_covariance=cov(x,y),method='same_mesocosm_endpoint_covariance_plus_shared_controls',programme_id='Penryn2012_2016') for a,x in components.items() for b,y in components.items()]
pd.DataFrame(covs).to_csv(O/'programme_covariance_including_anic2022.csv',index=False)
interactions=[]
edf=pd.DataFrame(effects).set_index('effect_id')
for light in ['Amber','White']:
 a=f'Bennie2015_aphid_count_3min_T3_{light}';b=f'Bennie2015_aphid_count_3min_T2_{light}'
 interactions.append(dict(endpoint='aphid_count_3min',light=light,contrast='ALAN lnRR with predators minus ALAN lnRR without predators',log_ratio_of_ratios=edf.loc[a,'lnRR']-edf.loc[b,'lnRR'],variance=cov(components[a],components[a])+cov(components[b],components[b])-2*cov(components[a],components[b]),n_per_cell=6,timepoint='2013-09-12',programme_id='Penryn2012_2016',status='single_final_date_not_published_repeated_model'))
pd.DataFrame(interactions).to_csv(O/'final_predator_interactions.csv',index=False)
ic=[]
for light1 in ['Amber','White']:
 for light2 in ['Amber','White']:
  value=0.
  for t1,s1 in [('T3',1),('T2',-1)]:
   for t2,s2 in [('T3',1),('T2',-1)]:
    value+=s1*s2*cov(components[f'Bennie2015_aphid_count_3min_{t1}_{light1}'],components[f'Bennie2015_aphid_count_3min_{t2}_{light2}'])
  ic.append(dict(light1=light1,light2=light2,sampling_covariance=value))
pd.DataFrame(ic).to_csv(O/'predator_interaction_covariance.csv',index=False)
manifest=dict(paper_doi='10.1098/rstb.2014.0131',source_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4375371/supplementaryFiles',file=file.name,sha256=sha,programme_id='Penryn2012_2016',mesocosms=54,matched_anic_mesocosms=54,repeated_rows=len(raw),blank_aphid_plant_only_not_zero=int(raw.value.isna().sum()),programme_independent_count=1,no_pseudocount=True,independent_verification='pending')
(O/'source_audit.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print(pd.DataFrame(effects)[['endpoint','trophic','light','mean_lit','mean_dark','lnRR','sampling_variance']].to_string(index=False));print(pd.DataFrame(interactions).to_string(index=False))
