"""Inspect deposited cohort/cage support without treating missing observations as death."""
from pathlib import Path
import pandas as pd
import hashlib,json
OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/gurholt2026'
raw=(CACHE/'original.csv').read_bytes()
assert hashlib.md5(raw).hexdigest()=='e666e41643f7357c22fbea8345d83574'
d=pd.read_csv(CACHE/'original.csv')
assert len(d)==294 and not d.duplicated(['Cohort','Individual']).any()
d['cage_prefix']=d.Cage.str[:2]
conflict=d[d.Treatment!=d.cage_prefix]
conflict[['Cohort','Individual','Treatment','Cage','Species']].to_csv(OUT/'treatment_conflicts.csv',index=False)
endpoints=['Initial Mass','Day 7 Mass','Adult wet mass','Adult dry mass','Date of pupation','Date of eclosion','Date of death']
rows=[]
for keys,g in d.groupby(['Cohort','Cage','Species']):
 r=dict(zip(['cohort','cage','host'],keys));r['rows']=len(g)
 r['treatment_labels']=';'.join(sorted(g.Treatment.unique()))
 for ep in endpoints:r[ep+'_observed']=int(g[ep].notna().sum())
 rows.append(r)
pd.DataFrame(rows).to_csv(OUT/'cohort_cage_support.csv',index=False)
observed=d.groupby('Cohort')[endpoints].count()
observed.to_csv(OUT/'cohort_endpoint_support.csv')
empty=['Initial Mass','Day 7 Mass','Date of pupation','Date of eclosion','Date of death','Adult wet mass']
noobs=d[d[empty].isna().all(axis=1)]
noobs[['Cohort','Individual','Treatment','Cage','Species']].to_csv(OUT/'no_primary_observation_rows.csv',index=False)
paired=d[['Adult wet mass','Adult wet mass.1']].dropna()
ratio=paired['Adult wet mass.1']/paired['Adult wet mass']
a=dict(rows=len(d),individual_labels=d.Individual.nunique(),cohort_individual_key_unique=True,
       cohort_cage_combinations=d[['Cohort','Cage']].drop_duplicates().shape[0],
       reused_cage_labels=d.Cage.nunique(),treatment_conflicts=len(conflict),
       rows_missing_all_six_selected_observations=len(noobs),
       adult_mass_paired_rows=len(paired),adult_mass_scale_ratio_min=float(ratio.min()),
       adult_mass_scale_ratio_max=float(ratio.max()),full_methods_read=True,
       missing_as_death=False,effects_extracted=0)
(OUT/'audit.json').write_text(json.dumps(a,indent=2)+'\n')
print(json.dumps(a,indent=2))
