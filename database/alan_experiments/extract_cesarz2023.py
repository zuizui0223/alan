"""Original Ecotron data: final phase, replicated extreme-dose contrast.

Separate extraction pending eligibility review; no pooling or wet/dry label.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'soil_sources'
OUT.mkdir(exist_ok=True)
src=ROOT.parents[2]/'gbif_annual/alan_experiment_sources/cesarz2023/rstb20220366_si_002.csv'
sha=hashlib.sha256(src.read_bytes()).hexdigest()
assert sha=='6fb6269122e721316e7e67c0c7f61dd6a606b3a8b946f33c727d48a471a0f18c'
d=pd.read_csv(src)
assert d.EcoUnit.nunique()==12 and len(d)==36
assert not d.duplicated(['EcoUnit','PhaseNr']).any()
assert d.groupby('EcoUnit').Lux.nunique().eq(1).all()
final=d[d.PhaseNr=='Phase 3'].copy()
lo,hi=final.Lux.min(),final.Lux.max()
endpoints={'bas':'soil_basal_respiration','cmic':'soil_microbial_biomass',
           'n.g':'total_nematode_density','BF_g':'bacterial_feeder_density',
           'HF_g':'fungal_feeder_density','PF_g':'plant_feeder_density',
           'OM_PR_g':'omnivore_predator_density','plant.biomass':'plant_biomass'}
arms={k:final[final.Lux==v].set_index('EcoUnit')[list(endpoints)] for k,v in [('low',lo),('high',hi)]}
assert all(len(v)==2 for v in arms.values())
assert all(v.notna().all().all() for v in arms.values())
means={k:v.mean() for k,v in arms.items()}
cov=sum(v.cov().to_numpy()/2/np.outer(means[k],means[k]) for k,v in arms.items())
rows=[]
for i,(col,name) in enumerate(endpoints.items()):
 x,y=arms['high'][col],arms['low'][col]
 assert x.mean()>0 and y.mean()>0
 rows.append(dict(effect_id='Cesarz2023_final_'+col,study_id='Cesarz2023',
    programme_id='iDiv_Ecotron2020_shared_with_Bucher2023',endpoint=name,
    mean_lit=x.mean(),sd_lit=x.std(),n_lit=2,mean_dark=y.mean(),sd_dark=y.std(),n_dark=2,
    high_lux=hi,low_lux=lo,lnRR=np.log(x.mean()/y.mean()),sampling_variance=cov[i,i],
    control='lowest_manipulated_light_not_dark',status='extracted_holdout_eligibility_and_overlap_review',
    final_phase='Phase 3',source_column=col))
pd.DataFrame(rows).to_csv(OUT/'cesarz2023_extreme_final_effects.csv',index=False)
pd.DataFrame([dict(effect_id_i=rows[i]['effect_id'],effect_id_j=rows[j]['effect_id'],covariance=cov[i,j])
              for i in range(len(rows)) for j in range(len(rows))]).to_csv(OUT/'cesarz2023_covariance.csv',index=False)
final[['EcoUnit','Date','PhaseNr','Lux','Light.type']+list(endpoints)].to_csv(OUT/'cesarz2023_final_units.csv',index=False)
manifest=dict(doi='10.1098/rstb.2022.0366',data_collection='10.6084/m9.figshare.c.6858054',
 download_url='https://ndownloader.figshare.com/files/42504285',sha256=sha,
 total_units=12,repeated_phases=3,extreme_units_per_arm=2,
 covariance_eigenvalue_min=float(np.linalg.eigvalsh(cov).min()),
 warning='8 outcomes from 4 chambers, not 8 studies; intermediate doses not in extreme contrast; indoor mesocosm scope and shared programme require review; no pooled estimate')
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(pd.DataFrame(rows)[['endpoint','lnRR','sampling_variance']].to_string(index=False))
