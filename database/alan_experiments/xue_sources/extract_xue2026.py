"""Extract original ScienceDB V1 data; no bucket-level pseudo-replication.

Only monoculture provides Zone identifiers verified against primary Methods.
Community data are descriptive and gated from meta-analysis pending lamp layout.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
RAW = OUT.parents[2].parent / 'gbif_annual/alan_experiment_sources/xue2026'
files = json.loads((RAW/'files_root_V1.json').read_text(encoding='utf8'))['data']
manifest=[]
for f in files:
    p=RAW/f['fileName']
    if not p.exists():
        continue
    b=p.read_bytes()
    assert hashlib.md5(b).hexdigest()==f['md5'], f['fileName']
    manifest.append(dict(filename=f['fileName'],file_id=f['id'],md5=f['md5'],sha256=hashlib.sha256(b).hexdigest(),data_doi='10.57760/sciencedb.35715',version='V1',url='https://www.scidb.cn/api/sdb-download-service/downloadFileMole?fileId='+f['id']))
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.csv',index=False)
d=pd.read_excel(RAW/'Monoculture_Growth.xlsx')
d['Species']=d.Species.str.strip()
assert len(d)==77 and d.ID.nunique()==77 and d.Rand.nunique()==77
assert d.groupby('Zone').Light.nunique().eq(1).all()
units=d.groupby(['Light','Zone','Species']).agg(biomass_g=('Total_biomass','mean'),n_surviving_cylinders=('ID','size')).reset_index()
units.to_csv(OUT/'monoculture_zone_units.csv',index=False)
species=sorted(d.Species.unique())
frames={light:units[units.Light.eq(light)].pivot(index='Zone',columns='Species',values='biomass_g').reindex(columns=species) for light in ['ALAN','No ALAN']}
assert all(x.shape==(3,4) and x.notna().all().all() for x in frames.values())
means={k:x.mean().to_numpy() for k,x in frames.items()}
cov=sum(x.cov().to_numpy()/len(x)/np.outer(means[k],means[k]) for k,x in frames.items())
effects=[]
for i,sp in enumerate(species):
    a,b=frames['ALAN'].iloc[:,i],frames['No ALAN'].iloc[:,i]
    effects.append(dict(study_id='xue2026',experiment_id='xue2026_monoculture',effect_id='xue2026_monoculture_'+sp,species=sp.replace('_',' '),response='total_dry_biomass',unit='g',endpoint_days=50,mean_light=a.mean(),sd_light=a.std(),n_light=3,mean_control=b.mean(),sd_control=b.std(),n_control=3,lnRR=np.log(a.mean()/b.mean()),variance=cov[i,i],independent_unit='lighting_zone',aggregation='equal_zone_weight_surviving_cylinder_means',dependency_block='xue2026_monoculture_shared_zones',habitat='freshwater_outdoor_mesocosm',extraction_status='original_data_extracted_pending_independent_check',survivor_selection=True,digitized=0))
pd.DataFrame(effects).to_csv(OUT/'monoculture_effects.csv',index=False)
pd.DataFrame([dict(effect_id_i=effects[i]['effect_id'],effect_id_j=effects[j]['effect_id'],covariance=cov[i,j],method='matched_zone_empirical_delta_method') for i in range(4) for j in range(4)]).to_csv(OUT/'monoculture_covariance.csv',index=False)
assert np.allclose(cov,cov.T) and np.linalg.eigvalsh(cov).min()>-1e-12
c=pd.read_excel(RAW/'Community.xlsx')
assert len(c)==57 and c.ID.nunique()==57 and c.Rand.nunique()==57
assert np.allclose(c.Invasive_proportional_biomass,c.Invasive_total_biomass/(c.Invasive_total_biomass+c.Native_total_biomass))
out=[]
for (light,herb),g in c.groupby(['Light','Herbivore']):
    for endpoint in ['Native_total_biomass','Invasive_total_biomass']:
        out.append(dict(light=light,herbivore=herb,endpoint=endpoint,mean_g=g[endpoint].mean(),sd_across_buckets_g=g[endpoint].std(),observed_buckets=len(g),independent_lighting_n='unresolved',eligible_for_meta=False,gate='community_lamp_layout_and_zone_mapping_missing',selection_note='three_omitted_mesocosms_all_ALAN_native_snail'))
pd.DataFrame(out).to_csv(OUT/'community_descriptive_only.csv',index=False)
print(pd.DataFrame(effects)[['species','lnRR','variance']].to_string(index=False))
print('Verified: 77 cylinders; 6 zones; 24 zone-species means; 57 community buckets; 5 source MD5 checks.')
