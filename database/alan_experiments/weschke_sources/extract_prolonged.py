"""Rebuild final night guild richness from author MaxN workbook, not abundance.

Keep all six guilds, missing lnRR for zero-arm means, and discrepancies with
the author derived table. No automatic admission to the main effect catalogue.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

OUT = Path(__file__).resolve().parent
RAW = OUT.parents[2].parent / 'gbif_annual/alan_experiment_sources/weschke2024'
manifest = json.loads((RAW/'manifest.json').read_text(encoding='utf-8'))
for item in manifest:
    assert hashlib.sha256((RAW/item['file']).read_bytes()).hexdigest() == item['sha256']
mapping = {'Piscivore':'Piscivores','Herbivore/detritivore':'Herbivores/detritivores',
           'Invertivore (mobile)':'Invertivores (mobile)',
           'Invertivore (sessile)':'Invertivores (sessile)',
           'Omnivore':'Omnivores','Planktivore':'Planktivores'}
d = pd.read_excel(RAW/'Raw_fish_video_survey_data.xlsx',sheet_name='Night Surveys',header=None)
assert set(d.iloc[7:,4]) == set(mapping)
rows = []
for c in range(5,d.shape[1]):
    if d.iloc[3,c] != 'Prolonged': continue
    counts = pd.to_numeric(d.iloc[7:,c],errors='raise')
    assert counts.notna().all() and counts.ge(0).all()
    row = dict(source_excel_column=c+1,site=str(d.iloc[0,c]).strip(),Treatment=d.iloc[2,c])
    row.update({dest:int(counts[d.iloc[7:,4].eq(src)].gt(0).sum()) for src,dest in mapping.items()})
    rows.append(row)
units = pd.DataFrame(rows)
assert len(units)==20 and units.site.is_unique
assert units.groupby('Treatment').size().to_dict()=={'ALAN':10,'Control':10}
units.to_csv(OUT/'prolonged_guild_units.csv',index=False)
derived = pd.read_csv(RAW/'Trophic_guild_Prolonged.csv')
assert len(derived)==20 and units.Treatment.equals(derived.Treatment)
diff = []
for i in range(20):
    for endpoint in mapping.values():
        if units.loc[i,endpoint] != derived.loc[i,endpoint]:
            diff.append(dict(derived_csv_row=i+2,raw_excel_site=units.loc[i,'site'],
                             endpoint=endpoint,raw_presence_count=int(units.loc[i,endpoint]),
                             author_derived_count=int(derived.loc[i,endpoint]),
                             row_link='same_order_candidate_not_independent_site_key'))
pd.DataFrame(diff).to_csv(OUT/'source_discrepancies.csv',index=False)
estimates=[]
finite=[]
for endpoint in mapping.values():
    x=units.loc[units.Treatment.eq('ALAN'),endpoint].to_numpy(float)
    y=units.loc[units.Treatment.eq('Control'),endpoint].to_numpy(float)
    mx,my=x.mean(),y.mean()
    row=dict(endpoint=endpoint,mean_lit=mx,mean_control=my,n_lit=10,n_control=10,
             sd_lit=x.std(ddof=1),sd_control=y.std(ddof=1),
             author_derived_mean_lit=derived.loc[derived.Treatment.eq('ALAN'),endpoint].mean(),
             author_derived_mean_control=derived.loc[derived.Treatment.eq('Control'),endpoint].mean(),
             status='zero_arm_mean_lnRR_undefined' if min(mx,my)==0 else 'conditional_raw_workbook_not_adopted')
    if min(mx,my)>0:
        parts=np.array([x.var(ddof=1)/(10*mx*mx),y.var(ddof=1)/(10*my*my)])
        v=parts.sum(); df=v*v/np.sum(parts*parts/9)
        lnrr=np.log(mx/my);half=t.ppf(.975,df)*np.sqrt(v)
        row.update(lnRR=lnrr,variance=v,welch_df=df,ratio=mx/my,
                   ci_low=np.exp(lnrr-half),ci_high=np.exp(lnrr+half))
        finite.append(endpoint)
    estimates.append(row)
pd.DataFrame(estimates).to_csv(OUT/'prolonged_guild_effects.csv',index=False)
cov=np.zeros((len(finite),len(finite)))
for arm in ['ALAN','Control']:
    a=units.loc[units.Treatment.eq(arm),finite]
    cov+=a.cov().to_numpy()/10/np.outer(a.mean(),a.mean())
assert np.linalg.eigvalsh(cov).min()>=-1e-12
pd.DataFrame(cov,index=finite,columns=finite).to_csv(OUT/'prolonged_guild_covariance.csv',index_label='endpoint')
audit=dict(source_doi='10.5281/zenodo.14266460',verified_small_files=len(manifest),
           source_rows_compared=120,discrepant_cells=len(diff),site_count=20,
           finite_lnRR=len(finite),undefined_zero_arm=6-len(finite),
           response='number of detected taxon labels per trophic guild; not individuals or predation',
           baseline_adjusted=False,final_eligibility=False,catalogue_added=0)
(OUT/'extraction_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(pd.DataFrame(estimates).to_string(index=False))
print(json.dumps(audit))
