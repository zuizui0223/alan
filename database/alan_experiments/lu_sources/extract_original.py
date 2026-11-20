"""Lu ECE 2026 original-sheet audit and paired-block descriptive effects.

Run from repository root. Does not populate approved meta-analysis tables.
"""
from pathlib import Path
import hashlib
import json
import math
import pandas as pd
import numpy as np

SRC = Path('../gbif_annual/alan_experiment_sources/lu2026')
OUT = Path(__file__).resolve().parent
EXPECTED = 'ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f'
hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in SRC.iterdir() if p.is_file()}
assert hashes['Date.xlsx'] == EXPECTED
df = pd.read_excel(SRC/'Date.xlsx', sheet_name='visitation')
cols = ['abundance','D_abundance','H_abundance','L_abundance','C_abundance','O_abundance']
keys = ['sampling_time','treatment','block','time','repeat']
assert len(df)==400 and not df[keys+cols].isna().any().any()
assert (df[cols]>=0).all().all()
df.insert(0,'source_excel_row',range(2,len(df)+2))
df.to_csv(OUT/'visitation_original_rows.csv',index=False)
duplicates=df[df.duplicated(keys,keep=False)]
duplicates.to_csv(OUT/'duplicate_session_keys.csv',index=False)
bad=df[df.abundance != df[cols[1:]].sum(axis=1)]
bad.to_csv(OUT/'guild_sum_discrepancies.csv',index=False)
df['three_orders_DHC']=df[['D_abundance','H_abundance','C_abundance']].sum(axis=1)
cols += ['three_orders_DHC']
groups=df.groupby(['sampling_time','treatment','block'])
block=groups[cols].sum().reset_index()
block['n_source_rows']=groups.size().values
block['nominal_survey_minutes_per_row']=20
block['implied_minutes_if_rows_distinct']=block.n_source_rows*20
block['planned_minutes_per_block_diel_arm']=400
block.to_csv(OUT/'block_totals.csv',index=False)
effects=[]
for diel in ['Diurnal','Nocturnal']:
 for metric in ['source_row_mean','source_block_total']:
  b=block[block.sampling_time==diel].copy()
  if metric=='source_row_mean':
   b[cols]=b[cols].div(b.n_source_rows,axis=0)
  for col in cols:
   a=b.pivot(index='block',columns='treatment',values=col).sort_index()
   x,y=a.Lit.to_numpy(),a.Unlit.to_numpy()
   n=len(x); assert n==5
   mx,my=x.mean(),y.mean(); sx,sy=x.std(ddof=1),y.std(ddof=1)
   cov=np.cov(x,y,ddof=1)[0,1]
   v=sx*sx/(n*mx*mx)+sy*sy/(n*my*my)-2*cov/(n*mx*my)
   status=('daytime_effort_verified_guild_inconsistency_flag_no_pool' if diel=='Diurnal' else 'holdout_pending_session_and_guild_audit')
   effects.append(dict(diel=diel,metric=metric,endpoint=col,n_pairs=n,mean_lit=mx,sd_lit=sx,mean_unlit=my,sd_unlit=sy,paired_sample_covariance=cov,lnRR=math.log(mx/my),paired_delta_variance=v,ratio=mx/my,status=status))
pd.DataFrame(effects).to_csv(OUT/'paired_block_descriptive_effects.csv',index=False)
# Cross-endpoint paired delta covariance, retaining shared plots/guilds.
day=block[block.sampling_time=='Diurnal'].copy()
day[cols]=day[cols].div(day.n_source_rows,axis=0)
wide=day.pivot(index='block',columns='treatment',values=cols)
influence=pd.DataFrame({c:wide[(c,'Lit')]/wide[(c,'Lit')].mean()-wide[(c,'Unlit')]/wide[(c,'Unlit')].mean() for c in cols})
(influence.cov()/len(influence)).to_csv(OUT/'daytime_effect_covariance.csv',index_label='endpoint')
audit=dict(doi='10.1002/ece3.74368',source_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles',hashes=hashes,rows=len(df),duplicate_key_rows=len(duplicates),unique_session_keys=len(df.drop_duplicates(keys)),guild_sum_discrepancy_rows=len(bad),block_row_counts=sorted(block.n_source_rows.unique().tolist()),totals_by_diel=df.groupby('sampling_time')[cols].sum().to_dict('index'),no_imputation=True,no_cross_sheet_merging=True,independent_blocks=5,nominal_session_minutes=20,meta_analysis_eligible=False)
(OUT/'extraction_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps(audit,indent=2))
print(pd.DataFrame(effects).query("diel=='Diurnal' and metric=='source_row_mean'")[['endpoint','mean_lit','mean_unlit','lnRR','paired_delta_variance']].to_string(index=False))
