"""Verify the saved V1 support audit against the downloaded workbook."""
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parent
CACHE=ROOT.parents[3]/'gbif_annual/alan_experiment_sources/grognuz2026'
d=pd.read_excel(CACHE/'23-02_OUTPUTPHENOLOGY.xlsx')
assert not d.duplicated(['STUDY.SITE','SPECIES','DAY']).any()
last=d.sort_values('DAY').groupby(['STUDY.SITE','SPECIES']).tail(1)
a=last.groupby(['SPECIES','TREATMENT']).agg(total_germinated=('SEEDLING.GERMINATED','sum'),sites_with_germination=('SEEDLING.GERMINATED',lambda v:(v>0).sum()))
b=pd.read_csv(ROOT/'germination_v1_support.csv').set_index(['SPECIES','TREATMENT'])
assert np.allclose(a,b.loc[a.index,a.columns])
passes=((a.total_germinated>=5)&(a.sites_with_germination>=3)).groupby(level=0).all()
assert passes.sum()==9 and len(d)==6336 and len(last)==176
sites=pd.read_csv(ROOT/'slug_programme_site_crosswalk.csv')
original=pd.read_csv(ROOT.parents[3]/'gbif_annual/alan_experiment_sources/slug2025/2023_STUDY SITES.csv',header=None)
assert set(sites['STUDY.SITE'])==set(original[0]) and len(sites)==16
print('Verified 176 site-species rows, nine qualifying codes and 16 shared site labels; no effect estimation.')
