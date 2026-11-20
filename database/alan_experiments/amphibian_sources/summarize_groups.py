"""Descriptive shared-light group summaries; no adoption or inferential effects."""
from pathlib import Path
import hashlib,json
import pandas as pd
OUT=Path(__file__).resolve().parent
SRC=OUT.parents[3]/'gbif_annual/alan_experiment_sources/amphibian2018/original.xlsx'
assert hashlib.md5(SRC.read_bytes()).hexdigest()=='4c54998fd893d7a766b34ffb9ea2af8b'
lookup=pd.read_excel(SRC,'AmphibianPipeSample')[['Mesocosm','Block','LightTreat']]
assert lookup.Mesocosm.is_unique and len(lookup)==36
rows=[]
for sheet,variables in [('InvertPipeSample',['PredatorTotal','NonPredatorTotal']),
                        ('toadmesocosmwaterquality',['Leaf2July','Water2July'])]:
    original=pd.read_excel(SRC,sheet)
    assert original.Mesocosm.is_unique and set(original.Mesocosm)==set(lookup.Mesocosm)
    frame=original.merge(lookup,on='Mesocosm',validate='one_to_one',suffixes=('_source',''))
    assert frame.LightTreat_source.replace({'Light':'Present','No Light':'Absent'}).eq(frame.LightTreat).all()
    for block,g in frame.groupby('Block'):
        assert len(g)==4 and len(g[['LidTreat','ToadTreat']].drop_duplicates())==4
        assert g.LidTreat.nunique()==2 and g.ToadTreat.nunique()==2
        assert g.LightTreat.nunique()==1
        for variable in variables:
            assert g[variable].notna().all() and g[variable].ge(0).all()
            rows.append(dict(sheet=sheet,endpoint=variable,block=block,
                             light=g.LightTreat.iloc[0],mesocosms=4,mean=g[variable].mean(),
                             mesocosm_ids=';'.join(map(str,sorted(g.Mesocosm)))))
groups=pd.DataFrame(rows)
groups.to_csv(OUT/'descriptive_group_means.csv',index=False)
summary=[]
for ep,g in groups.groupby('endpoint'):
    lit=g[g.light=='Present']['mean'];dark=g[g.light=='Absent']['mean']
    assert len(lit)==4 and len(dark)==5 and dark.mean()>0
    summary.append(dict(endpoint=ep,lit_groups=4,dark_groups=5,lit_mean=lit.mean(),
                        dark_mean=dark.mean(),descriptive_ratio=lit.mean()/dark.mean(),
                        meta_adopted=False,interval_estimated=False))
pd.DataFrame(summary).to_csv(OUT/'descriptive_arm_means.csv',index=False)
print(pd.DataFrame(summary).to_string(index=False))
