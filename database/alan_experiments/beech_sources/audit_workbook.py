"""Preserve additional spectral treatment and link leaves to tree-level chemistry."""
from pathlib import Path
import hashlib
import json
import pandas as pd

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/beech2025'
raw=(CACHE/'raw_data_ceppi.xlsx').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='b7c0c7a81315db06951cabcfcb070f926236a42d1f05c0b93b65333cad88e78d'
c=pd.read_excel(CACHE/'raw_data_ceppi.xlsx',sheet_name='Chemistry-leaf_physilogy')
l=pd.read_excel(CACHE/'raw_data_ceppi.xlsx',sheet_name='Leaf_area_consumed')
treatment='Treatment (color_intensity)'
target={'control','4000_30','4000_50','4000_100'}
assert len(c)==56 and c.Tree_ID.is_unique and len(l)==855
assert set(c[treatment])==target
assert l.Perc_Herbivory.notna().all() and l.Perc_Herbivory.ge(0).all()
assert l.groupby('Tree_ID').size().eq(15).all()
assert l.groupby('Tree_ID')[['Site','Site_ID',treatment]].nunique().eq(1).all().all()
trees=l.groupby(['Site','Site_ID',treatment,'Tree_ID']).agg(
    measured_leaves=('Perc_Herbivory','size'),
    damaged_leaves=('Perc_Herbivory',lambda x:int(x.gt(0).sum())),
    mean_leaf_damage=('Perc_Herbivory','mean')).reset_index()
trees['chemistry_present']=trees.Tree_ID.isin(c.Tree_ID)
trees['stated_experiment_treatment']=trees[treatment].isin(target)
assert trees[trees.stated_experiment_treatment].chemistry_present.all()
assert set(trees.loc[~trees.stated_experiment_treatment,'Tree_ID'])=={'Alpthal63'}
trees.to_csv(OUT/'tree_support_and_spectral_scope.csv',index=False)
g=trees.groupby(['Site','Site_ID',treatment,'stated_experiment_treatment']).agg(
    trees=('Tree_ID','size'),leaves=('measured_leaves','sum'),damaged_leaves=('damaged_leaves','sum')).reset_index()
g.to_csv(OUT/'installation_support.csv',index=False)
assert g[g.stated_experiment_treatment].trees.eq(7).all()
audit=dict(data_doi='10.16904/envidat.628',sha256=hashlib.sha256(raw).hexdigest(),
    chemistry_rows=56,herbivory_rows=855,herbivory_tree_ids=57,
    target_tree_ids=56,target_leaf_rows=840,extra_tree='Alpthal63',
    extra_treatment=trees.loc[~trees.stated_experiment_treatment,treatment].iloc[0],
    extra_leaf_rows=15,target_installation_ids=8,sites=2,effects_added=0,
    caveats=['Dose treatment labels are instrument settings, not measured lux.',
             'Trees and leaves share installation; do not use n56 or n840 for lighting replication.',
             'Damage probability and severity conditional on damage are separate endpoints.',
             'Extra spectral treatment retained separately rather than silently discarded.'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
