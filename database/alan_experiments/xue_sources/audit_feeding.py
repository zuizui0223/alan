"""Separate structural blanks, missing observations and signed mass differences."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/xue2026'
raw=(CACHE/'Feeding_Preferences.xlsx').read_bytes()
assert hashlib.sha256(raw).hexdigest()=='50b04e4fff99e4429938eba360d6f587f40fedcd43ec76984f4ba401ef3e3e06'
d=pd.read_excel(CACHE/'Feeding_Preferences.xlsx')
d['Species']=d.Species.str.strip()
assert len(d)==100 and d.ID.is_unique and d.Rand.is_unique
assert d.groupby('Zone').Light.nunique().eq(1).all()
mapping={'V.natans':'Vallisneria_natans','H.verticillata':'Hydrilla_verticillata',
         'M.spicatum':'Myriophyllum_spicatum','M.aquaticum':'Myriophyllum_aquaticum'}
parts=[];structural=0
for suffix,species in mapping.items():
    assigned=d.Species.isin([species,'Mixed'])
    col='Consumed_biomass_'+suffix;res='Residual_biomass_'+suffix
    assert d.loc[~assigned,[col,res]].isna().all().all()
    structural+=int((~assigned).sum())
    g=d.loc[assigned,['ID','Zone','Light','Herbivore','Species',col,res]].copy()
    g=g.rename(columns={'Species':'offered_assemblage',col:'signed_mass_difference',res:'residual_mass'})
    g['plant']=species;g['assay']=np.where(g.offered_assemblage.eq('Mixed'),'choice','no_choice')
    g['missing_observation']=g.signed_mass_difference.isna()
    g['implied_initial_mass']=g.signed_mass_difference+g.residual_mass
    assert g.implied_initial_mass.dropna().max()-g.implied_initial_mass.dropna().min()<1e-10
    parts.append(g)
long=pd.concat(parts,ignore_index=True)
long.to_csv(OUT/'feeding_assigned_measurements.csv',index=False)
cells=long.groupby(['assay','Herbivore','plant','Light','Zone']).agg(
    assigned=('ID','size'),observed=('signed_mass_difference','count'),
    mean_signed_difference=('signed_mass_difference','mean'),
    negative=('signed_mass_difference',lambda x:int(x.lt(0).sum())))
cells.to_csv(OUT/'feeding_zone_support.csv')
coverage=cells.reset_index().groupby(['assay','Herbivore','plant','Light']).agg(
    observed_zones=('observed',lambda x:int(x.gt(0).sum())),assigned_zones=('Zone','nunique'))
coverage.to_csv(OUT/'feeding_arm_zone_coverage.csv')
audit=dict(buckets=100,structural_unassigned_cells=structural,assigned_plant_measurements=len(long),
    observed_measurements=int(long.signed_mass_difference.notna().sum()),
    missing_assigned_measurements=int(long.missing_observation.sum()),
    negative_signed_differences=int(long.signed_mass_difference.lt(0).sum()),
    light_zones=sorted(d.Zone.unique().tolist()),
    min_observed_zones_per_arm=int(coverage.observed_zones.min()),
    max_observed_zones_per_arm=int(coverage.observed_zones.max()),
    effects_added=0,reason='Signed initial-minus-residual mass; no clipping, pseudocount or lnRR. Uneven zone support and baseline uncertainty require outcome-specific extraction.')
assert structural==240 and len(long)==160 and audit['missing_assigned_measurements']==3
(OUT/'feeding_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
