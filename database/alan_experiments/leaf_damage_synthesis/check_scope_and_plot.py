"""Audit macro-cohort coverage and render the already-computed leaf-damage result."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EIV = REPO.parent/'gbif_annual/EIV_Europe.csv'
COHORT = REPO/'results/all_species_forest/wet_reproduction/species.csv'
eiv = pd.read_csv(EIV, encoding='latin1', index_col=0)
cohort = pd.read_csv(COHORT)
assert len(cohort)==363 and cohort.candidate.sum()==58
assert eiv.index.is_unique
rows=[]
for programme,name in [('Cieraad2023','Quercus robur'),('Cieraad2023','Vaccinium myrtillus'),
                       ('Ceppi2025','Fagus sylvatica'),('Guo2025','Psychotria asiatica'),
                       ('Guo2025','Aporosa')]:
    resolved = name!='Aporosa'
    matched = cohort[cohort.accepted.eq(name)] if resolved else cohort.iloc[:0]
    indicator = eiv.loc[name,'MOIST'] if resolved and name in eiv.index else None
    numeric = pd.to_numeric(indicator, errors='coerce')
    rows.append(dict(programme=programme,source_taxon=name,species_resolved=resolved,
                     eiv_exact_match=resolved and name in eiv.index,moist_as_reported=indicator,
                     moist_numeric=numeric,
                     candidate_label_matches=int(matched.candidate.sum()) if resolved else None,
                     comparison_label_matches=int((~matched.candidate).sum()) if resolved else None,
                     match_status='exact_accepted_name_only' if resolved else 'unresolved_not_absence'))
support=pd.DataFrame(rows)
support.to_csv(HERE/'plant_scope.csv',index=False)
audit=dict(experimental_taxon_units=len(rows),resolved_binomials=int(support.species_resolved.sum()),
           unresolved_genus_units=int((~support.species_resolved).sum()),
           resolved_binomials_matching_wet_cohort=int(((support.candidate_label_matches+support.comparison_label_matches)>0).sum()),
           numeric_moisture_taxa=int(support.moist_numeric.notna().sum()),
           candidate_labels=58,comparison_labels=305,accepted_cohort_taxa=int(cohort.accepted.nunique()),
           scope='exact accepted-binomial matches; no synonym resolution or missing-to-absence inference',
           directly_tests_wet_decline_candidates=False,
           input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [EIV,COHORT]},
           conclusion='General experimental leaf-damage evidence; cannot reject or establish a wet-plant mechanism')
(HERE/'plant_scope_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')

d=pd.read_csv(HERE/'programme_estimates.csv')
p=pd.read_csv(HERE/'pooled_estimate.csv').iloc[0]
fig,ax=plt.subplots(figsize=(10,5.5))
fig.subplots_adjust(left=.26,right=.72,top=.79,bottom=.32)
labels=['Netherlands: Cieraad 2023','Switzerland: Ceppi 2025','Hong Kong: Guo 2025']
for i,r in d.iterrows():
    y=3-i
    ax.plot([r.lower,r.upper],[y,y],color='#436880',lw=1.8)
    ax.scatter(r.ratio,y,s=42,color='#436880',zorder=3)
    ax.text(1.03,y,f'{r.ratio:.2f} [{r.lower:.3f}, {r.upper:.2f}]',
            transform=ax.get_yaxis_transform(),va='center',fontsize=10)
ax.plot([p.lower,p.upper],[0,0],color='#a3422d',lw=2)
ax.scatter(p.ratio,0,marker='D',s=65,color='#a3422d',zorder=3)
ax.text(1.03,0,f'{p.ratio:.2f} [{p.lower:.3f}, {p.upper:.2f}]',
        transform=ax.get_yaxis_transform(),va='center',fontsize=10)
ax.axvline(1,color='0.45',ls='--',lw=1)
ax.axhline(.5,color='0.85',lw=.8)
ax.set_xscale('log');ax.set_xlim(.005,300);ax.set_ylim(-.6,3.6)
ax.set_xticks([.01,.1,1,10,100]);ax.set_xticklabels(['0.01','0.1','1','10','100'])
ax.set_yticks([3,2,1,0]);ax.set_yticklabels([f'{label}\n{int(n)} cluster pairs' for label,n in zip(labels,d.clusters)]+['Exploratory pooled\n3 programmes'])
ax.set_xlabel('Leaf-damage response ratio (illuminated / dark); logarithmic scale')
ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
fig.suptitle('Terminal leaf damage: three independent experimental programmes',fontsize=13,y=.96)
fig.text(.26,.85,'Estimate and approximate 95% interval',fontsize=10)
fig.text(.03,.055,'11 dependent contrasts; species and lighting menus averaged within programmes.\nPaired delta t intervals within programmes; REML + modified Knapp–Hartung t2 for pooling.\nNo verified overlap with the wet-decline cohort among four resolved plant names.',fontsize=9)
target=REPO/'manuscript/biotic_synthesis/figures';target.mkdir(exist_ok=True)
for suffix in ['png','pdf','svg']:fig.savefig(target/f'leaf_damage_meta.{suffix}',dpi=180)
plt.close(fig)
print(json.dumps(audit,indent=2))
