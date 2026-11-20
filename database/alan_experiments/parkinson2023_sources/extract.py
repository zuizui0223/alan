"""Original Figshare trap data; raw captures and selected-family biomass kept distinct."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
CACHE = ROOT.parent/'gbif_annual/alan_experiment_sources/parkinson2023'
FILE = CACHE/'rstb20220361_si_004.xlsx'
x = pd.ExcelFile(FILE)
raw = pd.read_excel(x, 'raw site data')
summary = pd.read_excel(x, 'overall abundance ')
top = pd.read_excel(x, 'top5families')
bio = pd.read_excel(x, 'total biomass')
nmds = pd.read_excel(x, 'June3NMDS')
aquatic = set('caenidae hydroptillidae tipulidae chironomidae coenagrionidae leptoceridae phryganeidae polycentropodidae simuliidae rhyacophilidae corydalidae culicidae gerrdiae notonectidae psychomiidae helicopsyche hydropsychidae'.split())
terrestrial = set('miturgidae cicadellidae noctuoidae dolichiopodidae muscidae tettigoniidae ephydridae formicidae hahnidae arctiinae erebidae cosmopterigidae rhopalosomatidae miridae diapriidae geometridae scarabaeidae tetragnathidae'.split())
spiders = set('miturgidae hahnidae tetragnathidae'.split())
five = ['caenidae','hydroptillidae','leptoceridae','cicadellidae','chironomidae']
cols = ['number caenidae','number hydroptillidae','number leptoceridae ','number cicadellidae','num chir']
starts = list(raw.index[raw.habitat.notna()])
assert len(starts) == 30
units=[]; taxa=[]; discrepancies=[]
for k,start in enumerate(starts):
    end = starts[k+1] if k+1<len(starts) else len(raw)
    block=raw.iloc[start:end].dropna(subset=['family','abundance']).copy()
    block['family']=block.family.str.strip().str.lower()
    counts=block.groupby('family').abundance.sum().to_dict()
    assert set(counts) <= aquatic|terrestrial, set(counts)-(aquatic|terrestrial)
    habitat={'aquatic':'littoral','terrestrial':'riparian'}[raw.loc[start,'habitat']]
    colour=raw.loc[start,'color'].strip()
    assert habitat == summary.iloc[k].habitat == bio.iloc[k].habitat == top.iloc[k].habitat
    assert colour == summary.iloc[k]['color'] == bio.iloc[k]['color'] == top.iloc[k]['color'] == nmds.iloc[k]['color']
    assert raw.loc[start,'habitat'] == nmds.iloc[k].habitat
    for family,col in zip(five,cols):
        if counts.get(family,0)!=top.iloc[k][col]:
            discrepancies.append(dict(source_row=start+2,summary_row=k+2,field=family,raw_value=counts.get(family,0),summary_value=top.iloc[k][col]))
    if sum(counts.get(f,0) for f in five)!=summary.iloc[k].abundance:
        discrepancies.append(dict(source_row=start+2,summary_row=k+2,field='five_family_count',raw_value=sum(counts.get(f,0) for f in five),summary_value=summary.iloc[k].abundance))
    # Source row key is primary; do not rely on cross-sheet IDs given discrepancies.
    unit=f'raw_row_{start+2}'
    u=dict(study_id='Parkinson2023',programme_id='SevenPonds2022',unit_id=unit,
           source_plot_id_unverified=nmds.iloc[k]['Plot ID'],source_excel_row=start+2,
           summary_excel_row=k+2,habitat=habitat,colour=colour,date='2022-06-03',
           all_invertebrate_capture=sum(counts.values()),
           aquatic_origin_capture=sum(v for f,v in counts.items() if f in aquatic),
           terrestrial_origin_capture=sum(v for f,v in counts.items() if f in terrestrial),
           spider_capture=sum(v for f,v in counts.items() if f in spiders),
           cicadellidae_capture=counts.get('cicadellidae',0),
           five_family_capture=sum(counts.get(f,0) for f in five))
    units.append(u)
    for f,n in counts.items():
        taxa.append(dict(unit_id=unit,family_source_spelling=f,count=n,
                         origin_as_published_figure3='aquatic' if f in aquatic else 'terrestrial',
                         spider=f in spiders))
u=pd.DataFrame(units); pd.DataFrame(taxa).to_csv(OUT/'taxon_counts.csv',index=False)
assert u.all_invertebrate_capture.sum()==4699
assert u.five_family_capture.sum()==4475
pd.DataFrame(discrepancies).to_csv(OUT/'cross_sheet_discrepancies.csv',index=False)
assert (u.groupby(['habitat','colour']).size()==3).all()
u.to_csv(OUT/'trap_units.csv',index=False)
endpoints=['all_invertebrate_capture','aquatic_origin_capture','terrestrial_origin_capture',
           'spider_capture','cicadellidae_capture']
effects=[]
for habitat in ['littoral','riparian']:
    control=u[(u.habitat==habitat)&(u.colour=='dark')]
    for colour in ['green','purple','red','white']:
        lit=u[(u.habitat==habitat)&(u.colour==colour)]
        for endpoint in endpoints:
            a,b=lit[endpoint],control[endpoint]
            valid=a.mean()>0 and b.mean()>0
            effects.append(dict(effect_id=f'P2023_{habitat}_{colour}_{endpoint}',study_id='Parkinson2023',
                programme_id='SevenPonds2022',experiment_id='SevenPonds_one_night_20220603',
                control_group_id=f'P2023_{habitat}_dark',habitat=habitat,colour=colour,endpoint=endpoint,
                n_lit=len(a),n_control=len(b),mean_lit=a.mean(),mean_control=b.mean(),
                sd_lit=a.std(ddof=1),sd_control=b.std(ddof=1),
                lnRR=np.log(a.mean()/b.mean()) if valid else np.nan,
                variance=a.var(ddof=1)/len(a)/a.mean()**2+b.var(ddof=1)/len(b)/b.mean()**2 if valid else np.nan,
                status='extractable_delta_lnRR' if valid else 'zero_mean_no_lnRR',
                interpretation='trap_capture_not_population_size' if endpoint!='five_family_estimated_biomass_mg' else 'allometry_derived_selected_five_families'))
e=pd.DataFrame(effects);e.to_csv(OUT/'effects.csv',index=False)
v=e[e.status=='extractable_delta_lnRR'];cov=[]
for _,a in v.iterrows():
    for _,b in v.iterrows():
        value=0.
        if a.habitat==b.habitat:
            for colour in ['dark']+([a.colour] if a.colour==b.colour else []):
                g=u[(u.habitat==a.habitat)&(u.colour==colour)]
                value+=g[[a.endpoint,b.endpoint]].cov().iloc[0,1]/len(g)/g[a.endpoint].mean()/g[b.endpoint].mean()
        cov.append(dict(effect_id_1=a.effect_id,effect_id_2=b.effect_id,covariance=value))
c=pd.DataFrame(cov);c.to_csv(OUT/'covariance.csv',index=False)
matrix=c.pivot(index='effect_id_1',columns='effect_id_2',values='covariance').reindex(index=v.effect_id,columns=v.effect_id).values
assert np.allclose(np.diag(matrix),v.variance)
assert np.linalg.eigvalsh(matrix).min()>-1e-9
mods=[]
for colour in ['green','purple','red','white']:
    for endpoint in endpoints:
        g=e[(e.colour==colour)&(e.endpoint==endpoint)].set_index('habitat')
        valid=g.lnRR.notna().all()
        mods.append(dict(colour=colour,endpoint=endpoint,
            lnRR_littoral_minus_riparian=g.loc['littoral','lnRR']-g.loc['riparian','lnRR'] if valid else np.nan,
            variance=g.variance.sum() if valid else np.nan,
            status='within_pond_placement_contrast_not_wet_dry_biome' if valid else 'zero_mean_gate'))
pd.DataFrame(mods).to_csv(OUT/'habitat_moderation.csv',index=False)
mc=[]
validmods=[m for m in mods if np.isfinite(m['lnRR_littoral_minus_riparian'])]
lookup={(r.effect_id_1,r.effect_id_2):r.covariance for r in c.itertuples()}
for a in validmods:
    for b in validmods:
        cv=0.
        for h1,s1 in [('littoral',1),('riparian',-1)]:
            for h2,s2 in [('littoral',1),('riparian',-1)]:
                cv+=s1*s2*lookup[(f"P2023_{h1}_{a['colour']}_{a['endpoint']}",f"P2023_{h2}_{b['colour']}_{b['endpoint']}")]
        mc.append(dict(moderator_id_1=a['colour']+'_'+a['endpoint'],moderator_id_2=b['colour']+'_'+b['endpoint'],covariance=cv))
pd.DataFrame(mc).to_csv(OUT/'habitat_moderation_covariance.csv',index=False)
manifest=dict(doi='10.1098/rstb.2022.0361',collection='10.6084/m9.figshare.c.6856475',
    source_file=FILE.name,sha256=hashlib.sha256(FILE.read_bytes()).hexdigest(),
    retrieval_date='2026-11-20',source_api='https://api.figshare.com/v2/articles/24212174',
    raw_total=4699,raw_five_family_total=int(u.five_family_capture.sum()),summary_five_family_total=4465,top5_sheet_total=4448,trap_units=30,
    contrasts=len(e),valid_lnRR=len(v),covariance_rows=len(c),
    covariance_min_eigenvalue=float(np.linalg.eigvalsh(matrix).min()),
    caveat='Small n=3; delta method exploratory; biomass held owing to cross-sheet discrepancies; no independent second extraction.')
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
print(pd.DataFrame(mods).to_string(index=False))
