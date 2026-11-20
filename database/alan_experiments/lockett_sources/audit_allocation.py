"""Original-data allocation audit, preserving conflicting labels; no effect fit."""
from pathlib import Path
import json
import hashlib
import pandas as pd

OUT=Path(__file__).resolve().parent
REPO=OUT.parents[2]
CACHE=REPO.parent/'gbif_annual/alan_experiment_sources/lockett2022'
metadata=json.loads((CACHE/'files.json').read_text(encoding='utf-8'))
manifest=[]
for f in metadata:
    path=CACHE/f['filename']
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    assert sha==f['content_details']['sha256_hash']
    assert path.stat().st_size==f['size']
    manifest.append(dict(filename=path.name,sha256=sha,bytes=path.stat().st_size,
                         url=f['content_details']['download_url']))
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.csv',index=False)
book=next(CACHE.glob('*Data (trees).xlsx'))
biomass=pd.read_excel(book,'Biomass and SLA')
growth=pd.read_excel(book,'Growth')
final=growth[growth.week.eq(38)]
pca=pd.read_csv(next(CACHE.glob('*Data (PCA).csv')))
main=pd.read_csv(next(CACHE.glob('*Data (main).csv')))
assert len(biomass)==23 and biomass.tree.is_unique and biomass['plot'].is_unique
assert len(final)==144 and final.tree.is_unique and pca.tree.is_unique
x=biomass[['tree','plot','light_treatment','ALAN_treatment','total_dryweight']].merge(
    final[['tree','plot1','plot38','light_treatment1','light_treatment38','swaptype']],
    on='tree',validate='one_to_one').merge(pca[['tree','treatment']],on='tree',validate='one_to_one')
assert len(x)==23 and x['plot'].eq(x.plot38).all()
x['biomass_growth_label_conflict']=x.light_treatment.str.lower().ne(x.light_treatment38.str.lower())
x['pca_light_label']=x.treatment.map({'ALAN':'light','No ALAN':'dark'})
assert x.pca_light_label.notna().all()
x['biomass_pca_label_conflict']=x.light_treatment.str.lower().ne(x.pca_light_label)
x['growth_pca_label_conflict']=x.light_treatment38.str.lower().ne(x.pca_light_label)
x['shared_origin_among_harvested']=x.plot1.duplicated(keep=False)
x.to_csv(OUT/'harvest_allocation_crosswalk.csv',index=False)
x[x.biomass_growth_label_conflict].to_csv(OUT/'treatment_label_discrepancies.csv',index=False)
x[x.shared_origin_among_harvested][['tree','plot1','plot38','light_treatment1','light_treatment38']].to_csv(OUT/'shared_origin_trees.csv',index=False)
assert x.biomass_growth_label_conflict.sum()==4
assert set(x.loc[x.biomass_growth_label_conflict,'tree'])=={7,14,60,91}
assert x.plot1.nunique()==18 and x.shared_origin_among_harvested.sum()==10
assert not x.growth_pca_label_conflict.any()
# Metadata preservation only: no NA -> absence, no guessing a treatment.
columns=['ID','tree','bag','plot1','plot2','ALANprior','ALANtrmt','Use4ALAN','FAILED','MISTIMED']
main[columns].to_csv(OUT/'psyllid_allocation_and_flags.csv',index=False)
main.groupby(['FAILED','Use4ALAN'],dropna=False).size().rename('bags').reset_index().to_csv(OUT/'failure_by_analysis_flag.csv',index=False)
main[main.tree.isin(x.loc[x.biomass_growth_label_conflict,'tree'])][columns].to_csv(OUT/'disputed_tree_psyllid_rows.csv',index=False)
audit=dict(doi='10.1016/j.envpol.2022.119803',data_doi='10.17632/pfrrkgs878.1',
    source_files_verified=len(manifest),harvest_trees=len(x),final_harvest_plots=x.plot38.nunique(),
    initial_harvest_plots=x.plot1.nunique(),shared_initial_plot_groups=x.loc[x.shared_origin_among_harvested,'plot1'].nunique(),
    conflicting_biomass_treatment_rows=int(x.biomass_growth_label_conflict.sum()),
    growth_and_pca_labels_agree=bool(~x.growth_pca_label_conflict.any()),
    psyllid_rows=len(main),psyllid_trees=main.tree.nunique(),
    psyllid_current_light_missing=int(main.ALANtrmt.isna().sum()),
    psyllid_use4alan_yes=int(main.Use4ALAN.eq('Yes').sum()),
    main_catalogue_changed=False,effects_extracted=0,
    decision='hold_treatment_labels_and_crossed_exposure_units',
    constraints=['Labels preserved as reported, no majority-vote correction',
        'One tree per harvest plot does not remove shared initial-plot exposure',
        'Psyllid bags nested in trees and crossed plots are not independent light replicates',
        'Missing current treatment or omitted reproductive trials are not zero responses',
        'Lerp production is not pollination, predation, or plant reproductive output'])
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps(audit,indent=2))
