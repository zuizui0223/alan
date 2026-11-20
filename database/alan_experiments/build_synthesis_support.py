"""Count biological estimands and programmes, not just finite effect rows.

This is a candidate-support manifest, not an eligibility override or a pooled
analysis. Supplementary extractions outside the main catalogue stay outside.
"""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent
d=pd.read_csv(ROOT/'effect_catalogue.csv').fillna('')
assert d.effect_id.is_unique
def classify(r):
    e=r.endpoint
    if e in ['aquatic_origin_capture','terrestrial_origin_capture','all_invertebrate_capture','cicadellidae_capture']:return e
    if e=='spider_abundance':return 'spider_abundance'
    if e=='web_prey_count':return 'web_prey_count_not_predation_rate'
    if e in ['aphid_count_3min','aphid_abundance']:return 'aphid_abundance'
    if e=='parasitoid_mummy_abundance':return 'standing_mummy_count_not_parasitism_rate'
    if e.endswith('_herbivory'):return 'plant_herbivory'
    if r.study_id in ['Giavi2021','Lu2026']:
        if e in ['three_orders_DHC','daytime_selected21_visits_all_three_orders']:return 'day_visits_three_order_sum'
        for full,short in [('Diptera','D'),('Hymenoptera','H'),('Coleoptera','C')]:
            if e.endswith('_'+full) or e==short+'_abundance':return 'day_visits_'+full
        if e=='abundance':return 'day_visits_source_total_unresolved'
        return 'day_visits_additional_taxa'
    if e in ['community_total_dry_biomass_g','grassland_biomass_live_total']:return 'community_living_dry_biomass'
    if e=='grassland_biomass_total':return 'community_dry_biomass_including_dead'
    if e=='grassland_biomass_Dead':return 'dead_plant_material'
    if 'biomass' in e:return 'species_dry_biomass'
    raise ValueError(f'Unclassified endpoint: {r.study_id} {e}')
d['candidate_group']=d.apply(classify,axis=1)
d['supplementary_time_summary']=d.analysis_role.eq('supplementary')
d['source_total_unresolved']=d.candidate_group.eq('day_visits_source_total_unresolved')
d['candidate_primary_row']=~(d.supplementary_time_summary|d.source_total_unresolved)
# Sum/parts and alternative illumination levels remain dependent even here.
d['candidate_primary_row_is_approval']=False
moist=pd.read_csv(ROOT/'site_context/plant_moisture_effect_support.csv')
moist=moist[moist.finite_effect_with_moisture==True][['effect_id','source_species','EIV_MOIST']]
assert moist.effect_id.is_unique
d=d.merge(moist,on='effect_id',how='left',validate='one_to_one')
d[['effect_id','study_id','programme_id','endpoint','candidate_group','analysis_role',
   'supplementary_time_summary','source_total_unresolved','candidate_primary_row',
   'candidate_primary_row_is_approval','approval_status','source_table','covariance_table',
   'source_species','EIV_MOIST']].to_csv(ROOT/'synthesis_support_effects.csv',index=False)
rows=[]
for group,g in d.groupby('candidate_group',sort=True):
    p=g[g.candidate_primary_row]
    rows.append(dict(candidate_group=group,catalogue_rows=len(g),candidate_rows=len(p),
                     paper_labels=p.study_id.nunique(),programme_labels=p.programme_id.nunique(),
                     programmes=';'.join(sorted(p.programme_id.unique())),
                     excluded_supplementary_time_rows=int(g.supplementary_time_summary.sum()),
                     source_unresolved_rows=int(g.source_total_unresolved.sum()),
                     rows_with_species_moisture=int(p.EIV_MOIST.notna().sum()),
                     taxa_with_species_moisture=p.loc[p.EIV_MOIST.notna(),'source_species'].nunique(),
                     programmes_with_species_moisture=p.loc[p.EIV_MOIST.notna(),'programme_id'].nunique(),
                     adjudication='pending; candidate does not mean accepted',pooled_estimate_available=False))
summary=pd.DataFrame(rows)
assert summary.catalogue_rows.sum()==len(d)==108
assert summary.candidate_rows.sum()==105
summary.to_csv(ROOT/'synthesis_support_groups.csv',index=False)
assert int(d.EIV_MOIST.notna().sum())==10
audit=dict(catalogue_rows=len(d),candidate_rows_before_final_adjudication=105,
           excluded_duplicate_time_summaries=2,source_total_inconsistency_rows=1,
           candidate_groups=len(summary),programme_labels=d.programme_id.nunique(),
           finite_effects_with_species_moisture=10,moisture_taxa=5,moisture_programmes=3,
           fully_adjudicated_effects=0,pooled_estimates=0,
           exclusions_are_not_literature_exclusions=True,
           warning='105 rows are not105 independent studies; sum/parts and doses overlap. Programme labels still require adjudication. Species moisture is not site wetness.')
(ROOT/'synthesis_support_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(summary[['candidate_group','catalogue_rows','candidate_rows','programme_labels']].to_string(index=False))
