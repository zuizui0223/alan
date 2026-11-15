"""P1 geographic coverage and exploratory day-only/open functional contrast.
Authentic Kendall & Nicholson 2024-09-13 author source, SHA guarded.
Not a causal island or ALAN treatment estimate. No source coordinates exported.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.alan_functional_dependence_20261115 import SOURCE_SHA,match_panels,endpoint_subtype

# Source-validated exact study / country / location. Never classify all Japan,
# Great Britain or Malaysia as either continental mainland or small oceanic island.
LOOKUP={
'Barrios_2012':('Cuba','Pan de Matanzas','caribbean','Barrios'),
'Barrios_2022':('Cuba','Pan de Matanzas','caribbean','Barrios'),
'Jaca_2020':('Canary Islands','Tenerife','oceanic','Jaca'),
'McMullen_2007':('Ecuador','Pinta Island','oceanic','McMullen'),
'McMullen_2009':('Ecuador','Galapagos','oceanic','McMullen'),
'McMullen_2011':('Ecuador','Santa Cruz Island','oceanic','McMullen'),
'McMullen_2012':('Ecuador','Galapagos (Santa Cruz)','oceanic','McMullen'),
'Nakajima_2018':('Japan','Izu islands','oceanic','Nakajima'),
'Chapurlat_2015':('Sweden','Skane','shelf_small','Chapurlat'),
'Chapurlat_2019':('Sweden','Skane','shelf_small','Chapurlat'),
'Noone_2022':('Canada','Kent Island','shelf_small','Noone'),
'Suetsugu_2018':('Japan','Yakushima Island','shelf_small','Suetsugu'),
'Yumoto_2000':('Malaysia','Lambir Hills National Park','large_island','Yumoto'),
'Macgregor_2019':('United Kingdom','Yorkshire','large_island','Macgregor'),
'Funamoto_2017':('Japan','Tsukuba','large_island','Funamoto'),
'Funamoto_2020':('Japan','Ashiu','large_island','Funamoto'),
'Hattori_2020':('Japan','Kaida','large_island','Hattori')}
CONTINENTS={'Argentina','Australia','Belgium','Bolivia','Brazil','Canada',
'Chile','China','Colombia','Costa Rica','Ecuador','Germany','India','Italy',
'Kenya','Mexico','Norway','Panama','South Africa','Spain','Sweden','Switzerland','United States'}
SMALL={'oceanic','caribbean','shelf_small'}

def classify(p):
    p=p.copy()
    typ=[];family=[];geo_error=[]
    for r in p.itertuples(index=False):
        if r.study_ID in LOOKUP:
            country,place,group,author=LOOKUP[r.study_ID]
            if (r.country,r.locality)!=(country,place):
                raise ValueError('Study locality drift '+r.study_ID)
        else:
            if r.country not in CONTINENTS:
                raise ValueError('Unverified island/mainland study '+r.study_ID+' '+r.country)
            group='continent';author=r.study_ID.split('_')[0]
        typ.append(group);family.append(author)
        is_pinta=r.study_ID=='McMullen_2007'
        if is_pinta and not float(r.coordinates_lat)<0:
            raise ValueError('Pinta coordinate sign has changed; re-audit')
        geo_error.append(is_pinta)
    p['geo_class']=typ;p['source_family']=family;p['geo_warning']=geo_error
    p['subtype']=[endpoint_subtype(m,u) for m,u in zip(p.treatment_effectiveness_metric,p.effectiveness_units)]
    p['log_day_open']=np.nan
    good=p.day_mean.gt(0)
    p.loc[good,'log_day_open']=np.log(p.loc[good,'day_mean']/p.loc[good,'open_mean'])
    return p

def compare(p,subtype,classes=('oceanic',),seed=20261115):
    a=p[(p.subtype==subtype)&p.geo_class.isin(classes)&(~p.geo_warning)]
    b=p[(p.subtype==subtype)&p.geo_class.eq('continent')]
    x=a.groupby('study_ID').night_loss.median().sort_index()
    y=b.groupby('study_ID').night_loss.median().sort_index()
    out={'subtype':subtype,'island_groups':list(classes),
         'island_panels':len(a),'island_studies':len(x),'island_author_families':a.source_family.nunique(),
         'mainland_panels':len(b),'mainland_studies':len(y),
         'excluded_Pinta_geocoding_warning':True,
         'island_study_effects':[{'study':k,'night_loss':float(v)} for k,v in x.items()],
         'same_study_paired':False,'confirmatory':False}
    if len(x) and len(y):
        out['island_study_median']=float(x.median())
        out['mainland_study_median']=float(y.median())
        out['exploratory_island_minus_mainland_night_loss']=float(x.median()-y.median())
        if len(x)>=2 and len(y)>=2:
            rng=np.random.default_rng(seed);B=12000
            xs=rng.choice(x.to_numpy(float),(B,len(x)),replace=True)
            ys=rng.choice(y.to_numpy(float),(B,len(y)),replace=True)
            out['descriptive_paper_bootstrap95']=(np.quantile(np.median(xs,axis=1)-np.median(ys,axis=1),[.025,.975])).tolist()
            out['bootstrap_ignores_shared_authors_and_phylogeny']=True
        q=p[(p.subtype==subtype)&p.geo_class.isin(list(classes)+['continent'])&(~p.geo_warning)]
        v=q[q.geo_class.isin(classes)].groupby('study_ID').log_day_open.median().dropna()
        w=q[q.geo_class.eq('continent')].groupby('study_ID').log_day_open.median().dropna()
        out['protocol_log_D_over_O_comparison']={
            'island_studies':len(v),'mainland_studies':len(w),
            'island_minus_mainland_log_median':float(v.median()-w.median()) if len(v) and len(w) else None,
            'zero_day_means_excluded_without_pseudocount':int(q.day_mean.eq(0).sum())}
    return out

def analyze(raw):
    matched,audit=match_panels(raw)
    p=matched[matched.treatment_effectiveness_metric.isin(['fruit set','seed set'])].copy()
    p=classify(p)
    if len(raw)!=1652 or len(p)!=340 or p.study_ID.nunique()!=113:
        raise ValueError('Source design count changed')
    count={}
    for cls,d in p.groupby('geo_class'):
        count[cls]={'n_panels':len(d),'n_studies':d.study_ID.nunique(),
                    'n_species':d.plant_species.nunique(),
                    'n_author_families':d.source_family.nunique(),
                    'n_geocoding_warning':int(d.geo_warning.sum())}
    expected={'oceanic':6,'caribbean':2,'shelf_small':4,'large_island':5,'continent':96}
    if {k:v['n_studies'] for k,v in count.items()}!=expected:
        raise ValueError('Source geography classification drift')
    overlap=p.groupby('study_ID').geo_class.nunique()
    if (overlap>1).any():raise ValueError('Within-paper island/mainland pairing changed')
    a=p[(p.geo_class.isin(SMALL))&(~p.geo_warning)]
    b=p[p.geo_class.eq('continent')]
    shared=sorted(set(a.plant_species)&set(b.plant_species))
    exact=[]
    for sp in shared:
        aa=a[a.plant_species.eq(sp)];bb=b[b.plant_species.eq(sp)]
        for ep in sorted(set(aa.subtype)&set(bb.subtype)):
            x=aa[aa.subtype.eq(ep)].groupby('study_ID').night_loss.median()
            y=bb[bb.subtype.eq(ep)].groupby('study_ID').night_loss.median()
            exact.append({'species':sp,'subtype':ep,'island_studies':list(x.index),
                          'mainland_studies':list(y.index),
                          'difference_night_loss':float(x.median()-y.median()),
                          'same_study_or_same_year_comparison':False})
    L=p[p.study_ID.eq('Nakajima_2018')]
    return {'schema':'ALAN_P1_island_mainland_diel_exclusion_source_20261115_v1',
      'source_doi':'10.1111/ele.70036','source_sha256':SOURCE_SHA,
      'source_rows':len(raw),'source_reproductive_panels':len(p),'source_studies':p.study_ID.nunique(),
      'geographic_groups':count,'small_sea_island_n_panels':int(p.geo_class.isin(SMALL).sum()),
      'small_sea_island_n_studies':p[p.geo_class.isin(SMALL)].study_ID.nunique(),
      'source_locality_mislabeled_as_Skane_but_primary_Oland':['Chapurlat_2015','Chapurlat_2019'],
      'within_study_island_mainland_overlap':int((overlap>1).sum()),
      'oceanic_versus_continent':{ep:compare(p,ep) for ep in ('fruit_proportion','seed_per_fruit','seed_per_ovule')},
      'all_small_islands_sensitivity':{ep:compare(p,ep,tuple(sorted(SMALL))) for ep in ('fruit_proportion','seed_per_fruit')},
      'same_species_cross_study_diagnostic':exact,
      'pinta_geocoding_warning':'Pinta Island archived latitude -0.542222 vs official +0.586 N, not corrected or used in geographic contrast',
      'supported_prior_nocturnal_shift_source':{'study_ID':'Nakajima_2018',
        'source_doi':'10.1111/1442-1984.12222',
        'island_night_loss_fruit_proportion':float(L[L.subtype.eq('fruit_proportion')].night_loss.iloc[0]),
        'island_night_loss_seed_per_ovule':float(L[L.subtype.eq('seed_per_ovule')].night_loss.iloc[0]),
        'no_corresponding_mainland_same_study_exclusion_row':True},
      'directional_hypothesis_H4_confirmed':False,
      'P1_island_shift_x_ALAN_effect_estimable':False,
      'reason':'no paired island/mainland study arms, endpoint island clusters <=2 research families, prior channel shift not documented for all island studies',
      'no_ALAN_manipulation_in_this_panel':True,
      'missing_not_zero':True,'not_an_inverse_variance_or_causal_meta_analysis':True,
      'study_family_dependence_not_repaired_by_paper_bootstrap':True}

def main():
    z=argparse.ArgumentParser()
    z.add_argument('--source',type=Path,required=True)
    z.add_argument('--out',type=Path,required=True)
    a=z.parse_args()
    raw=a.source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError('Authentic source SHA mismatch')
    result=analyze(pd.read_csv(a.source))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('ALAN_P1_ISLAND_SOURCE_QA',json.dumps({'n':result['source_reproductive_panels'],
      'papers':result['source_studies'],'classes':result['geographic_groups'],
      'oceanic_seed_fruit':result['oceanic_versus_continent']['seed_per_fruit']},ensure_ascii=False))

if __name__=='__main__':main()
