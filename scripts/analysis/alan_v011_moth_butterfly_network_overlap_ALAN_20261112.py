#!/usr/bin/env python3
"""Same 10,698 original ALAN photographic sites; GloBI moth/butterfly positive network redundancy.

This code reads PRIVATE source-generated photo×accepted-species rows but
publishes aggregate-only results. GloBI not-recorded is UNKNOWN biological use.
Both documented guilds is NOT proof of daytime functional compensation.
"""
import argparse,gzip,hashlib,io,json,math,zipfile
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
from alan_globi_crossland_light_robust_20261111 import cleaned,prepare,fit

GLOBI_SOURCE_SHA='670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3'
FLAGS_SHA='cf750192b4d641a2d904f4f74d5b87d41fc97e04822a67d513d9bad74ece8e01'

def require(ok,msg):
    if not ok:raise ValueError(msg)

def run(site_file,site_species_file,full_flags_zip):
    sites=cleaned(site_file)
    sites=prepare(sites,3,2)
    sites=sites.replace([np.inf,-np.inf],np.nan).dropna(
        subset=['bio5','lograin','elevation_km','logpop','loglight','logphotos','logspecies'])
    sites=sites.loc[sites.groupby('region_year').land_class.transform('nunique').eq(2)].copy()
    require(len(sites)==10698 and int(sites.island.sum())==933 and
      sites.region5.nunique()==78 and sites.region_year.nunique()==195,
      'Frozen 2019/21/23/24/25 source photo geography changed')
    old=fit(sites,'all','moth',balance=True)
    require(abs(old['beta_island_minus_mainland_slope_pp_per_log1p']+1.901720622055596)<1e-9,
      'Original observed positive source result is not reproduced')
    with zipfile.ZipFile(full_flags_zip) as zip_file:
        member='GloBI_v011_positive_species_flags.csv.gz'
        source=zip_file.read(member)
        require(hashlib.sha256(source).hexdigest()==FLAGS_SHA,'GloBI full positive-flag SHA changed')
        receipt=json.loads(zip_file.read('source_receipt.json'))
        require(receipt['sha256']==GLOBI_SOURCE_SHA and receipt['all_interaction_rows']==24453954,
                'GloBI is not the SHA-verified complete v0.11 source')
        flags=pd.read_csv(io.BytesIO(gzip.decompress(source)))
    sp=pd.read_csv(site_species_file)
    require({'target_id','wcvp_accepted_name'}<=set(sp),'Missing canonical observed plant identity')
    require(not sp.duplicated(['target_id','wcvp_accepted_name']).any(),'Duplicate recorded accepted plant')
    sp=sp.loc[sp.target_id.isin(sites.target_id)].copy()
    sp=sp.merge(flags,left_on='wcvp_accepted_name',right_on='plant_name',
                how='left',validate='many_to_one')
    for col in flags.columns:
        if col.endswith('_any_family'):
            sp[col]=sp[col].fillna(0).astype(int)
    require(len(sp)==int(sites.n_species.sum()),'Actual photographed accepted denominator changed')
    labels={
      'both':sp.all_moth_any_family.eq(1)&sp.all_butterfly_any_family.eq(1),
      'moth_recorded_no_butterfly_record':sp.all_moth_any_family.eq(1)&sp.all_butterfly_any_family.eq(0),
      'butterfly_recorded_no_moth_record':sp.all_moth_any_family.eq(0)&sp.all_butterfly_any_family.eq(1),
      'moth_any':sp.all_moth_any_family.eq(1),
      'butterfly_any':sp.all_butterfly_any_family.eq(1),
      'non_iNat_trophiCH_both':sp.non_inat_trophich_moth_any_family.eq(1)&sp.non_inat_trophich_butterfly_any_family.eq(1),
      'non_iNat_trophiCH_moth_recorded_no_butterfly_record':sp.non_inat_trophich_moth_any_family.eq(1)&sp.non_inat_trophich_butterfly_any_family.eq(0),
    }
    for k,v in labels.items():sp[k]=v.astype(int)
    count=sp.groupby('target_id')[list(labels)].sum()
    frame=sites.join(count,on='target_id',validate='one_to_one')
    require(frame[list(labels)].notna().all().all(),'Unmatched accepted species source')
    result={}
    for col in labels:
        frame['all_moth']=frame[col]
        fit_result=fit(frame,'all','moth',balance=True)
        positive=frame.loc[frame[col].gt(0)]
        result[col]={
          'documented_photo_plant_names':int(sp.loc[sp[col].eq(1),'wcvp_accepted_name'].nunique()),
          'positive_photo_island_site_years':int(positive.island.sum()),
          'positive_photo_mainland_site_years':int((1-positive.island).sum()),
          'beta':fit_result['beta_island_minus_mainland_slope_pp_per_log1p'],
          'CI95':fit_result['95CI_cluster_5deg'],
          'p_nominal':fit_result['p_nominal']}
    contrasts={
      'both_minus_moth_recorded_no_butterfly_record':('both','moth_recorded_no_butterfly_record'),
      'both_minus_butterfly_recorded_no_moth_record':('both','butterfly_recorded_no_moth_record'),
      'all_moth_minus_all_butterfly':('moth_any','butterfly_any'),
      'strict_both_minus_strict_moth_recorded_no_butterfly_record':
            ('non_iNat_trophiCH_both','non_iNat_trophiCH_moth_recorded_no_butterfly_record')
    }
    direct={}
    for k,(a,b) in contrasts.items():
        frame['all_moth']=frame[a]-frame[b]
        r=fit(frame,'all','moth',balance=True)
        direct[k]={'beta':r['beta_island_minus_mainland_slope_pp_per_log1p'],
                   'CI95':r['95CI_cluster_5deg'],'p_nominal':r['p_nominal']}
    qs=false_discovery_control([v['p_nominal'] for v in direct.values()],method='bh')
    for v,q in zip(direct.values(),qs):v['BH_q_four_postexploratory_tests']=float(q)
    return {
      'schema':'ALAN_GloBI_full_v011_multichannel_2019_2025_10698_photo_siteyears_v1',
      'GloBI_full_source_sha256':GLOBI_SOURCE_SHA,
      'original_legacy_moth_spatial_positive_beta_frozen':-1.901720622055596,
      'n_original_photo_site_years':len(sites),'original_photographed_plants':int(sp.wcvp_accepted_name.nunique()),
      'n_island_photo_site_years':int(sites.island.sum()),
      'n_independent_5degree_geo_clusters':int(sites.region5.nunique()),
      'guild_positive_groups':result,'direct_group_contrasts':direct,
      'no_GloBI_source_claim_is_NOT_true_guild_absence':True,
      'both_GloBI_guild_links_not_demonstrated_effective_pollination_compensation':True,
      'photo_observation_fraction_not_true_plant_abundance':True,
      'frozen_README_PL64_nominal_positive_still_reported_separately':True,
      'no_private_original_photo_or_species_site_rows_exported':True}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--private-sites',type=Path,required=True)
    p.add_argument('--private-site-species',type=Path,required=True)
    p.add_argument('--globi-flags-zip',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise ValueError('Refuse overwrite')
    r=run(a.private_sites,a.private_site_species,a.globi_flags_zip)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('GLOBI_MULTI_CHANNEL_REAL_SOURCE',json.dumps(r['direct_group_contrasts']))
