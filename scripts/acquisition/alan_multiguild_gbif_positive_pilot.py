#!/usr/bin/env python3
"""Pilot independent GBIF POSITIVE-only species records for GloBI predator guilds.

No geo-privacy-sensitive source points, no species abundance estimates, no
absence from no GBIF observation, no assumption of local interaction.
Source failures and provider caps remain UNRESOLVED, not zero.
"""
from __future__ import annotations
import argparse,json,re
from pathlib import Path
from scripts.acquisition import alan_gecko_gbif_positive_pilot as gbif

TARGETS={
 'Glossophaga soricina':'bat_dual_role_database_name',
 'Antrozous pallidus':'bat_dual_role_database_name',
 'Mystacina tuberculata':'bat_dual_role_database_name',
 'Thereuopoda clunifera':'centipede_dual_role_database_name',
 'Scolopendra subspinipes':'centipede_flower_visit_database_name',
 'Vespa velutina':'wasp_dual_role_database_name',
 'Argiope bruennichi':'spider_dual_role_database_name',
 'Mantis religiosa':'mantis_dual_role_database_name',
}

def require(ok,msg):
 if not ok:raise ValueError(msg)

def run(root,years=(2019,2024),cap=12000,selected=None):
 root=Path(root)
 require(not root.exists(),'Refuse overwrite of source run')
 require(all(y in (2019,2024) for y in years) and years,'Unexpected source years')
 choose=list(TARGETS) if selected is None else selected
 require(choose and set(choose)<=set(TARGETS),'Unknown source-positive taxon')
 require(cap>=500,'No silent unreasonable provider cap')
 gbif.GROUPS.update(TARGETS)
 root.mkdir(parents=True)
 rows=[]
 for name in choose:
  output=root/re.sub(r'[^A-Za-z0-9]+','_',name)
  try:
   receipt=gbif.run(years,[name],output,cap=cap)
   proof=receipt['species_audit'][0]
   counts=proof['years']
   require(len(counts)==len(years),'GBIF source-year count drift')
   source_status='verified_positive_observation_audit'
   error=None
  except (ValueError,RuntimeError,KeyError) as exc:
   source_status='UNRESOLVED_source_incomplete_or_taxonomy'
   error=str(exc)[:240]
   counts=[]
  rows.append({'scientific_name':name,'GloBI_role_prior':TARGETS[name],
   'source_status':source_status,
   'years':[{k:r[k] for k in ('year','source_rows','qualified_direct_records',
                             'occupied_1degree_source_grids','source_datasets','exclusions')}
             for r in counts],
   'error':error})
 report={
  'schema':'alan_multiguild_GBIF_positive_direct_observations_2019_2024_pilot_v1',
  'source':'GBIF live source API species/match and occurrence/search, observational HUMAN or MACHINE records with <=25km precision',
  'candidate_species_chosen_from_source_verified_full_GloBI_v011_role_names':True,
  'GloBI_v011_SHA256':'670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3',
  'years':list(years),'maximum_untruncated_source_records_per_species_year':cap,
  'species_audit':rows,
  'false_absence_from_no_record':False,
  'confirmed_predator_island_absence':False,
  'actual_predator_on_moth_interaction_observed_at_GBIF_site':False,
  'verified_pollination_service_or_bat_diet':False,
  'original_gecko_or_bat_coordinates_exported':False,
  'GBIF_datasets_licensed_source_references_must_be_preserved':True,
 }
 (root/'PUBLIC_source_support_receipt.json').write_text(
  json.dumps(report,indent=2,ensure_ascii=False)+'\n')
 return report

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--out',type=Path,required=True)
 p.add_argument('--years',nargs='+',type=int,default=[2019,2024])
 p.add_argument('--cap',type=int,default=12000)
 a=p.parse_args()
 r=run(a.out,years=a.years,cap=a.cap)
 print('MULTIGUILD_GBIF_POSITIVE_FEASIBILITY',json.dumps(
   [{'name':s['scientific_name'],'status':s['source_status'],
     'n':[x['qualified_direct_records'] for x in s['years']]}
    for s in r['species_audit']],sort_keys=True),flush=True)
