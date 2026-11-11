#!/usr/bin/env python3
"""GBIF gecko positive-detection pilot. LIVE API; NO inferred gecko absences.

Output intentionally contains ONLY common-gecko 1-degree year cells and
aggregate rare day-gecko counts. No rare gecko point coordinates, no plant
photo site identities, no native/invasive assignments are published.
"""
import argparse,csv,json,math,time
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
import requests

BASE='https://api.gbif.org/v1'
GROUPS={
 'Hemidactylus frenatus':'house_nocturnal',
 'Hemidactylus mabouia':'house_nocturnal',
 'Hemidactylus turcicus':'house_nocturnal',
 'Lepidodactylus lugubris':'house_nocturnal',
 'Phelsuma ornata':'day_flower_visitor',
 'Phelsuma cepediana':'day_flower_visitor',
 'Gehyra dubia':'native_gecko_reference_Australia'}
OBS={'HUMAN_OBSERVATION','MACHINE_OBSERVATION','OBSERVATION'}

def need(x,msg):
 if not x:raise ValueError(msg)

def api(s,path,params):
 for retry in range(4):
  try:
   r=s.get(BASE+path,params=params,timeout=55)
   if r.status_code in (429,500,502,503,504):raise requests.RequestException(f'HTTP {r.status_code}')
   r.raise_for_status();return r.json()
  except (requests.RequestException,ValueError) as e:
   if retry==3:raise RuntimeError(f'GBIF API request failed path={path}, query={params}') from e
   time.sleep(2**retry)
 raise AssertionError()

def taxon(s,name):
 r=api(s,'/species/match',{'name':name,'kingdom':'Animalia'})
 need(r.get('canonicalName')==name and r.get('rank')=='SPECIES','Ambiguous GBIF species match: '+name)
 need(r.get('status') in ('ACCEPTED','SYNONYM'),'Unresolved GBIF species')
 k=r.get('acceptedUsageKey') if r.get('status')=='SYNONYM' else r.get('usageKey')
 need(isinstance(k,int) and k>0,'Missing GBIF taxon usage key')
 return k,r

def record(r,name,year):
 status=str(r.get('occurrenceStatus','DETECTED')).upper()
 if status in ('ABSENT','NOT_DETECTED'):return None,'not_present'
 need(status in ('PRESENT','DETECTED'),f'Unsupported GBIF occurrence status: {status}')
 if r.get('basisOfRecord') not in OBS:return None,'not_direct_observation'
 if r.get('year')!=year:return None,'not_exact_year'
 if not str(r.get('species') or r.get('acceptedScientificName') or '').startswith(name):
  return None,'taxon_not_exact'
 try:lat,lon=float(r['decimalLatitude']),float(r['decimalLongitude'])
 except (TypeError,ValueError,KeyError):return None,'not_georeferenced'
 if not (math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180):
  return None,'invalid_coordinates'
 if r.get('hasGeospatialIssue') is True:return None,'geospatial_issue'
 try:unc=float(r['coordinateUncertaintyInMeters'])
 except (TypeError,ValueError,KeyError):return None,'precision_unknown'
 if not (math.isfinite(unc) and 0<=unc<=25000):return None,'uncertainty_gt25km'
 return (math.floor(lat),math.floor(lon),str(r.get('datasetKey') or 'unknown')),None

def year_source(s,k,year,cap=50000):
 q={'taxonKey':k,'year':str(year),'hasCoordinate':'true','hasGeospatialIssue':'false','limit':300,'offset':0}
 first=api(s,'/occurrence/search',q);count=first.get('count')
 need(isinstance(count,int),'GBIF total missing')
 need(count<=min(cap,100000),f'GBIF source requires per-year split, never truncate: {count}')
 rows=list(first.get('results',[]));i=len(rows)
 while i<count:
  q['offset']=i
  page=api(s,'/occurrence/search',q).get('results',[])
  need(0<len(page)<=300,'GBIF pagination incomplete')
  rows.extend(page);i+=len(page)
 need(len(rows)==count,'GBIF source changed count during pagination')
 keys=[r.get('key') for r in rows]
 need(all(keys) and len(set(keys))==len(keys),'Duplicate GBIF occurrence keys')
 return rows

def run(years,names,out,*,cap=50000,session=None):
 out=Path(out)
 need(not out.exists(),'Refuse source overwrite')
 need(len(set(years))==len(years) and years and all(2012<=y<=2025 for y in years),'Invalid years')
 need(names and len(set(names))==len(names) and all(n in GROUPS for n in names),'Unknown taxon')
 s=session or requests.Session()
 s.headers.update({'User-Agent':'ALAN-GloBI gecko GBIF positive-only source audit/1.0'})
 diagnostics=[];grid=defaultdict(set)
 for name in names:
  key,tax=taxon(s,name);audit=[]
  for year in years:
   rows=year_source(s,key,year,cap=cap)
   valid=0;excluded=Counter();grids=set();datasets=set()
   for r in rows:
    obs,why=record(r,name,year)
    if why:excluded[why]+=1;continue
    valid+=1;lat,lon,dataset=obs;datasets.add(dataset);grids.add((lat,lon))
    if GROUPS[name]=='house_nocturnal':grid[(year,lat,lon)].add(name)
   need(valid+sum(excluded.values())==len(rows),'GBIF count reconciliation failed')
   audit.append({'year':year,'source_rows':len(rows),'qualified_direct_records':valid,
                 'occupied_1degree_source_grids':len(grids),
                 'source_datasets':len(datasets),'exclusions':dict(excluded)})
  diagnostics.append({'species':name,'functional_group_prior':GROUPS[name],
    'accepted_match_taxonKey':key,'GBIF_match_type':tax.get('matchType'),
    'GBIF_match_confidence':tax.get('confidence'),'years':audit})
 out.mkdir(parents=True)
 # For rare Phelsuma we deliberately release only aggregate species/year counts.
 with (out/'house_gecko_positive_1degree_cells.csv').open('w',newline='') as f:
  w=csv.writer(f)
  w.writerow(['year','grid_lat_floor','grid_lon_floor','documented_common_gecko_names'])
  for (year,lat,lon),names in sorted(grid.items()):
   w.writerow([year,lat,lon,len(names)])
 receipt={'schema':'alan_GBIF_gecko_positive_source_pilot_v1',
  'source':'GBIF live public occurrence/search, exploratory not DOI-frozen',
  'retrieved_utc':datetime.now(timezone.utc).isoformat(),
  'years':years,'species_audit':diagnostics,
  'documented_house_gecko_gridyear_counts':{str(y):sum(x[0]==y for x in grid) for y in years},
  'house_gecko_positive_grid_is_not_gecko_absence':True,
  'locally_native_or_invasive_gecko_status_not_estimated':True,
  'flower_pollination_predation_not_estimated':True,
  'rare_flower_gecko_locations_suppressed':True,
  'GBIF_observation_effort_not_normalized':True}
 (out/'source_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 return receipt

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--years',nargs='+',type=int,default=[2019,2024])
 p.add_argument('--taxa',nargs='+',default=list(GROUPS))
 p.add_argument('--cap',type=int,default=50000)
 p.add_argument('--out',type=Path,required=True)
 a=p.parse_args()
 result=run(a.years,a.taxa,a.out,cap=a.cap)
 print(json.dumps({'counts':result['documented_house_gecko_gridyear_counts'],
   'audited_species_years':sum(len(x['years']) for x in result['species_audit'])}))
