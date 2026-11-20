"""Extract original printed data table, not graph digitization or meta table.
Requires pymupdf,numpy. Source PDF hash and table/cell counts fixed.
Usage: python extract_anic2022.py CACHE_DIRECTORY
"""
import csv,hashlib,json,re,sys
from pathlib import Path
import numpy as np
import pymupdf
OUT=Path(__file__).resolve().parent
HASH='baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0'
def write(name,rows):
 with (OUT/name).open('w',encoding='utf8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 p=Path(sys.argv[1])/'anic2022_data.pdf';raw=p.read_bytes()
 assert hashlib.sha256(raw).hexdigest()==HASH,'original source changed'
 doc=pymupdf.open(stream=raw,filetype='pdf');text='\n'.join(page.get_text() for page in doc)
 section=text.split('Values in bold: not included in the shoot/root ratio analyses')[1].split('Nitrogen or Carbon content')[0]
 parsed=re.findall(r'([A-H][1-7])\s*\n(L[123])\s*\n(T[123])\s*\n([\d.]+)',section)
 assert len(parsed)==54 and len({r[0] for r in parsed})==54
 records=[dict(mesocosm=m,light=l,trophic=t,total_biomass_g=float(y)) for m,l,t,y in parsed]
 write('anic2022_community_units.csv',records)
 effects=[]
 for trophic in ['T1','T2','T3']:
  dark=np.array([r['total_biomass_g'] for r in records if r['trophic']==trophic and r['light']=='L1'])
  assert len(dark)==6
  for light in ['L2','L3']:
   lit=np.array([r['total_biomass_g'] for r in records if r['trophic']==trophic and r['light']==light]);assert len(lit)==6
   effects.append(dict(effect_id=f'Anic2022_community_{trophic}_{light}',study_id='Anic2022',programme_id='Penryn2012_2016',experiment_id='Penryn2012_2016',control_group_id=f'Penryn_{trophic}_dark',endpoint='community_total_dry_biomass_g',trophic=trophic,light=light,independent_unit='mesocosm',timepoint='2016_final_harvest',mean_lit=float(lit.mean()),sd_lit=float(lit.std(ddof=1)),n_lit=6,mean_dark=float(dark.mean()),sd_dark=float(dark.std(ddof=1)),n_dark=6,lnRR=float(np.log(lit.mean()/dark.mean())),sampling_variance=float(lit.var(ddof=1)/6/lit.mean()**2+dark.var(ddof=1)/6/dark.mean()**2),source_location='Figshare20304039 original PDF printed table Plant biomass per mesocosm',digitized=0,independent_verification='pending'))
 write('anic2022_effects.csv',effects)
 cov=[]
 for a in effects:
  for b in effects:
   v=a['sampling_variance'] if a['effect_id']==b['effect_id'] else (a['sd_dark']**2/6/a['mean_dark']**2 if a['control_group_id']==b['control_group_id'] else 0)
   cov.append(dict(effect_id_1=a['effect_id'],effect_id_2=b['effect_id'],sampling_covariance=v,method='shared_control_delta'))
 write('anic2022_covariance.csv',cov)
 write('anic2022_source_manifest.csv',[dict(source_url='https://ndownloader.figshare.com/files/36265164',source_doi='10.6084/m9.figshare.20304039',paper_doi='10.1002/ece3.9157',sha256=HASH,source_location='original PDF table; treatment codes in original supplementary Table S2',units=54,independent_programmes=1)])
 print(json.dumps(effects,indent=2))
if __name__=='__main__':main()
