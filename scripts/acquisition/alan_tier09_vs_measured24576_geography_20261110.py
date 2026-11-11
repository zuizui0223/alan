#!/usr/bin/env python3
"""SHA-fixed tier09 light-support overlap with already measured tier08 24,576 originals.

Source-only comparison; no new RGB, no causal ALAN or moth-obligacy inference.
The attempted/ROI9 universe here is 24,576, NOT the full pooled 112,428
and NOT the 133-island multisite *regression* support.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

SCHEMA = 'ALAN_tier09_vs_measured24576_source_geography_v1'

def require(ok, why):
    if not bool(ok):
        raise ValueError('ALAN_TIER09_GEO: ' + why)

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):
            h.update(b)
    return h.hexdigest()

def score(old, new, pairs):
    require(old.photo_id.is_unique and new.photo_id.is_unique,
            'source photo ID duplicated')
    require(set(old.photo_id).isdisjoint(new.photo_id),
            'tier08 and tier09 original photo ID duplicated')
    require(old.strict9.isin((True,False)).all() and
            new.viirs_eligible.isin((True,False)).all(),
            'numeric ROI or same-year VIIRS status missing')
    contrasts=pairs[pairs.cohort.eq('GloBI_moth_record_positive')]
    ids=set(contrasts.land_id)
    roi=old.loc[old.strict9]
    old['cell_id']=old.target_id.str.split('|',regex=False,n=2).str[1]
    roi=old.loc[old.strict9]
    qualified=(new[new.viirs_eligible].groupby(['land_id','year','cell_id'])
               .agg(photos=('photo_id','size'),taxa=('species','nunique'))
               .reset_index())
    qualified=qualified[qualified.photos.ge(5)&qualified.taxa.ge(3)&qualified.land_id.isin(ids)]
    q=set(zip(qualified.land_id,qualified.year,qualified.cell_id))
    attempted=set(zip(old.land_id,old.year,old.cell_id))
    numeric=set(zip(roi.land_id,roi.year,roi.cell_id))
    area=new.groupby('land_id').area_km2.first().astype(float)
    require(ids.issubset(set(area.index)), 'missing signed GSHHG island areas')
    small=[]
    for i in sorted(ids):
        if area[i]>20:continue
        island_cells={x for x in q if x[0]==i}
        old_island_cells={x for x in numeric if x[0]==i}
        small.append(dict(land_id=i,area_km2=float(area[i]),
          tier09_photo_ids=int(new.land_id.eq(i).sum()),
          tier09_qualified_cells=len(island_cells),
          tier08_photo_attempts=int(old.land_id.eq(i).sum()),
          tier08_numeric_ROI9=int(roi.land_id.eq(i).sum()),
          tier09_qualified_cells_not_previously_numeric=len(island_cells-old_island_cells)))
    pair_keys=set(zip(contrasts.land_id,contrasts.year))
    old_years=set(zip(old.land_id,old.year))
    roi_years=set(zip(roi.land_id,roi.year))
    return dict(tier09_light_supported_islands=len(ids),
      tier09_light_supported_island_years=len(pair_keys),
      tier08_attempted_islands=old.land_id.nunique(),
      tier08_genuine_numeric_ROI9_islands=roi.land_id.nunique(),
      tier09_supported_islands_already_attempted=len(ids&set(old.land_id)),
      tier09_supported_islands_with_prior_numeric_ROI9=len(ids&set(roi.land_id)),
      tier09_supported_islands_without_prior_numeric_ROI9=len(ids-set(roi.land_id)),
      tier09_supported_island_years_already_attempted=len(pair_keys&old_years),
      tier09_supported_island_years_with_prior_numeric_ROI9=len(pair_keys&roi_years),
      tier09_qualified_existingVIIRS_cells_in_56_islands=len(q),
      tier09_qualified_cells_already_attempted=len(q&attempted),
      tier09_qualified_cells_with_prior_numeric_ROI9=len(q&numeric),
      tier09_qualified_cells_not_previously_attempted=len(q-attempted),
      tier09_qualified_cells_without_prior_numeric_ROI9=len(q-numeric),
      tier09_small_islands_le20km2=small)

def run(waves,wave_receipt,tier09,tier09_receipt,pairs,out):
    out=Path(out)
    require(not out.exists(),'do not overwrite report')
    w=json.loads(Path(wave_receipt).read_text())
    t=json.loads(Path(tier09_receipt).read_text())
    require(w.get('schema')=='ALAN_real_original_tier08_additional24576_RGB_ROI9_3waves_v1' and
       w.get('real_original_image_photo_attempts')==24576 and
       w.get('actual_numeric_original_9palette_flowerROI_success')==8803 and
       w.get('real_original_technical_unknown_photos')==15773,
       'frozen 24,576 genuine original ROI9 source receipt drift')
    require(t.get('schema')=='ALAN_tier09_moth_sameyear_VIIRS_existingcell_source_gate_v1' and
       t.get('source_roster_photo_count')==78214 and
       t.get('new_RGB_downloads')==0 and t.get('new_numeric_ROI9')==0,
       'frozen tier09 same-year source gate receipt drift')
    require(sha(tier09)==t['candidate_status_SHA256'],'tier09 original roster SHA drift')
    frames=[]
    for i,path in enumerate(waves,1):
        require(sha(path)==w['wave_output_SHA256'][str(i)],
                'real original RGB wave SHA drift')
        d=pd.read_csv(path,compression='gzip',usecols=['photo_id','land_id','target_id','year','strict9'],
           dtype={'photo_id':str,'land_id':str,'target_id':str,'year':str})
        require(len(d)==8192 and d.photo_id.is_unique,'8192 source worker photo IDs changed')
        frames.append(d)
    old=pd.concat(frames,ignore_index=True)
    require(len(old)==24576 and old.photo_id.is_unique and int(old.strict9.sum())==8803,
            '24,576 real measured photo ID/ROI9 count changed')
    n=pd.read_csv(tier09,compression='gzip',
       usecols=['photo_id','land_id','year','cell_id','species','area_km2','viirs_eligible'],
       dtype={'photo_id':str,'land_id':str,'year':str,'cell_id':str,'species':str})
    require(len(n)==78214 and n.photo_id.is_unique,'tier09 photo IDs changed')
    p=pd.read_csv(pairs,dtype={'land_id':str,'year':str})
    require(p.cohort.eq('GloBI_moth_record_positive').sum()==137 and
        p.loc[p.cohort.eq('GloBI_moth_record_positive'),'land_id'].nunique()==56,
        'tier09 verified VIIRS 56-island 137-island-year source changed')
    m=score(old,n,p)
    require(m['tier09_supported_islands_with_prior_numeric_ROI9']==56 and
        m['tier09_supported_islands_without_prior_numeric_ROI9']==0,
        'independent island overlap drift')
    report=dict(schema=SCHEMA,source_only=True,additional_RGB_attempts=0,
      additional_numeric_ROI9=0,full_133_multisite_island_model_NOT_recomputed=True,
      no_causal_ALAN_effect_or_native_white_flower_frequency_claim=True,
      no_moth_obligacy_inferred_from_GloBI=True,
      original_photo_ID_join_only=True,
      tier09_original_status_SHA256=sha(tier09),
      tier08_original_wave_SHA256={str(i):sha(x) for i,x in enumerate(waves,1)},
      result=m)
    out.mkdir(parents=True)
    (out/'tier09_vs_true24576_original_measurements.json').write_text(
      json.dumps(report,sort_keys=True,indent=2)+'\n')
    print('ALAN_TIER09_MEASURED_GEO_OVERLAP',json.dumps(m,sort_keys=True),flush=True)
    return report

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for key in ('wave1','wave2','wave3','wave-receipt','tier09','tier09-receipt','pairs','out'):
        a.add_argument('--'+key,required=True,type=Path,dest=key.replace('-','_'))
    z=a.parse_args()
    run([z.wave1,z.wave2,z.wave3],z.wave_receipt,z.tier09,z.tier09_receipt,z.pairs,z.out)
