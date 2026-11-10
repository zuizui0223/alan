#!/usr/bin/env python3
"""Original photograph months for temporal ALAN controls, source-locked and PRIVATE.

Creates a private year×0.05°cell×landmass seasonal covariate file, and a
public aggregate-only audit. The original source is one authentic photograph
per observation, not independent species abundance or pollination events.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

PHOTO_SHA = 'f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c'
PHOTO_MEMBER = 'canonical_observations_climate.csv.gz'
PHOTO_YEARS = (2019, 2021, 2023, 2024, 2025)

def digest(buf):
    return hashlib.sha256(buf).hexdigest()

def run(canonical_zip, private_site_year, private_out, public_summary,
        *, expected_sha=PHOTO_SHA, synthetic_fixture=False):
    private_out, public_summary = Path(private_out), Path(public_summary)
    if private_out.exists() or public_summary.exists():
        raise ValueError('Refuse overwrite of photographic provenance output')
    if expected_sha != PHOTO_SHA and not synthetic_fixture:
        raise ValueError('Do not substitute the signed source SHA')
    with zipfile.ZipFile(canonical_zip) as z:
        raw = z.read(PHOTO_MEMBER)
    if digest(raw) != expected_sha:
        raise ValueError('Original one-photo-per-observation source SHA mismatch')
    source = pd.read_csv(io.BytesIO(gzip.decompress(raw)),
        usecols=['photo_id', 'observed_on', 'source_subarm',
                 'latitude', 'longitude'],
        dtype={'photo_id':str, 'observed_on':str}, low_memory=False)
    if source.photo_id.duplicated().any():
        raise ValueError('Repeated original photograph source ID')
    sites = pd.read_csv(private_site_year,
        usecols=['target_id','year','photo_attempts'])
    if sites.target_id.duplicated().any() or set(sites.year) != set(PHOTO_YEARS):
        raise ValueError('Original year×photographic site population drift')
    source = source.loc[
        source.source_subarm.eq('historical_main') &
        source.observed_on.str[:4].isin([str(y) for y in PHOTO_YEARS]) &
        source.latitude.between(-89.999, 89.999) &
        source.longitude.between(-179.999, 179.999)].copy()
    n_source_photos = len(source)
    source['year'] = source.observed_on.str[:4].astype(int)
    source['cell'] = (
        np.floor((source.longitude+180)*20).astype(int).astype(str) + '_' +
        np.floor((source.latitude+90)*20).astype(int).astype(str))
    source['join_key'] = source.year.astype(str) + '|' + source.cell
    sites['join_key'] = sites.target_id.str.split('|',regex=False).map(
        lambda x:'|'.join(x[:2]) if len(x)==3 else '')
    if sites.join_key.eq('').any() or sites.join_key.duplicated().any():
        raise ValueError('Duplicated or invalid source year×original cell')
    source = source.loc[source.join_key.isin(set(sites.join_key))].copy()
    source['date'] = pd.to_datetime(source.observed_on,errors='coerce',utc=True)
    if source.date.isna().any():
        raise ValueError('Photo date could not be independently verified')
    month = source.date.dt.month
    source['sin_month'] = np.sin(2*np.pi*(month.astype(float)-.5)/12)
    source['cos_month'] = np.cos(2*np.pi*(month.astype(float)-.5)/12)
    south = source.latitude.lt(0)
    source['south_sin'] = source.sin_month*south
    source['south_cos'] = source.cos_month*south
    source['month'] = month
    agg = source.groupby('join_key',sort=False).agg(
        month_source_photos=('photo_id','size'),
        valid_dates=('date','count'),
        sin_month=('sin_month','mean'),
        cos_month=('cos_month','mean'),
        south_sin=('south_sin','mean'),
        south_cos=('south_cos','mean'),
        mean_month=('month','mean'),
        number_distinct_photo_months=('month','nunique'))
    joined = sites.merge(agg,left_on='join_key',right_index=True,
                         how='left',validate='one_to_one')
    if joined.month_source_photos.isna().any() or not (
        joined.month_source_photos==joined.photo_attempts).all():
        raise ValueError('Photographic month/ROI opportunity source counts disagree')
    private_out.parent.mkdir(parents=True,exist_ok=True)
    joined[['target_id','valid_dates','sin_month','cos_month',
        'south_sin','south_cos','mean_month',
        'number_distinct_photo_months']].to_csv(
        private_out,index=False,compression={'method':'gzip','mtime':0})
    public = {
        'schema':'alan_photo_season_same_original_5years_v1',
        'canonical_original_member_sha256':expected_sha,
        'private_derived_source_sha256':digest(private_out.read_bytes()),
        'original_eligible_source_photos':n_source_photos,
        'matched_original_photo_attempts':len(source),
        'site_years_exact_photo_attempts_reconciled':len(joined),
        'sites_with_more_than_one_month':int(
            joined.number_distinct_photo_months.gt(1).sum()),
        'source_photo_years':list(PHOTO_YEARS),
        'no_original_coordinates_or_species_in_public_receipt':True,
        'seasonality_is_observer_photograph_timing_NOT_pollination_timing':True,
        'synthetic_fixture':bool(synthetic_fixture)}
    public_summary.parent.mkdir(parents=True,exist_ok=True)
    public_summary.write_text(json.dumps(public,ensure_ascii=False,
        indent=2)+'\n',encoding='utf-8')
    return public

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('canonical-zip','private-site-year','private-out','public-summary'):
        p.add_argument('--'+arg,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.canonical_zip,a.private_site_year,
        a.private_out,a.public_summary),ensure_ascii=False))
if __name__=='__main__':main()
