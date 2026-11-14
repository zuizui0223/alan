#!/usr/bin/env python3
"""Public primary-source CSV availability and ecological grain, not a SEM refit."""
from __future__ import annotations
import argparse, csv, hashlib, io, json, math, statistics
from pathlib import Path
from scripts.analysis.alan_crossstudy_source_intake_20261114 import get

URL='https://static-content.springer.com/esm/art%3A10.1007%2Fs00114-021-01783-5/MediaObjects/114_2021_1783_MOESM1_ESM.csv'

def main(out):
    raw=get(URL,1_000_000)
    enc=None
    for candidate in ('utf-8-sig','cp1252'):
        try: decoded=raw.decode(candidate);enc=candidate;break
        except UnicodeDecodeError: pass
    if enc is None or '<html' in decoded.lower(): raise ValueError('Not original CSV')
    dialect=csv.Sniffer().sniff(decoded[:5000],delimiters=',;\t')
    reader=csv.DictReader(io.StringIO(decoded),dialect=dialect)
    rows=list(reader)
    if not rows or len(rows)>10000: raise ValueError('Unbounded or empty source')
    schema={}
    for key in reader.fieldnames:
        vals=[r[key] for r in rows if r[key] is not None and r[key].strip() not in {'','NA','NaN','N/A'}]
        nums=[]
        try:
            nums=[float(v) for v in vals]
            if not all(math.isfinite(v) for v in nums): raise ValueError('nonfinite')
        except ValueError: nums=[]
        sensitive=any(s in key.lower() for s in ('lat','lon','coord','site','location','localidad','name'))
        item={'nonmissing_rows':len(vals),'distinct_values':len(set(vals)),
              'numeric':bool(nums),'potential_identifier_or_location':sensitive}
        if nums and not sensitive:
            item.update({'minimum':min(nums),'maximum':max(nums),'mean':statistics.mean(nums),
                         'zeros':sum(x==0 for x in nums)})
        schema[key]=item
    r={'schema':'ALAN_Ceiba2022_primary_csv_intake_v1',
       'article_doi':'10.1007/s00114-021-01783-5',
       'source_url':URL,'bytes':len(raw),'source_sha256':hashlib.sha256(raw).hexdigest(),
       'encoding':enc,'delimiter':dialect.delimiter,'source_rows':len(rows),
       'original_column_summary':schema,'raw_rows_exported':False,
       'study_design':'observational light and noise gradients; not randomized illumination',
       'source_row_ecological_independence_verified':False,'new_model_fitted':False}
    path=Path(out);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():raise ValueError('Refuse overwrite')
    path.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n')
    print('CEIBA_ORIGINAL_CSV_INTAKE',json.dumps(r,ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    main(p.parse_args().out)
