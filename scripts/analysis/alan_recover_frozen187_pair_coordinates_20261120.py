#!/usr/bin/env python3
"""Recover genuine original site coordinates for the frozen 187 ALAN pairs.

Uses only original published site-coordinate rows with Git blob pins; never
reconstructs coordinates from the 100-km block or 1-km nearest-water distance.
Original authoritative mapping lives in local gbif_annual/annual.duckdb /geo.
Outputs aggregated completeness, not the potentially sensitive survey geometry.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from alan_light_development_overlap_20261119 import read_source, pair_sites

ROOT=Path(__file__).resolve().parents[2]
SOURCES={
    "results/all_species_forest/wet_annual/habitats.csv":"6c464c2c2d2a9ddc71a95cab912175c263a90ff0",
    "results/all_species_forest/wet_shore/nearest_wells.csv":"142bafb45550f73124e6342d694d2871ad3b5f15",
    "results/all_species_forest/seven_moth/nearest_wells.csv":"f36def0fefdff9ea3c4486183cc0804a84e1116a",
}

def git_blob(path):
    data=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\x00"+data).hexdigest()

def validate_coordinates(df):
    required=["unit_id","longitude","latitude"]
    if not set(required)<=set(df):
        raise ValueError("GEO_SOURCE_REQUIRED_FIELDS_MISSING")
    x=df[required].copy().dropna()
    for col in required:
        x[col]=pd.to_numeric(x[col],errors="raise")
    if x.duplicated("unit_id").any():
        # Identical multi-series coordinates are valid, conflicting IDs are not.
        spread=x.groupby("unit_id")[["latitude","longitude"]].agg(lambda s:s.max()-s.min())
        if (spread>1e-7).any().any():
            raise ValueError("GEO_UNIT_ID_CONFLICT_STOP")
        x=x.drop_duplicates("unit_id")
    if not x.latitude.between(-90,90).all() or not x.longitude.between(-180,180).all():
        raise ValueError("GEO_COORDINATE_RANGE_STOP")
    return x.sort_values("unit_id").reset_index(drop=True)

def known_repository_coordinates():
    tables=[]
    for path,sha in SOURCES.items():
        source=ROOT/path
        if git_blob(source)!=sha:
            raise ValueError("SOURCE_GEOMETRY_SHA_MISMATCH:"+path)
        x=pd.read_csv(source)[["unit_id","longitude","latitude"]].copy()
        x["source"]=path
        tables.append(x)
    u=pd.concat(tables,ignore_index=True)
    # Preserve all contributing sources and validate a single exact original
    # coordinate per survey ID rather than accepting discordant assignments.
    validated=validate_coordinates(u)
    validated["source_count"]=u.groupby("unit_id").source.nunique().reindex(validated.unit_id).to_numpy(int)
    return validated

def matched_ids():
    df,sha=read_source()
    pairs,_,_=pair_sites(df,False)
    if len(pairs)!=187 or pairs.block.nunique()!=30:
        raise ValueError("FROZEN_MATCHING_CHANGED")
    ids=pd.concat([pairs.bright_site_id,pairs.comparison_site_id])
    if ids.nunique()!=374:
        raise ValueError("UNEXPECTED_187_PAIR_SITE_REUSE")
    return pairs,sha

def pair_coverage(pairs,geo):
    if geo.unit_id.duplicated().any():
        raise ValueError("SOURCE_GEO_DUPLICATE_STOP")
    present=set(geo.unit_id)
    b=pairs.bright_site_id.isin(present)
    c=pairs.comparison_site_id.isin(present)
    return {"frozen_pairs":int(len(pairs)),"frozen_blocks":int(pairs.block.nunique()),
      "unique_original_matched_sites":int(pd.unique(pd.concat([
        pairs.bright_site_id,pairs.comparison_site_id])).size),
      "available_original_coordinate_sites":int(len(geo)),
      "matched_bright_arm_coordinates":int(b.sum()),
      "matched_comparison_arm_coordinates":int(c.sum()),
      "matched_pairs_both_coordinates":int((b&c).sum()),
      "matched_pairs_both_coordinate_blocks":int(pairs.loc[b&c,"block"].nunique())}

def from_original_duckdb(path):
    if not path or not Path(path).is_file():
        return None
    import duckdb
    db=duckdb.connect(str(Path(path)),read_only=True)
    try:
        z=db.sql("SELECT DISTINCT unit_id, longitude, latitude FROM geo").df()
    finally:
        db.close()
    return validate_coordinates(z)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-geo-db",default=None,
      help="Local original annual.duckdb only; never upload original survey coordinates")
    ap.add_argument("--out",default="")
    args=ap.parse_args()
    pairs,source_sha=matched_ids()
    known=known_repository_coordinates()
    baseline=pair_coverage(pairs,known)
    original=from_original_duckdb(args.original_geo_db)
    report={
      "source_frozen_panel_sha256":source_sha,
      "recovery_source_git_blobs":SOURCES,
      "publicly_recoverable_coordinate_coverage":baseline,
      "original_local_geo_source_supplied":original is not None,
      "original_local_geo_coverage":pair_coverage(pairs,original) if original is not None else None,
      "original_187_pair_vegetation_possible": (
        bool(original is not None and
             pair_coverage(pairs,original)["matched_pairs_both_coordinates"]==187)),
      "status":"GEO_SOURCE_RECOVERED" if original is not None and
          pair_coverage(pairs,original)["matched_pairs_both_coordinates"]==187
          else "ORIGINAL_GEO_DB_REQUIRED_FOR_FULL_187_PAIR_SPATIAL_OVERLAY",
      "derived_coordinates_from_block_or_nearest_distance":False,
      "survey_coordinates_published_in_result":False
    }
    s=json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)
    print("ALAN_FROZEN187_GEO_RECOVERY_AUDIT\n"+s,flush=True)
    if args.out:
        path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(s+"\n")
if __name__=="__main__":
    main()
