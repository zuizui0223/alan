#!/usr/bin/env python3
"""Official CORINE change polygons (CHA 2000-06, 06-12, 12-18) at ALAN wetlands.

Never infer land-cover change by subtracting CLC 2000 and 2018 status maps.
The 17 surveyed mapped wetlands are NOT the 187 original matched ALAN pairs.
Outputs aggregate site IDs only, never publishes original survey coordinates.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import pandas as pd
import requests
from pyproj import Transformer
from shapely.geometry import Point, shape
from shapely.ops import transform

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT / "results/all_species_forest/wet_annual"
HAB = DIR / "habitats.csv"
COAST = DIR / "coastal_screen.csv"
PIN_HAB = "6c464c2c2d2a9ddc71a95cab912175c263a90ff0"
EPOCHS = (
    ("2000_2006", "code_00", "code_06"),
    ("2006_2012", "code_06", "code_12"),
    ("2012_2018", "code_12", "code_18"),
)
BASE = "https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine"
WET = {"411", "412"}
RADIUS_M = 1000
OUT = ROOT / "results/all_species_forest/wet_falsification"
TF = Transformer.from_crs(4326, 3035, always_xy=True)
BACK = Transformer.from_crs(4326, 3035, always_xy=True)

def git_blob(path):
    data = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\x00" + data).hexdigest()

def site_source():
    if git_blob(HAB) != PIN_HAB:
        raise RuntimeError("SOURCE_HABITAT_PIN_FAILED")
    p = pd.read_csv(HAB)
    c = pd.read_csv(COAST)
    if len(p) != 18 or p.unit_id.duplicated().any():
        raise RuntimeError("HABITAT_SOURCE_SITE_COUNT_STOP")
    p = p.merge(c, on="unit_id", how="left", validate="1:1")
    p = p[p.coastal_classes_in_2km_square.isna()].copy()
    if len(p) != 17 or set(p.habitat.astype(str)) - WET:
        raise RuntimeError("INLAND_HABITAT_SOURCE_STOP")
    if (~p.latitude.between(-90, 90)).any() or (~p.longitude.between(-180, 180)).any():
        raise RuntimeError("INVALID_ORIGINAL_COORDINATES")
    return p.sort_values("unit_id")

def query_parameters(lon, lat, radius_m=RADIUS_M):
    x, y = TF.transform(lon, lat)
    if not all(math.isfinite(z) for z in (x,y)):
        raise ValueError("INVALID_PROJECTED_COORDS")
    buffer = Point(x, y).buffer(radius_m, resolution=24)
    ring = list(buffer.exterior.coords)
    geom = {"rings": [[list(z) for z in ring]],"spatialReference":{"wkid":3035}}
    params = dict(where="1=1",geometry=json.dumps(geom,separators=(",",":")),
        geometryType="esriGeometryPolygon",inSR=3035,spatialRel="esriSpatialRelIntersects",
        returnGeometry="true",outFields="*",outSR=4326,f="geojson",resultRecordCount=1000)
    return params, buffer

def extract_transition(feature, epoch):
    """Consume actual from/to fields; missing or invalid field is a SOURCE STOP."""
    a,b = epoch[1:]
    props = {str(k).lower():v for k,v in feature["properties"].items()}
    if a not in props or b not in props:
        raise RuntimeError("CORINE_REQUIRED_TRANSITION_FIELDS_MISSING_"+epoch[0])
    src, dst = str(props[a]).strip(),str(props[b]).strip()
    if len(src)!=3 or len(dst)!=3 or not src.isdigit() or not dst.isdigit():
        raise RuntimeError("CORINE_TRANSITION_CODE_INVALID_"+epoch[0])
    return src,dst

def summarize_features(features, circle, epoch):
    """Real intersected polygon area in ETRS89 LAEA, not feature-wide Area_ha."""
    record = dict(epoch=epoch[0],mapped_changes_count=0,
        changed_area_ha=0.0,inland_wetland_loss_ha=0.0,
        inland_wetland_gain_ha=0.0,inland_wetland_to_open_water_ha=0.0,
        inland_wetland_to_artificial_ha=0.0,
        inland_wetland_to_grass_or_heath_ha=0.0,transitions={})
    for f in features:
        src,dst = extract_transition(f,epoch)
        g=f.get("geometry")
        if not g:
            raise RuntimeError("CORINE_MISSING_FEATURE_GEOMETRY")
        projected=transform(BACK.transform,shape(g))
        if not projected.is_valid:
            projected=projected.buffer(0)
        if projected.is_empty:
            raise RuntimeError("CORINE_EMPTY_FEATURE_GEOMETRY")
        clipped=projected.intersection(circle)
        if clipped.is_empty or clipped.area < 1e-7:
            continue
        area=clipped.area/10000.0
        record["mapped_changes_count"]+=1
        record["changed_area_ha"]+=area
        if src in WET and dst not in WET:
            record["inland_wetland_loss_ha"]+=area
            if dst in {"511","512"}: record["inland_wetland_to_open_water_ha"]+=area
            if dst.startswith("1"): record["inland_wetland_to_artificial_ha"]+=area
            if dst in {"231","321","322"}: record["inland_wetland_to_grass_or_heath_ha"]+=area
        if src not in WET and dst in WET:
            record["inland_wetland_gain_ha"]+=area
        key=src+"→"+dst
        record["transitions"][key]=record["transitions"].get(key,0)+area
    area_ha=circle.area/10000
    if record["changed_area_ha"]>area_ha*1.02:
        raise RuntimeError("SOURCE_POLYGON_OVERLAP_AREA_STOP")
    record["mapped_net_inland_wetland_change_ha"]=(
        record["inland_wetland_gain_ha"]-record["inland_wetland_loss_ha"])
    record["survey_circle_area_ha"]=area_ha
    return record

def query_one(args):
    unit,lon,lat,epoch=args
    params,circle=query_parameters(lon,lat)
    url=BASE+"/CHA"+epoch[0]+"_WM/MapServer/0/query"
    err=None
    for attempt in range(3):
        try:
            response=requests.get(url,params=params,timeout=(15,75))
            response.raise_for_status()
            data=response.json()
            if "error" in data or data.get("type")!="FeatureCollection":
                raise RuntimeError("CORINE_API_INVALID_RESPONSE_"+str(data)[:200])
            if data.get("exceededTransferLimit") or len(data.get("features",[])) >= 1000:
                raise RuntimeError("CORINE_QUERY_PAGINATION_STOP")
            rec=summarize_features(data["features"],circle,epoch)
            rec.update(unit_id=int(unit),source_url=url,
                       source_features_sha256=hashlib.sha256(
                           json.dumps(data,sort_keys=True,separators=(",",":")).encode()).hexdigest())
            return rec
        except (requests.exceptions.RequestException,ValueError) as ex:
            err=ex
            if attempt<2: time.sleep(2*(attempt+1))
    raise RuntimeError("CORINE_API_DOWNLOAD_FAILED unit="+str(unit)+" "+epoch[0]+" "+str(err))

def summarize(rows):
    """Areas are descriptive 1-km-circle intersections, not site-level plant change."""
    z=pd.DataFrame(rows)
    if len(z)!=17*len(EPOCHS) or z.duplicated(["unit_id","epoch"]).any():
        raise RuntimeError("CORINE_SITE_EPOCH_COVERAGE_STOP")
    wide=z.groupby("unit_id",sort=True)[["inland_wetland_loss_ha",
        "inland_wetland_gain_ha","mapped_net_inland_wetland_change_ha",
        "changed_area_ha","mapped_changes_count"]].sum().reset_index()
    wetloss=int((wide.inland_wetland_loss_ha>0).sum())
    wetgain=int((wide.inland_wetland_gain_ha>0).sum())
    return dict(site_count=int(len(wide)),site_epochs=int(len(z)),
        sites_with_mapped_wetland_loss=wetloss,
        sites_with_mapped_wetland_gain=wetgain,
        sum_wetland_loss_ha=float(wide.inland_wetland_loss_ha.sum()),
        sum_wetland_gain_ha=float(wide.inland_wetland_gain_ha.sum()),
        sites_with_any_clc_change=int((wide.changed_area_ha>0).sum())),wide

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default=str(OUT/"corine_official_cha_17wetland_20261120.json"))
    p.add_argument("--csv",default=str(OUT/"corine_official_cha_17wetland_20261120.csv"))
    args=p.parse_args()
    sites=site_source()
    jobs=[(r.unit_id,r.longitude,r.latitude,epoch) for r in sites.itertuples() for epoch in EPOCHS]
    with ThreadPoolExecutor(max_workers=4) as pool:
        result=list(pool.map(query_one,jobs))
    total,wide=summarize(result)
    doc={"source":"EEA Copernicus official CLC-Change polygons, not subtracting CLC status",
         "source_original_habitat_git_blob":PIN_HAB,"radius_m":RADIUS_M,
         "n_initial_2000_wetland_sites":17,
         "epochs":[e[0] for e in EPOCHS],
         "summary":total,"status":"SOURCE_VALID_17_INITIAL_INLAND_WETLAND_SITES",
         "not_187_pair_sensitivity":True,
         "no_plant_causal_identification":True,
         "scope_note":"5ha change polygon minimum mapping unit; does not measure small/fine wetland vegetation or hydrology",
         "site_rows":[{**{k:v for k,v in r.items() if k!="transitions"},
                  "transitions":r["transitions"]} for r in result]}
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(doc,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    wide.to_csv(args.csv,index=False)
    print("SOURCE_VERIFIED_CORINE_CHANGE_SUMMARY "+json.dumps(total,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
