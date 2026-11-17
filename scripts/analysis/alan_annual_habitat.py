from pathlib import Path
import sys,json
import requests
import pandas as pd
import numpy as np
import duckdb
from shapely.geometry import shape,Point
from concurrent.futures import ThreadPoolExecutor
from alan_wet_shore import frozen,cache,out as shore

out=shore.parent/"wet_annual"
out.mkdir(exist_ok=True)
db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
g=frozen("coordinate_trends.csv.gz").merge(db.sql("select distinct unit_id,longitude,latitude from geo").df(),validate="1:1")
pts=g[g.longitude.between(-25,45)&g.latitude.between(30,75)].sort_values(["longitude","latitude"])
url="https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2000_WM/MapServer/0/query"
def identify(part):
    args=dict(where="code_00 in ('411','412')",geometryType="esriGeometryMultipoint",geometry=json.dumps({"points":part[["longitude","latitude"]].values.tolist(),"spatialReference":{"wkid":4326}},separators=(",",":")),inSR=4326,outSR=4326,spatialRel="esriSpatialRelIntersects",outFields="code_00,OBJECTID",returnGeometry="true",f="geojson")
    res=requests.get(url,params=args,timeout=90);res.raise_for_status()
    j=res.json()
    assert "features" in j and not j.get("exceededTransferLimit")
    records=[]
    geoms=[(shape(f["geometry"]),f["properties"]) for f in j["features"]]
    for r in part.itertuples():
        for poly,prop in geoms:
            if poly.covers(Point(r.longitude,r.latitude)):
                records.append((r.unit_id,r.longitude,r.latitude,prop["code_00"],prop["OBJECTID"]))
    return records
with ThreadPoolExecutor(5) as pool:
    rows=[r for batch in pool.map(identify,[pts.iloc[i:i+35] for i in range(0,len(pts),35)]) for r in batch]
pd.DataFrame(rows,columns="unit_id longitude latitude habitat polygon".split()).drop_duplicates().to_csv(out/"habitats.csv",index=False)
print("screen",len(g),len(pts),len(rows),g[~g.unit_id.isin(pts.unit_id)][["longitude","latitude"]].describe())
