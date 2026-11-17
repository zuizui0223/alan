"""Environment Agency OGL groundwater observations near the seven focal GBIF plants."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,time
import requests
import duckdb
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
repo=Path(__file__).resolve().parents[2]
cache=repo.parent/"gbif_annual"
out=repo/"results/all_species_forest/seven_moth"
url="https://environment.data.gov.uk/hydrology/id/stations.json"
catalog=cache/"ea_groundwater_stations.json"
if not catalog.exists():
    response=requests.get(url,params={"observedProperty":"groundwaterLevel","_limit":10000},timeout=90)
    response.raise_for_status()
    catalog.write_bytes(response.content)
g=pd.DataFrame(json.loads(catalog.read_bytes())["items"]).dropna(subset=["lat","long"]).reset_index()
assert len(g)<10000
plants=pd.read_csv(out.parent/"wet_reproduction/species.csv")
ids=pd.read_parquet(cache/"all_species_rf_sites.parquet",columns=["sp","unit_id"])
ids=ids[ids.sp.isin(plants[plants.candidate & plants.documented_transport].sp)].unit_id.unique()
db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
xy=db.sql("select distinct unit_id,latitude,longitude from geo").df()
xy=xy[xy.unit_id.isin(ids)].copy()
dist,idx=BallTree(np.radians(g[["lat","long"]]),metric="haversine").query(np.radians(xy[["latitude","longitude"]]))
xy["km"]=dist.ravel()*6371.009
xy["station"]=g.notation.to_numpy()[idx.ravel()]
xy=xy[xy.km.le(1)]
measures=[]
for r in g[g.notation.isin(xy.station)].itertuples():
    for m in r.measures:
        link=m["@id"].replace("http:","https:")
        if "-gw-" in link and "mAOD-qualified" in link:
            measures.append((r.notation,link))
def years(item):
    station,link=item
    params={"mineq-date":"1992-01-01","max-date":"2025-01-01","_limit":15000}
    if "logged" in link:params["time"]="09:00:00"
    for retry in range(4):
        try:
            res=requests.get(link+"/readings.json",params=params,timeout=75)
            res.raise_for_status()
            items=res.json()["items"]
            assert len(items)<15000
            break
        except requests.RequestException:
            if retry==3:raise
            time.sleep(2**retry)
    a=pd.DataFrame(items)
    if a.empty:return []
    a=a[a.quality.eq("Good")].copy()
    a["year"]=pd.to_datetime(a.date).dt.year
    a["month"]=pd.to_datetime(a.date).dt.month
    monthly=a.groupby(["year","month"]).value.median()
    yearly=monthly.groupby(level=0).agg(["median","count"]).query("count>=9")
    return [dict(station=station,measure=link,y=i,level=v["median"],months=v["count"]) for i,v in yearly.iterrows()]
with ThreadPoolExecutor(7) as pool:
    data=pd.DataFrame([row for batch in pool.map(years,measures) for row in batch])
# One reproducible, longest qualified series per well; never stitch instruments.
choose=data.groupby(["station","measure"]).size().reset_index(name="n").sort_values(["station","n","measure"],ascending=[True,False,True]).drop_duplicates("station")
matched=xy.merge(data.merge(choose,on=["station","measure"]),on="station")
matched.to_csv(out/"groundwater_years.csv.gz",index=False)
xy.to_csv(out/"nearest_wells.csv",index=False)
