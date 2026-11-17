"""Terra MOD13Q1 QA-screened growing-season time series for mapped inland wetland points."""
import json,math,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import threading
from alan_wet_shore import repo,cache
out=repo/"results/all_species_forest/wet_annual"
folder=cache/"freshwater_modis";folder.mkdir(exist_ok=True)
g=pd.read_csv(out/"habitats.csv")
api="https://planetarycomputer.microsoft.com/api"
bands=["250m_16_days_NDVI","250m_16_days_pixel_reliability","250m_16_days_composite_day_of_the_year"]
local=threading.local()
def get(path,url,params):
    if path.exists():return json.loads(path.read_bytes())
    if not hasattr(local,"s"):
        local.s=requests.Session()
        local.s.mount("https://",HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.6,status_forcelist=[429,500,502,503,504])))
    res=local.s.get(url,params=params,timeout=(12,55));res.raise_for_status()
    j=res.json();path.write_bytes(res.content);return j
g["tile"]=[f"h{int((p.longitude*math.cos(math.radians(p.latitude))+180)//10):02d}v{int((90-p.latitude)//10):02d}" for p in g.itertuples()]
def catalog(job):
    p,y=job
    j=get(folder/f"{p.tile}_{y}.json",api+"/stac/v1/search",dict(collections="modis-13Q1-061",bbox=f"{p.longitude-.001},{p.latitude-.001},{p.longitude+.001},{p.latitude+.001}",datetime=f"{y}-01-01/{y}-12-31",limit=500,fields="id,links"))
    assert not any(l["rel"]=="next" for l in j["links"])
    ids=[f["id"] for f in j["features"] if f["id"].startswith(f"MOD13Q1.A{y}") and f["id"].split(".")[2]==p.tile and 80<=int(f["id"].split(".")[1][5:])<=305]
    assert len(ids)>=13
    return (p.tile,y),ids
with ThreadPoolExecutor(8) as pool:
    catalogs=dict(pool.map(catalog,[(r,y) for r in g.drop_duplicates("tile").itertuples() for y in range(2000,2024)]))
def value(job):
    p,y,item=job
    j=get(folder/f"{p.unit_id}_{item}.json",api+f"/data/v1/item/point/{p.longitude},{p.latitude}",dict(collection="modis-13Q1-061",item=item,assets=bands))
    assert j["band_names"]==[i+"_b1" for i in bands]
    ndvi,qa,day=j["values"]
    return p.unit_id,y,int(item.split(".")[1][5:]),day,qa,ndvi/10000 if ndvi is not None and -2000<=ndvi<=10000 and qa==0 and day is not None and 1<=day<=366 else None
jobs=[(p,y,i) for y in range(2000,2024) for p in g.itertuples() for i in catalogs[p.tile,y]]
with ThreadPoolExecutor(16) as pool:
    rows=list(pool.map(value,jobs))
pd.DataFrame(rows,columns="unit_id y start_doy day qa ndvi".split()).to_csv(out/"vegetation.csv.gz",index=False)
