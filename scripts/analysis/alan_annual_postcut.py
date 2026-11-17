"""Sentinel-2 L2A greenness at survey points with first-cut dates in CLMS."""
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import numpy as np
import requests
import xml.etree.ElementTree as ET
from alan_wet_shore import repo,cache
out=repo/"results/all_species_forest/wet_annual"
cuts=pd.read_csv(out/"mowing.csv").query("mask==1 and cut1>=1 and cut1<=366").merge(pd.read_csv(out/"habitats.csv"),on="unit_id")
api="https://planetarycomputer.microsoft.com/api"
jobs=[]
for r in cuts.itertuples():
    date=pd.Timestamp(r.y,1,1)+pd.Timedelta(days=r.cut1-1)
    lon,lat=r.longitude,r.latitude
    params=dict(collections="sentinel-2-l2a",bbox=f"{lon-.00001},{lat-.00001},{lon+.00001},{lat+.00001}",datetime=f"{date-pd.Timedelta(days=28):%Y-%m-%d}/{date+pd.Timedelta(days=42):%Y-%m-%d}",limit=250)
    j=requests.get(api+"/stac/v1/search",params=params,timeout=60).json()
    assert not any(l["rel"]=="next" for l in j["links"])
    jobs.extend((r,date,f) for f in j["features"])
def sample(job):
    r,date,f=job
    id=f["id"];props=f["properties"]
    when=pd.Timestamp(props["datetime"]).tz_convert(None).normalize()
    response=requests.get(api+f"/data/v1/item/point/{r.longitude},{r.latitude}",params=dict(collection="sentinel-2-l2a",item=id,assets=["B04","B08","SCL"]),timeout=60);response.raise_for_status()
    j=response.json();assert j["band_names"]==["B04_b1","B08_b1","SCL_b1"]
    red,nir,scl=j["values"];offset_red=offset_nir=0
    # Sentinel-2 MSI L2A DQR: BOA_ADD_OFFSET is -1000 for all spectral bands from PB 04.00.
    if float(props["s2:processing_baseline"])>=4:
        offset_red=offset_nir=-1000
    val=np.nan
    if scl in [4,5] and None not in [red,nir]:
        a=red+offset_red;b=nir+offset_nir
        if min(a,b)>=0 and a+b>0:val=(b-a)/(a+b)
    return r.unit_id,r.y,r.cut1,when.date(),(when-date).days,id,scl,val
with ThreadPoolExecutor(8) as pool:res=pd.DataFrame(pool.map(sample,jobs),columns="unit_id y cut1 date since_cut item scl ndvi".split())
res.to_csv(out/"postcut_sentinel.csv.gz",index=False)
d=res.dropna(subset="ndvi").groupby(["unit_id","y","since_cut"],as_index=False).ndvi.median()
d["window"]=pd.cut(d.since_cut,[-29,-1,14,42],labels=["before","early","late"])
summary=d.groupby(["unit_id","y","window"],observed=True).ndvi.agg(["median","count"])
print(summary)
summary.to_csv(out/"postcut_windows.csv")
# Resolve within-year ordering using the original survey day.
import duckdb
survey=pd.read_csv(out/"survey_series.csv.gz")
cut_sites=cuts.merge(survey[["unit_id","site"]].drop_duplicates(),on="unit_id")
db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
db.sql("""select c.unit_id,c.y,c.cut1,s.m,s.d from cut_sites c left join survey_days s using(site,y)
 order by c.y,s.m,s.d""").df().to_csv(out/"cut_survey_dates.csv",index=False)
