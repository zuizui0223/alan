"""Observe wet-to-dry and dry-to-wet pixel transitions in consecutive Landsat calendar months."""
from concurrent.futures import ThreadPoolExecutor
import time
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from alan_wet_shore import repo,cache
out=repo/"results/all_species_forest/wet_annual"
folder=cache/"freshwater_monthly";folder.mkdir(exist_ok=True)
g=pd.read_csv(out/"habitats.csv")
g["tile"]=np.where(g.longitude<0,"60N_010W","60N_000E")
assert g.latitude.between(50,60).all() and g.longitude.abs().lt(10).all()
def run(job):
    tile,y,points=job
    path=folder/f"{tile}_{y}.csv"
    if path.exists():return pd.read_csv(path)
    prev={};rows=[]
    dates=[(y-1,12)] if y>1999 else []
    for year,month in dates+[(y,m) for m in range(1,13)]:
        url=f"https://storage.googleapis.com/earthenginepartners-hansen/waterC2/individualMonthMaps/10deg/{tile}/{year}_{month:02}_percent.tif"
        for retry in range(3):
            try:
                current={}
                with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",GDAL_HTTP_TIMEOUT="45",CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"),rasterio.open("/vsicurl/"+url) as ds:
                    for p in points.itertuples():
                        rr,cc=ds.index(p.longitude,p.latitude)
                        ry=round(.5/(111.32*abs(ds.transform.e)));rx=round(.5/(111.32*np.cos(np.deg2rad(p.latitude))*ds.transform.a))
                        current[p.unit_id]=ds.read(1,window=Window(cc-rx,rr-ry,2*rx+1,2*ry+1))
                break
            except rasterio.errors.RasterioIOError as err:
                if "404" in str(err):current={};break
                if retry==2:raise
                time.sleep(2+retry)
        if year==y:
            for unit in points.unit_id:
                a=current.get(unit);b=prev.get(unit)
                valid=np.nan;cover=np.nan;joint=np.nan;gain=np.nan;loss=np.nan
                if a is not None:
                    v=a<=100;valid=v.mean()
                    if valid>=.8:cover=100*(a[v]>=50).mean()
                if a is not None and b is not None:
                    v=(a<=100)&(b<=100);joint=v.mean()
                    if joint>=.8:
                        gain=100*((a[v]>=50)&(b[v]<50)).mean()
                        loss=100*((a[v]<50)&(b[v]>=50)).mean()
                rows.append((unit,y,month,valid,cover,joint,gain,loss))
        prev=current
    d=pd.DataFrame(rows,columns="unit_id y m valid cover joint gain loss".split())
    d.to_csv(path,index=False);return d
jobs=[(t,y,s) for t,s in g.groupby("tile") for y in range(1999,2024)]
with ThreadPoolExecutor(10) as e:data=pd.concat(e.map(run,jobs))
assert len(data)==len(g)*25*12 and not data.duplicated(["unit_id","y","m"]).any()
data.sort_values(["unit_id","y","m"]).to_csv(out/"monthly_water.csv.gz",index=False)
