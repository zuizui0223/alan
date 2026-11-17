"""Baseline distance to a mapped open-water pixel, restricted to the paired survey panel."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import io
import subprocess
import duckdb
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
import requests

repo=Path(__file__).resolve().parents[2]
cache=repo.parent/"gbif_annual"
out=repo/"results/all_species_forest/wet_shore"
out.mkdir(exist_ok=True)
def frozen(name):
    data=subprocess.check_output(["git","show","6a9a484:results/all_species_forest/wet_water/"+name],cwd=repo)
    return pd.read_csv(io.BytesIO(data),compression="gzip" if name.endswith(".gz") else None)
def map_tile(item):
    (west,north),sites=item
    tile=f"{abs(north):02}{'N' if north>=0 else 'S'}_{abs(west):03}{'E' if west>=0 else 'W'}"
    path=cache/"glad_water"/(tile+"_1999_percent.tif")
    if not path.exists():
        url=f"https://storage.googleapis.com/earthenginepartners-hansen/waterC2/{tile}/1999_percent.tif"
        with requests.get(url,timeout=(20,180)) as res:
            if res.status_code==404:
                return [(i,np.nan,np.nan,"missing_tile") for i in sites.unit_id]
            res.raise_for_status()
            path.write_bytes(res.content)
    rows=[]
    with rasterio.open(path) as ds:
        for p in sites.itertuples():
            lat,lon=p.latitude,p.longitude
            rr,cc=ds.index(lon,lat)
            ky=111320; kx=ky*np.cos(np.deg2rad(lat))
            ry=int(np.ceil(1000/(ky*abs(ds.transform.e))))+1
            rx=int(np.ceil(1000/(kx*ds.transform.a)))+1
            if min(rr-ry,cc-rx)<0 or rr+ry>=ds.height or cc+rx>=ds.width:
                rows.append((p.unit_id,np.nan,np.nan,"tile_edge"));continue
            w=Window(cc-rx,rr-ry,2*rx+1,2*ry+1)
            z=ds.read(1,window=w)
            dx=(ds.transform.c+(cc+np.arange(-rx,rx+1)+.5)*ds.transform.a-lon)*kx
            dy=(ds.transform.f+(rr+np.arange(-ry,ry+1)+.5)*ds.transform.e-lat)*ky
            dist=np.hypot(dy[:,None],dx)
            disk=dist<=1000
            valid=(z<=100)&disk
            frac=valid.sum()/disk.sum()
            wet=valid&(z>=50)
            distance=dist[wet].min() if wet.any() else 1000.
            rows.append((p.unit_id,distance if frac>=.95 else np.nan,frac,"ok" if frac>=.95 else "low_coverage"))
    return rows

if __name__=="__main__":
    g=frozen("coordinate_trends.csv.gz")
    db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
    xy=db.sql("select distinct unit_id,longitude,latitude from geo").df().merge(g[["unit_id"]],validate="1:1")
    xy["west"]=(np.floor(xy.longitude/10)*10).astype(int)
    xy["north"]=((np.floor(xy.latitude/10)+1)*10).astype(int)
    with ThreadPoolExecutor(4) as pool:
        rows=[r for batch in pool.map(map_tile,xy.groupby(["west","north"])) for r in batch]
    result=pd.DataFrame(rows,columns=["unit_id","distance_m","valid_fraction","status"])
    result.to_csv(out/"distance.csv",index=False)
    print(result.status.value_counts())
