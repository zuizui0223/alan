"""Match GLAD annual Landsat water observations to long-running wet-plant surveys."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
import requests
import rasterio
from rasterio.io import MemoryFile
from rasterio.windows import Window


def stations(db,root):
    wet=pd.read_csv(root/"all_species_forest/decline_traits/moisture_evidence.csv.gz").query("indicator > 7")
    db.register("wet",wet[["sp","candidate"]])
    db.execute("""create temp table surveys as select * from year_surveys where y>=1999 and y<=2024
      qualify count(*) over(partition by site)>=5 and max(y) over(partition by site)-min(y) over(partition by site)>=20""")
    db.execute("""create temp table roster as select distinct site,sp,candidate from year_hits join wet using(sp)
                  semi join surveys using(site,y)""")
    db.execute("""create temp table both_sites as select site from roster group by site
                  having bool_or(candidate) and bool_or(not candidate)""")
    return db.sql("""select distinct g.unit_id,g.longitude,g.latitude,s.y from surveys s join both_sites using(site)
                     join geo g on g.site=s.site""").df()


def extract(key,g,cache):
    west,north,year=map(int,key)
    tile=f"{abs(north):02}{'N' if north>=0 else 'S'}_{abs(west):03}{'E' if west>=0 else 'W'}"
    file=cache/f"{tile}_{year}.csv"
    if file.exists():
        stored=pd.read_csv(file)
        if set(stored.unit_id)==set(g.unit_id):return stored
    url=f"https://storage.googleapis.com/earthenginepartners-hansen/waterC2/{tile}/{year}_percent.tif"
    # Isolate each download in memory. Read geographic windows; never load the full decompressed raster.
    for attempt in range(4):
        try:
            response=requests.get(url,timeout=(15,120))
            if response.status_code==404:
                result=g[["unit_id","y"]].assign(water_area=np.nan,water_freq=np.nan,valid=0.,status="tile_unavailable")
                break
            response.raise_for_status()
            with MemoryFile(response.content) as mem,mem.open() as ds:
                assert ds.nodata==255 and abs(ds.transform.c-west)<1e-6 and abs(ds.transform.f-north)<1e-6
                rows=[]
                for p in g.itertuples():
                    rr,cc=ds.index(p.longitude,p.latitude)
                    ry=round(.5/(111.32*abs(ds.transform.e)))
                    rx=round(.5/(111.32*np.cos(np.deg2rad(p.latitude))*ds.transform.a))
                    if rr-ry<0 or cc-rx<0 or rr+ry>=ds.height or cc+rx>=ds.width:
                        rows.append((p.unit_id,p.y,np.nan,np.nan,0.,"tile_edge"));continue
                    a=ds.read(1,window=Window(cc-rx,rr-ry,2*rx+1,2*ry+1))
                    v=a[a<=100]
                    frac=len(v)/a.size
                    rows.append((p.unit_id,p.y,100*(v>=50).mean() if frac>=.8 else np.nan,
                                 v.mean() if frac>=.8 else np.nan,frac,"ok" if frac>=.8 else "low_coverage"))
            result=pd.DataFrame(rows,columns=["unit_id","y","water_area","water_freq","valid","status"])
            break
        except (requests.RequestException,rasterio.errors.RasterioIOError):
            if attempt==3:raise
    result.to_csv(file,index=False)
    return result


def acquire(work,root):
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    points=stations(db,root)
    cache=work/"gbif_annual/glad_water";cache.mkdir(exist_ok=True)
    points["west"]=np.floor(points.longitude/10)*10
    points["north"]=(np.floor(points.latitude/10)+1)*10
    groups=list(points.groupby(["west","north","y"]))
    print(len(points),len(groups),"tile-years",flush=True)
    with ThreadPoolExecutor(6) as pool:
        rows=list(pool.map(lambda pair:extract(*pair,cache),groups))
    water=pd.concat(rows)
    assert not water.duplicated(["unit_id","y"]).any() and len(water)==len(points)
    water.to_parquet(cache/"survey_water.parquet",index=False)
    print(water.status.value_counts(),flush=True)


if __name__=="__main__":
    import sys
    acquire(*map(Path,sys.argv[1:3]))
