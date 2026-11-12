"""Read a single MOD13Q1.061 pixel from the public Planetary Computer COGs."""
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
import requests
from alan_long_ndvi_download import STATIONS

API="https://planetarycomputer.microsoft.com/api"
BANDS=["250m_16_days_NDVI","250m_16_days_pixel_reliability","250m_16_days_composite_day_of_the_year"]


def run(folder):
    cache=folder/"long_ndvi"
    stations=pd.read_parquet(folder/"npn/individual_visits.parquet").groupby("station").first()
    jobs=[]
    for s in STATIONS:
        lat,lon=stations.loc[s,["latitude","longitude"]]
        url=API+"/stac/v1/search"
        params={"collections":"modis-13Q1-061","bbox":f"{lon-.00001},{lat-.00001},{lon+.00001},{lat+.00001}",
                "datetime":"2012-01-01/2025-12-31","limit":500,"fields":"id,links"}
        items=[]
        for year in range(2012,2026):
            params["datetime"]=f"{year}-01-01/{year}-12-31"
            r=requests.get(url,params=params,timeout=60)
            r.raise_for_status()
            catalog=r.json()
            assert not any(link["rel"]=="next" for link in catalog["links"])
            items.extend(f["id"] for f in catalog["features"])
        tile=f"h{int((lon*math.cos(math.radians(lat))+180)//10):02d}v{int((90-lat)//10):02d}"
        items=[i for i in items if i.startswith("MOD13Q1.") and 2012<=int(i.split(".")[1][1:5])<=2025 and i.split(".")[2]==tile]
        items=list({i.split(".")[1]:i for i in sorted(items)}.values())
        print(s,tile,len(items),min(items),max(items),flush=True)
        assert len(items)==len(set(i.split(".")[1] for i in items))
        jobs.extend((s,lon,lat,i) for i in items)
    def fetch(job):
        s,lon,lat,i=job
        path=cache/f"ndvi_{s}_{i}.json"
        if path.exists(): return
        token=i.split(".")[1]
        doy=int(token[5:])
        chunk=1 if doy<=145 else 161 if doy<=305 else 321
        if (cache/f"modis_{s}_{token[1:5]}_{chunk}.json").exists():return
        for attempt in range(3):
            try:
                r=requests.get(API+f"/data/v1/item/point/{lon},{lat}",params={"collection":"modis-13Q1-061","item":i,"assets":BANDS},timeout=50)
                r.raise_for_status()
                v=r.json()
                assert v["band_names"]==[b+"_b1" for b in BANDS] and v["coordinates"]==[lon,lat]
                path.write_bytes(r.content)
                return
            except requests.RequestException:
                if attempt==2: raise
                time.sleep(2**attempt)
    with ThreadPoolExecutor(10) as pool:
        list(pool.map(fetch,jobs))


if __name__=="__main__":
    run(Path(sys.argv[1]))
