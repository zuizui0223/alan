"""Fetch Swiss phenological station observations and annual GIMMS3g SOS (method 4)."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import time

import pandas as pd
import requests

base="https://data.geo.admin.ch/ch.meteoschweiz.ogd-phenology/"
stac="https://data.geo.admin.ch/api/stac/v1/collections/ch.meteoschweiz.ogd-phenology/items"
remote="https://data.globalecology.unh.edu/data/GIMMS_NDVI3g_Phenology/"
cache=Path(sys.argv[1])

def retrieve(job):
    url,path=job
    if path.exists():return
    for attempt in range(4):
        try:
            r=requests.get(url,timeout=120)
            r.raise_for_status()
            if path.suffix==".tif" and r.content[:2] not in [b"II",b"MM"]:raise ValueError("not a TIFF")
            if path.suffix==".csv" and b";" not in r.content[:1024]:raise ValueError("not CSV")
            tmp=path.with_suffix(".partial")
            tmp.write_bytes(r.content)
            tmp.replace(path)
            return
        except requests.RequestException:
            if attempt==3:raise
            time.sleep(2**attempt)
if __name__=="__main__":
    cache.mkdir(exist_ok=True,parents=True)
    jobs=[(base+"ogd-phenology_meta_"+k+".csv",cache/(k+".csv")) for k in ["parameters","stations"]]
    url=stac+"?limit=100"
    while url:
        response=requests.get(url,timeout=30)
        response.raise_for_status()
        page=response.json()
        for item in page["features"]:
            for key,v in item["assets"].items():
                if key.endswith("_y.csv"):jobs.append((v["href"],cache/key))
        url=next((link["href"] for link in page["links"] if link["rel"]=="next"),None)
    for year in range(1992,2015):
        f=f"GIMMS3g_v1_SOS_{year}_Method_4_30_90.tif"
        jobs.append((remote+f,cache/f))
    with ThreadPoolExecutor(10) as pool:list(pool.map(retrieve,jobs))
    print("downloaded",len(jobs))
