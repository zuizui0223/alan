"""Retrieve public flowering and MODIS records for paired long-running NPN stations."""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

STATIONS = {19586:76, 2:76, 25009:312, 33137:312}
YEARS = range(2012,2026)


def download(folder):
    base = folder / "long_ndvi"
    base.mkdir(exist_ok=True)
    visits = pd.read_parquet(folder/"npn/individual_visits.parquet")
    phases = pd.read_json(folder/"npn/getPhenophases.json")
    phases = phases.loc[phases.phenophase_name.str.lower().str.strip().eq("open flowers"),"phenophase_id"].tolist()
    coords = visits.groupby("station")[["latitude","longitude"]].first()
    def fetch(job):
        station, year, product, first, last = job
        path = base / (f"{product}_{station}_{year}_{first}.json" if product=="modis" else f"{product}_{station}_{year}.json")
        if product=="modis" and all(any(base.glob(f"ndvi_{station}_MOD13Q1.A{year}{day:03}.*.json")) for day in range(first,last+1,16)):
            return
        if path.exists():
            return
        if product == "flower":
            params = {"request_src":"ALAN long term flowering research","start_date":f"{year}-01-01",
                      "end_date":f"{year}-12-31","station_id[1]":station,"species_id[1]":STATIONS[station]}
            params.update({f"phenophase_id[{i}]":v for i,v in enumerate(phases,1)})
            r = requests.post("https://services.usanpn.org/npn_portal/observations/getObservations.json",data=params,timeout=65)
            r.raise_for_status()
            records=r.json()
            assert isinstance(records,list) and all(int(x["site_id"])==station and int(x["species_id"])==STATIONS[station] for x in records)
        else:
            lat,lon=coords.loc[station]
            r=requests.get("https://modis.ornl.gov/rst/api/v1/MOD13Q1/subset",
                params=dict(latitude=lat,longitude=lon,startDate=f"A{year}{first:03}",endDate=f"A{year}{last:03}",kmAboveBelow=0,kmLeftRight=0),timeout=100)
            r.raise_for_status()
            records=r.json()
            assert records["nrows"]==records["ncols"]==1 and records["subset"]
        path.write_bytes(r.content)
    def attempt(job):
        try: fetch(job)
        except Exception as exc: print(job,str(exc)[:150],flush=True)
    with ThreadPoolExecutor(8) as pool:
        list(pool.map(attempt,[(s,y,p,a,b) for y in YEARS for s in STATIONS for p,a,b in [("flower",0,0),("modis",1,145),("modis",161,305),("modis",321,353)]]))


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("folder",type=Path)
    download(ap.parse_args().folder)
