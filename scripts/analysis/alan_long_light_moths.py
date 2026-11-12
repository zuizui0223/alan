"""Annual four-month VIIRS and nearby dated GBIF adult records of the linked moths."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from alan_long_ndvi_download import STATIONS, YEARS
from mois_global_multipoint_numeric_20261107 import MONTHS, VIIRS, api_params, daterange, extract_samples


def run(folder):
    cache = folder/"long_ndvi"
    visits = pd.read_parquet(folder/"npn/individual_visits.parquet")
    points = visits.groupby("station")[["longitude","latitude"]].first().loc[list(STATIONS)]
    points = points.rename(columns={"longitude":"query_longitude","latitude":"query_latitude"})
    rows = points.to_dict("records")
    def light(job):
        y,m=job
        params=api_params("viirs",m,rows)
        params["time"]=daterange(y,m)
        path=cache/f"light_{y}_{m}.json"
        if not path.exists():
            r=requests.get(VIIRS,params=params,timeout=60)
            r.raise_for_status()
            if r.json().get("error"):
                print(y,m,r.json()["error"],flush=True)
            else:
                extract_samples("viirs",m,rows,r.json(),y)
            path.write_bytes(r.content)
    def moth(station):
        spp=["Spoladea recurvalis","Pyreferra citrombra"] if STATIONS[station]==76 else ["Ctenucha rubroscapus"]
        lon,lat=points.loc[station]
        for name in spp:
            path=cache/f"moth_{station}_{name.replace(' ','_')}.json"
            if path.exists(): continue
            r=requests.get("https://api.gbif.org/v1/occurrence/search",params=dict(scientificName=name,
                decimalLongitude=f"{lon-.025},{lon+.025}",decimalLatitude=f"{lat-.025},{lat+.025}",
                year="2012,2025",hasCoordinate="true",lifeStage="Adult",limit=300),timeout=60)
            r.raise_for_status()
            assert r.json()["endOfRecords"]
            path.write_bytes(r.content)
    with ThreadPoolExecutor(5) as pool:
        list(pool.map(moth, STATIONS))
        list(pool.map(light,[(y,m) for y in YEARS for m in MONTHS]))


if __name__=="__main__":
    run(Path(sys.argv[1]))
