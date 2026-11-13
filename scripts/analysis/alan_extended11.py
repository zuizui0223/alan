from pathlib import Path
import pandas as pd, requests, json, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from alan_long_ndvi_cogs import API,BANDS
from alan_long_ndvi_download import YEARS
from mois_global_multipoint_numeric_20261107 import MONTHS,api_params,daterange,VIIRS,extract_samples
from alan_long_phenology import OUT
F=Path(sys.argv[1]); C=F/"extended11"; C.mkdir(exist_ok=True)
O=OUT.parent/"long_phenology_11";O.mkdir(exist_ok=True)
pairs=pd.read_csv(OUT/"frozen_flowering_pair_screen.csv")
v=pd.read_parquet(F/"npn/individual_visits.parquet")
e=pd.read_parquet(F/"npn_exposure.parquet")
frame=v.merge(e,on=["year","latitude","longitude"]).merge(pairs,on=["plant","region5"])
targets=frame[frame.status.eq(1)][["station","plant","land_class","region5","latitude","longitude"]].drop_duplicates()
assert not targets.duplicated(["plant","station"]).any()
coords=targets.drop_duplicates("station").set_index("station")
names=pd.read_csv(F/"npn/target_names.csv")
ph=pd.read_json(F/"npn/getPhenophases.json")
phases=ph[ph.phenophase_name.str.strip().str.lower().eq("open flowers")].phenophase_id.tolist()
def batch(fn,jobs,n):
    with ThreadPoolExecutor(n) as pool: list(pool.map(fn,jobs))
def get(path,url,params,post=False):
    if path.exists():return json.loads(path.read_bytes())
    for i in range(4):
        try:
            r=requests.post(url,data=params,timeout=180) if post else requests.get(url,params=params,timeout=90)
            r.raise_for_status();d=r.json();path.write_bytes(r.content);return d
        except requests.RequestException:
            if i==3:raise
            time.sleep(2**i)
def flower(job):
    plant,year=job
    pp=targets[targets.plant.eq(plant)]
    ids=names.loc[names.accepted_name.eq(plant),"species_id"].tolist()
    params=dict(start_date=f"{year}-01-01",end_date=f"{year}-12-31",request_src="ALAN flowering and NDVI")
    for key,values in [("station_id",pp.station),("species_id",ids),("phenophase_id",phases)]:
        params.update({f"{key}[{i}]":int(x) for i,x in enumerate(values,1)})
    data=get(C/f"{plant}_{year}.json","https://services.usanpn.org/npn_portal/observations/getObservations.json",params,True)
    assert isinstance(data,list) and all(x["site_id"] in set(pp.station) and x["species_id"] in ids and x["phenophase_id"] in phases for x in data)
def light(job):
    year,month,k=job
    rows=coords.iloc[k:k+80].rename(columns={"longitude":"query_longitude","latitude":"query_latitude"}).to_dict("records")
    params=api_params("viirs",month,rows);params["time"]=daterange(year,month)
    get(C/f"lights_{year}_{month}_{k}.json",VIIRS,params)
if __name__=="__main__":
    print(len(targets),len(coords),flush=True)
    batch(flower,[(p,y) for p in pairs.plant for y in YEARS],6)
    batch(light,[(y,m,k) for y in YEARS for m in MONTHS for k in range(0,len(coords),80)],6)
