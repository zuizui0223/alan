from alan_decadal_events import *
from concurrent.futures import ThreadPoolExecutor
def fruits(y):
    ids=PH.loc[PH.phenophase_name.str.strip().str.lower().eq("ripe fruits"),"phenophase_id"]
    params=dict(start_date=f"{y}-01-01",end_date=f"{y}-12-31",request_src="ALAN fruit bird phenology")
    for k,v in [("station_id",coords.index),("phenophase_id",ids),("species_id",sorted(IDS))]:
        params.update({f"{k}[{j}]":int(x) for j,x in enumerate(v,1)})
    get(C/f"npn_{y}_fruit.json","https://services.usanpn.org/npn_portal/observations/getObservations.json",params,True)
if __name__=="__main__":batch(fruits,range(2008,2025),5)
