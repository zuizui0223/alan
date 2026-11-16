from pathlib import Path
import concurrent.futures, requests, json
import sys
root=Path(sys.argv[1])
def get(y):
 p=dict(start_date=f"{y}-01-01",end_date=f"{y}-12-31",request_src="ALAN annual phenology")
 ids={201,205,210,494,501}
 p.update({f"phenophase_id[{i}]":k for i,k in enumerate(sorted(ids),1)})
 r=requests.get('https://services.usanpn.org/npn_portal/observations/getSummarizedData.json',params=p,timeout=160)
 r.raise_for_status()
 j=r.json()
 print(y,len(j),flush=True)
 assert all(v["phenophase_id"] in ids for v in j)
 (root/f"npn_{y}.json").write_bytes(r.content)
with concurrent.futures.ThreadPoolExecutor(3) as pool:list(pool.map(get,range(2009,2015)))
