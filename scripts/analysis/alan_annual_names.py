"""Resolve field names against the Catalogue of Life checklist used by the source snapshot."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import sys
import time
import requests
import pandas as pd

cache=Path(sys.argv[1])
folder=cache/"taxonomy"
folder.mkdir(exist_ok=True)
checklist="7ddf754f-d193-4cc9-b351-99906754a03b"
def resolve(name):
    path=folder/(hashlib.sha256(name.encode()).hexdigest()+".json")
    if not path.exists():
        for attempt in range(4):
            try:
                response=requests.get("https://api.gbif.org/v2/species/match",params={
                    "scientificName":name,"kingdom":"Plantae","checklistKey":checklist},timeout=20)
                response.raise_for_status()
                result=response.json();path.write_text(json.dumps(result));break
            except requests.RequestException:
                if attempt==3:return dict(field_name=name,match="service_unavailable")
                time.sleep(1+attempt)
    result=json.loads(path.read_text())
    usage=result.get("acceptedUsage") or result.get("usage",{})
    diag=result.get("diagnostics",{})
    match=diag.get("matchType","NONE")
    ok=match=="EXACT" and usage.get("rank")=="SPECIES"
    return dict(field_name=name,match=match,matched_key=usage.get("key") if ok else None,
                accepted_name=usage.get("canonicalName") if ok else None)
names=sorted(pd.read_parquet(cache/"flowering.parquet").species.unique())
with concurrent.futures.ThreadPoolExecutor(8) as pool:
    rows=list(pool.map(resolve,names))
pd.DataFrame(rows).to_csv(cache/"field_crosswalk.csv",index=False)
print(len(rows),flush=True)
