"""Read the public species accounts for the frozen cohort and retain the source URL."""
import io
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote

from bs4 import BeautifulSoup
import pandas as pd
import requests

repo=Path(__file__).resolve().parents[2]
cache=repo.parent/"gbif_annual/pladias_wet_accounts.csv"
f=pd.read_csv(io.BytesIO(subprocess.check_output(["git","show","b780370:results/all_species_forest/wet_network/species.csv"],cwd=repo)))
traits=["Pollination syndrome","Reproduction type","Generative reproduction type","Flower colour","Flower symmetry"]
def lookup(row):
    for name in dict.fromkeys([row.accepted,row.species]):
        if pd.isna(name):continue
        url="https://pladias.cz/en/taxon/data/"+quote(name)
        for attempt in range(3):
            try:
                response=requests.get(url,timeout=30)
                if response.status_code==404:break
                response.raise_for_status()
                page=BeautifulSoup(response.content,"html.parser")
                labels=[h.get_text(" ",strip=True) for h in page.select("h1")]
                if name not in labels:break
                values={s.get_text(strip=True):" ".join(s.parent.b.stripped_strings)
                        for s in page.select("span.font-italic") if s.get_text(strip=True) in traits and s.parent.b}
                return dict(sp=row.sp,url=url,account_name=name,status="verified",**values)
            except requests.RequestException:
                if attempt==2:return dict(sp=row.sp,url=url,status="request_failed")
                time.sleep(2**attempt)
    return dict(sp=row.sp,status="no_exact_account")
with ThreadPoolExecutor(6) as pool:
    pd.DataFrame(pool.map(lookup,f.itertuples())).to_csv(cache,index=False)
