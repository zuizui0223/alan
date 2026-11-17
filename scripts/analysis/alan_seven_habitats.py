"""Collect species-level narrow habitat optima from previously verified Pladias accounts."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import time
import requests
from bs4 import BeautifulSoup
import pandas as pd

repo=Path(__file__).resolve().parents[2]
root=repo/"results/all_species_forest"
species=pd.read_csv(root/"wet_reproduction/species.csv").dropna(subset=["url","account_name"])
def collect(row):
    for attempt in range(4):
        try:
            res=requests.get(row.url,timeout=40)
            res.raise_for_status()
            break
        except requests.RequestException:
            if attempt==3:raise
            time.sleep(2**attempt)
    soup=BeautifulSoup(res.content,"html.parser")
    assert row.account_name in [h.get_text(" ",strip=True) for h in soup.select("h1")]
    labels=[]
    for field in soup.select("span.font-italic"):
        title=field.get_text(" ",strip=True)
        if not re.match(r"^\d+[A-Z] ",title) or not field.parent.b:continue
        value=field.parent.b.get_text(" ",strip=True)
        if "optimum" in value:
            labels.append((row.sp,title,row.url))
    return labels
with ThreadPoolExecutor(6) as pool:
    records=[r for rows in pool.map(collect,species.itertuples()) for r in rows]
pd.DataFrame(records,columns=["sp","optimum","source"]).to_csv(root/"seven_moth/habitat_optima.csv.gz",index=False)
