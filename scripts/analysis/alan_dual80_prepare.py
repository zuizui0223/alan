"""Reconstruct the frozen flora cohort and resolve the names used to screen dual-role predators."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import time

import pandas as pd
import requests

from alan_gecko_primary_fiveyear import primary, ROOT, ID


def resolve(name):
    for attempt in range(5):
        try:
            r=requests.get("https://api.gbif.org/v1/species/match",params=dict(name=name,kingdom="Animalia",strict="true"),timeout=45)
            r.raise_for_status()
            t=r.json()
            break
        except requests.RequestException:
            if attempt==4: raise
            time.sleep(2**attempt)
    if t.get("rank")!="SPECIES" or t.get("matchType")!="EXACT" or t.get("canonicalName")!=name:
        return dict(name=name,status="unresolved",response=t)
    return dict(name=name,status=t["status"],key=t.get("acceptedUsageKey",t["usageKey"]),gbif_name=t["species"],response=t)


def run(inputs,signed,private):
    private.mkdir(parents=True,exist_ok=True)
    cohort=pd.read_csv(ROOT/"results/dual_role_80/cohort.csv")
    with ThreadPoolExecutor(max_workers=6) as pool:
        taxa=list(pool.map(resolve,cohort.animal_species))
    (private/"taxa.json").write_text(json.dumps(taxa,indent=2))
    s,digest,plants,legacy=primary(inputs,signed)
    parts=s.target_id.str.split("|",regex=False,expand=True)
    s["land_id"]=parts[2]
    grid=parts[1].str.split("_",expand=True).astype(int)
    s["gx"],s["gy"]=grid[0],grid[1]
    s["x1"],s["y1"]=s.gx//20-180,s.gy//20-90
    s[["key","target_id","year","land_class","land_id","gx","gy","x1","y1"]].to_parquet(private/"sites.parquet",index=False)
    plants[["key",ID,"wcvp_accepted_name"]].to_parquet(private/"plants.parquet",index=False)
    legacy[legacy.guild.eq("moth_family_visit")][["insect_name","plant_name"]].drop_duplicates().to_parquet(private/"visits.parquet",index=False)
    print("target set",digest,"taxa unresolved",[(t["name"],t["status"]) for t in taxa if t["status"]!="ACCEPTED"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs","signed","private"):p.add_argument(arg,type=Path)
    run(**vars(p.parse_args()))
