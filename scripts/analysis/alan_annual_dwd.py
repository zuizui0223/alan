"""Read DWD annual reporters' quality-checked wild-plant onset, one latest series per plant."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import re
import sys
import requests
import pandas as pd
from bs4 import BeautifulSoup

BASE="https://opendata.dwd.de/climate_environment/CDC/"
HIST="observations_germany/phenology/annual_reporters/wild/historical/"


def read_dwd(root):
    root.mkdir(exist_ok=True)
    def retrieve(relative):
        path=root/relative.split("/")[-1]
        if not path.exists():
            r=requests.get(BASE+relative,timeout=120);r.raise_for_status()
            p=path.with_suffix(".tmp");p.write_bytes(r.content);p.replace(path)
        return path
    u=requests.get(BASE+HIST,timeout=30);u.raise_for_status()
    files={}
    for a in BeautifulSoup(u.text,"html.parser").select("a[href]"):
        name=a["href"]
        match=re.fullmatch(r"(PH_Jahresmelder_Wildwachsende_Pflanze_.+)_(\d{4})_(\d{4})_hist.txt",name)
        if match and int(match[2])<=2014 and int(match[3])>=1992:
            key=(match[1],match[2])
            if key not in files or name>files[key]:files[key]=name
    def metadata(file):
        txt=retrieve("help/"+file).read_bytes().decode("latin1")
        rec=[[v.strip() for v in row.strip().split(";")] for row in txt.split(";eor;") if row.strip() and row.strip()!="\x1a"]
        return pd.DataFrame([r for r in rec[1:] if len(r)==len(rec[0])],columns=rec[0])
    plants=metadata("PH_Beschreibung_Pflanze.txt")
    stations=metadata("PH_Beschreibung_Phaenologie_Stationen_Jahresmelder.txt")
    plants=plants.rename(columns={"Objekt_ID":"Objekt_id","Objekt_latein":"species"})
    plants["Objekt_id"]=pd.to_numeric(plants.Objekt_id,errors="coerce")
    stations=stations.rename(columns={"geograph.Breite":"latitude","geograph.Laenge":"longitude","Stations_id":"station"})
    for c in ["latitude","longitude","station"]:stations[c]=pd.to_numeric(stations[c],errors="coerce")
    def read(name):
        d=pd.read_csv(retrieve(HIST+name),sep=";",encoding="latin1",low_memory=False)
        d.columns=d.columns.str.strip()
        for col in ["Phase_id","Referenzjahr","Eintrittsdatum_QB","Objekt_id","Stations_id","Jultag","Eintrittsdatum"]:
            d[col]=pd.to_numeric(d[col],errors="coerce")
        d=d[d.Phase_id.eq(5)&d.Referenzjahr.between(1992,2014)&d.Eintrittsdatum_QB.isin([1,2,3])]
        d=d[d.Eintrittsdatum.floordiv(10000).eq(d.Referenzjahr)]
        return d[["Stations_id","Referenzjahr","Objekt_id","Jultag"]]
    with ThreadPoolExecutor(5) as pool:d=pd.concat(pool.map(read,files.values())).drop_duplicates()
    d=d.rename(columns={"Stations_id":"station","Referenzjahr":"year","Jultag":"flower"})
    d=d.merge(plants[["Objekt_id","species"]],validate="m:1").merge(stations[["station","latitude","longitude"]].drop_duplicates(),validate="m:1")
    d=d[d.species.str.fullmatch("[A-Z][a-z]+ [a-z-]+")].assign(source="DWD")
    return d[["source","station","year","latitude","longitude","flower","species"]]


if __name__=="__main__":
    folder=Path(sys.argv[1]);d=read_dwd(folder);d.to_parquet(folder/"onset.parquet",index=False);print(d.groupby("species").size())
