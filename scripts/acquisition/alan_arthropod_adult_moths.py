"""Acquire adult moth records only for 1-degree cell-years in which their predator and plant qualify."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import sys
import threading

import requests
from shapely.geometry import box
from shapely.ops import unary_union
from shapely.geometry.polygon import orient
from shapely.geometry import MultiPolygon
from shapely.wkt import dumps

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"analysis"))
from alan_arthropod_prey_flora import candidates
from alan_gecko_gbif_fiveyear import fetch_query

local = threading.local()


def session():
    if not hasattr(local,"session"):
        local.session=requests.Session()
    return local.session


def bounds(cells):
    shape=unary_union([box(x-0.00001,y-0.00001,x+1.00001,y+1.00001) for x,y in cells])
    if shape.geom_type=="Polygon":
        shape=orient(shape)
    else:
        shape=MultiPolygon([orient(p) for p in shape.geoms])
    return dumps(shape,rounding_precision=5,trim=True)


def collect(job):
    folder,name,taxon,year,cells=job
    query=dict(year=year,taxonKey=taxon,basisOfRecord="HUMAN_OBSERVATION",hasCoordinate="true",lifeStage=["Adult","Imago"],
               coordinateUncertaintyInMeters="0,1000",geometry=bounds(cells),limit=300)
    label=hashlib.sha256(json.dumps(query,sort_keys=True).encode()).hexdigest()[:20]
    path=folder/(label+".json.gz")
    if not path.exists():
        rows=fetch_query(session(),query)
        raw=json.dumps(dict(query=query,rows=rows),sort_keys=True).encode()
        tmp=path.with_suffix(".tmp")
        with gzip.open(tmp,"wb") as stream: stream.write(raw)
        tmp.replace(path)
    with gzip.open(path,"rb") as stream: raw=stream.read()
    o=json.loads(raw)
    assert o["query"]==query
    return dict(moth=name,taxon=taxon,year=year,cells=cells,file=path.name,count=len(o["rows"]),sha256=hashlib.sha256(raw).hexdigest())


def run(private, edges):
    _,_,_,possible=candidates(private,edges)
    needed=possible[possible.scale.eq("1deg")][["moth","year","x1","y1"]].drop_duplicates()
    folder=private/"moth_queries"
    folder.mkdir(exist_ok=True)
    jobs=[]
    for name,g in needed.groupby("moth"):
        response=session().get("https://api.gbif.org/v1/species/match",params={"name":name,"kingdom":"Animalia","strict":"true"},timeout=45)
        response.raise_for_status()
        tax=response.json()
        assert tax.get("matchType")=="EXACT" and tax.get("rank")=="SPECIES" and tax.get("status")=="ACCEPTED" and tax.get("canonicalName")==name,(name,tax)
        for year,group in g.groupby("year"):
            cells=sorted(group[["x1","y1"]].itertuples(index=False,name=None))
            for start in range(0,len(cells),16):
                jobs.append((folder,name,tax["usageKey"],int(year),cells[start:start+16]))
    print("eligible species",needed.moth.nunique(),"cell searches",len(needed),"requests",len(jobs),flush=True)
    with ThreadPoolExecutor(max_workers=5) as pool:
        done=list(pool.map(collect,jobs))
    (private/"moth_queries.json").write_text(json.dumps(done,indent=2),encoding="utf8")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("private",type=Path)
    ap.add_argument("edges",type=Path)
    opts=ap.parse_args()
    run(opts.private,opts.edges)
