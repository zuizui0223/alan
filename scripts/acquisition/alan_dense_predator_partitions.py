"""Collect a high-volume species in complete bounded spatial partitions of the plant sampling frame."""
import argparse
import concurrent.futures
import gzip
import json
import threading
from pathlib import Path

import requests

from alan_gecko_gbif_fiveyear import URL, fetch_query, page

local = threading.local()


def session():
    if not hasattr(local, "session"):
        local.session = requests.Session()
    return local.session


def collect_dense(root, scope, year, taxon):
    root.mkdir(parents=True, exist_ok=True)
    dest = root / f"{year}_{taxon}.json.gz"
    if dest.exists():
        return
    cache = root / f"parts_{year}_{taxon}"
    cache.mkdir(exist_ok=True)
    params = dict(taxonKey=taxon, year=year, hasCoordinate="true", basisOfRecord="HUMAN_OBSERVATION", limit=300)
    cells = set(map(tuple,scope["one_degree"][str(year)]))
    level = [(x,y,5) for x,y in scope["blocks"][str(year)]]
    queries, counts, records = [], [], {}

    def part(box):
        x,y,side = box
        q = {**params,"decimalLongitude":f"{x:.8f},{x+side:.8f}","decimalLatitude":f"{y:.8f},{y+side:.8f}"}
        file = cache / f"{x:g}_{y:g}_{side:g}.json.gz"
        if file.exists():
            with gzip.open(file,"rt") as src:
                return q,json.load(src),[]
        total = page(session(),{**q,"limit":0})["count"]
        if total>2000 and side>1/128:
            if side==5:
                children = [(a,b,1) for a,b in sorted(cells) if x<=a<x+5 and y<=b<y+5]
            else:
                h=side/2
                children=[(x+i,y+j,h) for i in (0,h) for j in (0,h)]
            return q,None,children
        rows = fetch_query(session(),q) if total else []
        assert len(rows)==total
        with gzip.open(file.with_suffix(".tmp"),"wt") as src:
            json.dump(rows,src)
        file.with_suffix(".tmp").replace(file)
        return q,rows,[]

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        while level:
            next_level=[]
            for query, rows, children in pool.map(part,level):
                if children:
                    next_level.extend(children)
                else:
                    queries.append(query)
                    counts.append(len(rows))
                    for row in rows:
                        if row["key"] in records:
                            assert records[row["key"]]==row
                        records[row["key"]]=row
            print(year, "boxes",len(queries),"rows",len(records),"remaining boxes",len(next_level),flush=True)
            level=next_level
    with gzip.open(dest.with_suffix(".tmp"),"wt") as f:
        json.dump(dict(url=URL,query=params,queries=queries,query_counts=counts,count=len(records),
                       spatial_tree=True,rows=list(records.values())),f)
    dest.with_suffix(".tmp").replace(dest)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("private",type=Path)
    p.add_argument("scope",type=Path)
    p.add_argument("year",type=int,nargs="+")
    p.add_argument("--taxon",type=int,default=1311477)
    a=p.parse_args()
    d=json.loads(a.scope.read_text())
    assert d["target_hash"]=="37f3c2c519871b6905ffa81c5a14536b248611861ad143b1d4d77d95b7c6b00d"
    for year in a.year:
        collect_dense(a.private,d,year,a.taxon)
