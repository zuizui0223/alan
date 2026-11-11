"""Positive island screen and exhaustive candidate-cell occurrence searches on the original GSHHG frame."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from math import floor
from pathlib import Path
import sys

import geopandas as gpd
import pandas as pd
from shapely.geometry import box, Point, MultiPolygon
from shapely.geometry.polygon import orient
from shapely.ops import unary_union
from shapely.wkt import dumps

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"analysis"))
from alan_gecko_primary_fiveyear import prepare_land
from alan_arthropod_adult_moths import session
from alan_gecko_gbif_fiveyear import FIELDS,page

SCALES={"0.05deg":("gx","gy"),"1deg":("x1","y1")}


def geometry(cells,scale):
    rect=[]
    for a,b in cells:
        x,y,step=(a/20-180,b/20-90,.05) if scale=="0.05deg" else (a,b,1)
        rect.append(box(max(-180,x-.00001),max(-90,y-.00001),min(180,x+step+.00001),min(90,y+step+.00001)))
    g=unary_union(rect)
    g=orient(g) if g.geom_type=="Polygon" else MultiPolygon([orient(p) for p in g.geoms])
    return dumps(g,rounding_precision=5,trim=True)


class Geography:
    def __init__(self,coast):
        receipt=json.loads((coast/"source_receipt.json").read_text(encoding="utf8"))
        for file,h in receipt["extracted_sha256"].items():
            assert hashlib.sha256((coast/file).read_bytes()).hexdigest()==h
        land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
        self.land=land.sindex
        self.ids=list(land.land_id)
        self.lakes=gpd.read_file(coast/"GSHHS_h_L2.shp").sindex

    def match(self,r,scale,sites):
        x,y=r["decimalLongitude"],r["decimalLatitude"]
        cell=(floor((x+180)*20),floor((y+90)*20)) if scale=="0.05deg" else (floor(x),floor(y))
        if cell not in sites:return None
        pt=Point(x,y)
        idx=self.land.query(pt,predicate="within")
        assert len(idx)<=1
        if len(idx) and self.ids[idx[0]] in sites[cell] and not len(self.lakes.query(pt,predicate="intersects")):
            return self.ids[idx[0]],cell
        return None


def cached_page(folder,q,offset):
    token=hashlib.sha256(json.dumps([q,offset],sort_keys=True).encode()).hexdigest()
    path=folder/(token+".gz")
    if not path.exists():
        data=page(session(),q,offset)
        assert data["offset"]==offset
        data["results"]=[{k:r.get(k) for k in FIELDS} for r in data["results"]]
        tmp=path.with_suffix(".tmp")
        with gzip.open(tmp,"wt") as f:json.dump(dict(query=q,**data),f)
        tmp.replace(path)
    with gzip.open(path,"rt") as f:data=json.load(f)
    assert data["query"]==q and data["offset"]==offset
    return data,token


def query_cells(folder,geo,year,scale,sites,taxa,screen,adult=False,cell_presence=False):
    """Screen may stop searching a taxon after a positive. Full searches page to exhaustion."""
    if cell_presence:
        assert len(taxa)==1 and not screen
        sites={cell:set(ids) for cell,ids in sites.items()}
    keys=set(taxa)
    geo_query=geometry(sites,scale)
    hits={}
    receipts=[]
    offset=0
    previous=None
    while keys:
        q=dict(taxonKey=sorted(keys),year=int(year),hasCoordinate="true",basisOfRecord="HUMAN_OBSERVATION",
               coordinateUncertaintyInMeters="0,1000",geometry=geo_query,limit=300)
        if adult:q["lifeStage"]=["Adult","Imago"]
        data,token=cached_page(folder,q,offset)
        receipts.append(token)
        if previous is not None:assert previous==data["count"]
        found=set()
        occupied=set()
        for r in data["results"]:
            k=r["speciesKey"]
            if k not in keys or r["species"]!=taxa[k]["gbif_name"]:continue
            if not 0<(r["coordinateUncertaintyInMeters"] or 0)<=1000 or r["occurrenceStatus"] not in ("DETECTED","PRESENT"):continue
            assert r["year"]==year and r["basisOfRecord"]=="HUMAN_OBSERVATION"
            if adult and str(r["lifeStage"]).lower() not in ("adult","imago"):continue
            match=geo.match(r,scale,sites)
            if match:
                land,cell=match
                hits[(taxa[k]["name"],year,land,*cell)]=r["key"]
                found.add(k)
                occupied.add((cell,land))
        if data["endOfRecords"]:
            assert offset+len(data["results"])==data["count"]
            break
        if cell_presence and occupied:
            for cell,land in occupied:
                sites[cell].discard(land)
                if not sites[cell]: del sites[cell]
            if not sites: break
            geo_query=geometry(sites,scale)
            offset,previous=0,None
        elif screen and found:
            keys-=found
            offset,previous=0,None
        else:
            offset+=len(data["results"])
            previous=data["count"]
            if not data["results"] or offset>=100000:raise ValueError("incomplete search or GBIF ceiling")
    return hits,receipts


def chunks(frame,scale,size=12):
    a,b=SCALES[scale]
    for year,g in frame.groupby("year",sort=False):
        cells=g.groupby([a,b]).land_id.agg(set).to_dict()
        pairs=sorted(cells)
        for i in range(0,len(pairs),size):
            yield int(year),{xy:cells[xy] for xy in pairs[i:i+size]}


def screen(private,coast):
    taxa=json.loads((private/"taxa.json").read_text())
    assigned={t["key"]:t for t in taxa if "key" in t}
    assert len(assigned)==sum("key" in t for t in taxa)
    sites=pd.read_parquet(private/"sites.parquet").query("land_class=='island'").sort_values("year",ascending=False)
    geo=Geography(coast)
    folder=private/"pages"
    folder.mkdir(exist_ok=True)
    positives={}
    logs=[]
    for scale in SCALES:
        eligible=dict(positives.get("0.05deg",{}))
        queries=list(chunks(sites,scale))
        with ThreadPoolExecutor(max_workers=5) as pool:
            for i in range(0,len(queries),5):
                active={k:t for k,t in assigned.items() if t["name"] not in eligible}
                if not active:break
                tasks=[pool.submit(query_cells,folder,geo,y,scale,cells,active,True) for y,cells in queries[i:i+5]]
                for future in tasks:
                    hits,refs=future.result()
                    logs.extend(refs)
                    for (name,year,land,x,y),obs in hits.items():
                        eligible.setdefault(name,dict(year=year,land=land,x=x,y=y,occurrence=obs,scale=scale))
                print(scale,i+1,len(queries),"species documented",len(eligible),flush=True)
        positives[scale]=eligible
        (private/"island_screen.json").write_text(json.dumps(dict(positives=positives,pages=logs),indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("private",type=Path);p.add_argument("coast",type=Path)
    screen(**vars(p.parse_args()))
