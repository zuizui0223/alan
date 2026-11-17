"""Read raw CLMS mowing-date bands from its public WMTS at plant survey coordinates."""
import math,base64
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import requests
from rasterio.io import MemoryFile
from alan_wet_shore import repo,cache

out=repo/"results/all_species_forest/wet_annual"
g=pd.read_csv(out/"habitats.csv")
folder=cache/"freshwater_cuts";folder.mkdir(exist_ok=True)
api="https://land.copernicus.eu/cdse/gramd_europe_10m_yearly/"
code='//VERSION=3\nfunction setup(){return {input:["GRAMD1","GRAMD2","GRAMD3","GRAMD4","dataMask"],output:{bands:5,sampleType:"UINT16"}}} function evaluatePixel(s){return [s.GRAMD1,s.GRAMD2,s.GRAMD3,s.GRAMD4,s.dataMask]}'
def value(job):
    r,y=job
    z=15;n=2**z
    x=(r.longitude+180)/360*n
    v=(1-math.asinh(math.tan(math.radians(r.latitude)))/math.pi)/2*n
    path=folder/f"{int(x)}_{int(v)}_{y}.tif"
    if not path.exists():
        p=dict(SERVICE="WMTS",VERSION="1.0.0",REQUEST="GetTile",LAYER="Grassland Mowing Dates 1st event",TILEMATRIXSET="PopularWebMercator256",TILEMATRIX=z,TILECOL=int(x),TILEROW=int(v),TIME=f"{y}-01-01/{y}-12-31",FORMAT="image/tiff",STYLE="default",EVALSCRIPT=base64.b64encode(code.encode()).decode(),UPSAMPLING="NEAREST",DOWNSAMPLING="NEAREST")
        response=requests.get(api,params=p,timeout=75)
        response.raise_for_status()
        assert response.headers["Content-Type"].startswith("image/tiff"),response.text[:200]
        path.write_bytes(response.content)
    with MemoryFile(path.read_bytes()) as mem,mem.open() as ds:
        assert ds.count==5 and ds.dtypes[0]=="uint16"
        a=ds.read()[:,int(v%1*256),int(x%1*256)]
    return [r.unit_id,y,*a]
# Tile-year requests are grouped by year to avoid concurrent writes for near neighbours.
batches=[]
for y in range(2017,2024):
    with ThreadPoolExecutor(5) as pool:batches.extend(pool.map(value,[(p,y) for p in g.itertuples()]))
d=pd.DataFrame(batches,columns=["unit_id","y","cut1","cut2","cut3","cut4","mask"])
d.to_csv(out/"mowing.csv",index=False)
print(d.groupby(["mask","cut1"]).size().to_string())
