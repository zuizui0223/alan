from alan_extended11 import coords, F
import requests, hashlib, rasterio, re, sys
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
C=F/"decades";C.mkdir(exist_ok=True)
OUT=Path(__file__).resolve().parents[2]/"results/decadal_phenology";OUT.mkdir(exist_ok=True)
def fetch(f):
    y=int(re.search(r'19\d\d|20\d\d',f['name'])[0])
    if y==2013 and 'simVIIRS' in f['name']:return []
    path=C/f['name']
    if not path.exists():
        r=requests.get(f['download_url'],timeout=120);r.raise_for_status()
        assert hashlib.md5(r.content).hexdigest()==f['computed_md5']
        path.write_bytes(r.content)
    with rasterio.open(path) as src:
        xy=coords[['longitude','latitude']].to_numpy()
        return [dict(station=st,year=y,light_dn=float(x[0]) if not x.mask.any() else None,land=coords.loc[st,'land_class']) for st,x in zip(coords.index,src.sample(xy,masked=True))]
if __name__=='__main__':
    files=requests.get('https://api.figshare.com/v2/articles/9828827/versions/10',timeout=40).json()['files']
    with ThreadPoolExecutor(5) as ex: rows=sum(ex.map(fetch,files),[])
    d=pd.DataFrame(rows);assert not d.duplicated(['station','year']).any()
    d.to_csv(OUT/'light.csv',index=False)
