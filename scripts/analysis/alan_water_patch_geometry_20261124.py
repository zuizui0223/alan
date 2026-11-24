"""Local source-gated GLAD open-water patch extraction, NOT executed on real rasters.
python scripts/analysis/alan_water_patch_geometry_20261124.py --geo-db ../gbif_annual/annual.duckdb --raster-dir ../gbif_annual/glad_water --out LOCAL_DIR
Run --self-test to check same amount/distance different geometry, adjacency, missingness.
No coordinates or raw rasters are emitted; never derive them from block or distance.
"""
import argparse
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from scipy import ndimage

COORD_BLOB="1c025bcf6c28a6a3cd125ad05628169e126aef0d"
DIST_BLOB="c2bd22ab52c4674fedfb39a12f67cf847a5e2c3d"
RADII=(250,500,1000)
def patch_metrics(z, disk, pixel_area, threshold=50):
    valid=(z>=0)&(z<=100)&disk
    wet=valid&(z>=threshold)
    n=int(wet.sum()); nv=int(valid.sum()); nd=int(disk.sum())
    rec={'valid_fraction':nv/nd,'water_pixels':n,'water_fraction':n/nv if nv else None}
    for k in (4,8):
        pre=f'c{k}_'
        rec.update({pre+x:None for x in ('patches','largest_share','same_patch_probability','effective_patches','largest_area_m2')})
        if not n or nv!=nd: continue
        lab,nc=ndimage.label(wet,structure=ndimage.generate_binary_structure(2,1 if k==4 else 2))
        sizes=np.bincount(lab[wet],minlength=nc+1)[1:]
        p=sizes/n; prob=float(p@p)
        rec.update({pre+'patches':int(nc),pre+'largest_share':float(p.max()),pre+'same_patch_probability':prob,pre+'effective_patches':1/prob,pre+'largest_area_m2':float(sizes.max()*pixel_area)})
    return rec

def self_test():
    disk=np.ones((9,9),dtype=bool)
    a=np.zeros((9,9),int); b=a.copy()
    a[4,5:9]=80; a[5,5:9]=80
    for i,j in [(4,5),(2,2),(2,6),(0,0),(0,8),(8,0),(8,8),(6,2)]: b[i,j]=80
    ma=patch_metrics(a,disk,900); mb=patch_metrics(b,disk,900)
    assert ma['water_pixels']==mb['water_pixels']==8
    assert ma['water_fraction']==mb['water_fraction']
    assert ma['c4_same_patch_probability']==1 and mb['c4_same_patch_probability']==.125
    assert ma['c8_patches']==1 and mb['c8_patches']==8
    for z in (a,b):
        rr,cc=np.where(z>=50); assert np.hypot(rr-4,cc-4).min()==1
    z=np.zeros((3,3),int); z[0,0]=80;z[1,1]=80
    m=patch_metrics(z,np.ones_like(z,dtype=bool),900)
    assert m['c4_patches']==2 and m['c8_patches']==1
    z[0,1]=255;m=patch_metrics(z,np.ones_like(z,dtype=bool),900)
    assert m['c4_patches'] is None and m['c8_patches'] is None
    m=patch_metrics(np.zeros((2,2),int),np.ones((2,2),bool),900)
    assert m['water_pixels']==0 and m['c4_same_patch_probability'] is None
    print("Synthetic connectivity tests PASS")

def sha256(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()

def run(args):
    import pandas as pd
    import duckdb
    import rasterio
    from rasterio.windows import Window
    coord=pd.read_csv(io.BytesIO(gzip.decompress(subprocess.check_output(["git","show",COORD_BLOB]))))
    dist=pd.read_csv(io.BytesIO(subprocess.check_output(["git","show",DIST_BLOB]))).dropna(subset=["distance_m"])
    target=coord[["unit_id"]].merge(dist,on="unit_id",validate="1:1")
    if len(target)!=1706: raise ValueError("FROZEN_PANEL_COUNT_STOP")
    if not args.geo_db.is_file(): raise FileNotFoundError("ORIGINAL_GEO_DB_REQUIRED")
    with duckdb.connect(str(args.geo_db),read_only=True) as db:
        geo=db.sql("select distinct unit_id,longitude,latitude from geo").df()
    geo=geo[geo.unit_id.isin(target.unit_id)].dropna()
    spread=geo.groupby("unit_id")[["longitude","latitude"]].agg(lambda s:s.max()-s.min())
    if (spread>1e-7).any().any(): raise ValueError("COORDINATE_CONFLICT_STOP")
    geo=geo.drop_duplicates("unit_id")
    p=target.merge(geo,on="unit_id",validate="1:1")
    if len(p)!=1706: raise ValueError(f"ORIGINAL_COORDINATE_SUPPORT_STOP {len(p)}/1706")
    if not p.latitude.between(-90,90).all() or not p.longitude.between(-180,180).all():
        raise ValueError("COORDINATE_RANGE_STOP")
    p["west"]=(np.floor(p.longitude/10)*10).astype(int)
    p["north"]=((np.floor(p.latitude/10)+1)*10).astype(int)
    tasks=[]
    for (west,north),sites in p.groupby(["west","north"]):
        tile=f"{abs(north):02}{'N' if north>=0 else 'S'}_{abs(west):03}{'E' if west>=0 else 'W'}"
        path=args.raster_dir/f"{tile}_1999_percent.tif"
        tasks.append((west,north,tile,path,sites))
    if any(not path.is_file() for _,_,_,path,_ in tasks):
        raise FileNotFoundError("MISSING_CACHED_GLAD_TILES_STOP: "+",".join(tile for _,_,tile,path,_ in tasks if not path.is_file()))
    rows=[]; manifest={}
    for west,north,tile,path,sites in tasks:
        manifest[tile]={"sha256":sha256(path),"source":f"https://storage.googleapis.com/earthenginepartners-hansen/waterC2/{tile}/1999_percent.tif"}
        with rasterio.open(path) as ds:
            if ds.nodata!=255 or ds.crs is None or not ds.crs.is_geographic or abs(ds.transform.c-west)>1e-6 or abs(ds.transform.f-north)>1e-6 or ds.transform.b or ds.transform.d or ds.transform.a<=0 or ds.transform.e>=0:
                raise ValueError("RASTER_GRID_STOP")
            for s in sites.itertuples():
                rr,cc=ds.index(s.longitude,s.latitude)
                ky=111320; kx=ky*math.cos(math.radians(s.latitude))
                ry=math.ceil(1000/(ky*abs(ds.transform.e)))+1
                rx=math.ceil(1000/(kx*ds.transform.a))+1
                if min(rr-ry,cc-rx)<0 or rr+ry>=ds.height or cc+rx>=ds.width:
                    raise ValueError("WINDOW_EDGE_STOP")
                z=ds.read(1,window=Window(cc-rx,rr-ry,2*rx+1,2*ry+1))
                dx=(ds.transform.c+(cc+np.arange(-rx,rx+1)+.5)*ds.transform.a-s.longitude)*kx
                dy=(ds.transform.f+(rr+np.arange(-ry,ry+1)+.5)*ds.transform.e-s.latitude)*ky
                d=np.hypot(dy[:,None],dx); disk=d<=1000
                valid=(z<=100)&disk; wet=valid&(z>=50)
                nearest=float(d[wet].min()) if wet.any() else 1000.
                if valid.sum()/disk.sum()<.95 or abs(nearest-s.distance_m)>1e-3:
                    raise ValueError("BASELINE_DISTANCE_REPLAY_STOP")
                for radius in RADII:
                    for threshold in (50,80):
                        result=patch_metrics(z,d<=radius,abs(ds.transform.a*ds.transform.e)*kx*ky,threshold)
                        # Amount requires >=95% valid. Topology requires 100% valid; unknown is not land.
                        if result["valid_fraction"]<.95:
                            result["water_fraction"]=None
                        rows.append({"unit_id":int(s.unit_id),"year":1999,"radius_m":radius,"threshold":threshold,**result})
    f=pd.DataFrame(rows)
    if len(f)!=1706*len(RADII)*2 or f.duplicated(["unit_id","radius_m","threshold"]).any():
        raise ValueError("OUTPUT_COVERAGE_STOP")
    args.out.mkdir(parents=True,exist_ok=True)
    f.to_csv(args.out/"patch_metrics.csv",index=False)
    summary={"sites":1706,"metrics_rows":len(f),"distance_replayed":True,"raster_manifest":manifest,
             "valid_topology_rows":int(f.c4_same_patch_probability.notna().sum()),"coordinate_rows_published":False,
             "scope":"clipped open-water patch geometry, not hydrological-network or insect functional connectivity"}
    (args.out/"audit.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k!="raster_manifest"}))

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geo-db",type=Path);p.add_argument("--raster-dir",type=Path);p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test: self_test()
    elif not all((a.geo_db,a.raster_dir,a.out)):p.error("geo-db, raster-dir, out required")
    else:run(a)
