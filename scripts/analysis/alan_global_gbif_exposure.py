"""Match same-visit-year satellite brightness to permanent GBIF coordinates and classify land."""
from pathlib import Path
import sys
import duckdb
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
from pyproj import Transformer
from mois_global_2020_GSHHG_scope_20261107 import prepare_land,assign_points

cache,lights,coast,out=map(Path,sys.argv[1:5])
db=duckdb.connect(str(cache/"global_gbif.duckdb"),config={"threads":3,"memory_limit":"1400MB"})
places=db.sql("""select site, datasetkey,basisofrecord, decimallatitude as latitude,decimallongitude as longitude from repeat_places
 semi join repeated using(site)""").df()
coords=places[["longitude","latitude"]].drop_duplicates().reset_index(drop=True)
coords["unit_id"]=coords.index
land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
lakes=gpd.read_file(coast/"GSHHS_h_L2.shp")
assignment=assign_points(coords,land,lakes)
coords=coords.merge(assignment[["unit_id","land_id","land_class"]],validate="1:1")
east,north=Transformer.from_crs(4326,6933,always_xy=True).transform(coords.longitude,coords.latitude)
coords["block"]=(np.floor(east/100000)*1000+np.floor(north/100000)).astype(int)
places=places.merge(coords,validate="m:1")
# Always use the actual years of matched survey dates (not a long-term light climatology).
vis=db.sql("""select site,y,m,count(*) as n,avg(richness) as richness from eligible_visits
 semi join repeated using(site) group by all""").df().merge(places[["site","unit_id"]],validate="m:1")
surfaces=[]
for y,g in vis.groupby("y"):
    path=lights/f"Harmonized_DN_NTL_{y}_{'calDMSP' if y<=2013 else 'simVIIRS'}.tif"
    q=coords.loc[coords.unit_id.isin(g.unit_id.unique()),["unit_id","longitude","latitude"]].copy()
    with rasterio.open(path) as src:
        q=q[q.latitude.gt(src.bounds.bottom)&q.latitude.le(src.bounds.top)&q.longitude.ge(src.bounds.left)&q.longitude.lt(src.bounds.right)].copy()
        if not len(q):continue
        rr,cc=map(np.asarray,rasterio.transform.rowcol(src.transform,q.longitude,q.latitude))
        h,w=src.block_shapes[0]
        values=np.empty(len(q),dtype=float)
        for (i,j),indices in pd.Series(np.arange(len(q))).groupby([rr//h,cc//w]):
            win=src.block_window(1,int(i),int(j))
            grid=src.read(1,window=win,masked=True)
            positions=indices.to_numpy()
            pixels=grid[rr[positions]-int(win.row_off),cc[positions]-int(win.col_off)].astype(float)
            values[positions]=np.ma.filled(pixels,np.nan)
        q["dn"]=values
        surfaces.append(q[["unit_id","dn"]].assign(y=y))
vis=vis.merge(pd.concat(surfaces),how="left",validate="m:1")
# Require measurement of all eligible visit years so the exposure and outcome supports agree.
bad=vis.loc[~vis.dn.between(0,63),"site"].unique()
vis=vis[~vis.site.isin(bad)].copy()
vis["late"]=vis.y.ge(2020)
vis["weighted_light"]=vis.dn*vis.n
vis["weighted_richness"]=vis.richness*vis.n
monthly=vis.groupby(["site","late","m"])[["weighted_light","weighted_richness","n"]].sum()
for v in ["light","richness"]:monthly[v]=monthly["weighted_"+v]/monthly.n
area=monthly.groupby(["site","late"])[["light","richness"]].mean().unstack()
values=pd.DataFrame({"light0":area.light[False],"delta":area.light[True]-area.light[False],
                      "rich_change":np.log(area.richness[True])-np.log(area.richness[False])})
geo=places.merge(values,on="site",validate="1:1")
db.register("geotable",geo)
db.execute("create or replace table site_exposure as select * from geotable")
db.execute("checkpoint")
geo.groupby("land_class").agg(source_locations=("site","size"),coordinates=("unit_id","nunique"),blocks=("block","nunique"),landmasses=("land_id","nunique")).to_csv(out/"geography.csv")
print("light / survey locations",len(geo),len(places),flush=True)
