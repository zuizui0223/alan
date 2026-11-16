"""Build a private long-span panel from the MOIS register and sample contemporaneous satellite light."""
from pathlib import Path
import sys
import hashlib
import duckdb
import pandas as pd
import numpy as np
import rasterio

cache, wcvp, rasters, out = map(Path,sys.argv[1:5])
assert len(list(cache.glob("tier-*.parquet")))==224
lookup = wcvp / "mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"
assert hashlib.sha256(lookup.read_bytes()).hexdigest()=="dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2"
con=duckdb.connect(config={"threads":3,"memory_limit":"1200MB"})
taxa=pd.read_csv(lookup).query("wcvp_match_state == 'unique_exact_name_to_accepted_id'")[["inat_taxon_id","wcvp_accepted_plant_name_id"]].astype("int64")
assert taxa.inat_taxon_id.is_unique
con.register("names",taxa)
con.execute(f"create view obs as select * from read_parquet('{cache.as_posix()}/tier-*.parquet')")
# A cell must have a qualifying pre-2005 survey to reach 20 years by 2024.
con.execute("""create temp table early as select y,x from (
select y,x,year from obs join names on taxon=inat_taxon_id where year<=2004
group by all having sum(n)>=5) group by all""")
con.execute("""create temp table members as select year,y,x,wcvp_accepted_plant_name_id as sp,sum(n)::bigint as n
 from obs join early using (y,x) join names on taxon=inat_taxon_id group by all""")
con.execute("""create temp table siteyears as select year,y,x,sum(n)::bigint as effort from members
 group by all having sum(n)>=5""")
con.execute("""create temp table longcells as select y,x from siteyears group by all
 having max(year)-min(year)>=20 and count(*)>=3""")
sites=con.sql("select * from siteyears join longcells using (y,x) order by year,y,x").df()
sites["lon"]=(sites.x+.5)/20-180
sites["lat"]=(sites.y+.5)/20-90
sites["light"]=np.nan
for year,g in sites.groupby("year"):
    sensor="calDMSP" if year<=2013 else "simVIIRS"
    with rasterio.open(rasters/f"Harmonized_DN_NTL_{year}_{sensor}.tif") as src:
        valid=g[g.lat.between(src.bounds.bottom,src.bounds.top)&g.lon.between(src.bounds.left,src.bounds.right)]
        if valid.empty:continue
        rr,cc=map(np.asarray,rasterio.transform.rowcol(src.transform,valid.lon,valid.lat))
        h,w=src.block_shapes[0]
        for (i,j),idx in pd.Series(range(len(valid))).groupby([rr//h,cc//w]):
            window=src.block_window(1,int(i),int(j))
            block=src.read(1,window=window)
            k=idx.to_numpy()
            sites.loc[valid.index[k],"light"]=block[rr[k]-int(window.row_off),cc[k]-int(window.col_off)]
sites=sites.dropna(subset=["light"])
sites=sites.groupby(["y","x"]).filter(lambda s:len(s)>=3 and s.year.max()-s.year.min()>=20)
assert sites.light.between(0,63).all()
con.register("included",sites)
members=con.sql("select m.* from members m join included using (year,y,x)").df()
assert members.n.sum()==sites.effort.sum()
sites.to_parquet(cache/"long_sites.parquet",index=False)
members.to_parquet(cache/"long_members.parquet",index=False)
out.mkdir(parents=True,exist_ok=True)
sites.groupby("year").agg(sampled_cells=("effort","size"),observation_uuids=("effort","sum")).to_csv(out/"annual_coverage.csv")
print("rows, cells, taxa, observations",len(sites),len(sites.drop_duplicates(["x","y"])),members.sp.nunique(),members.n.sum())
