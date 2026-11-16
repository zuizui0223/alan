"""Build a worldwide annual database of reconstructed vascular-plant survey days."""
from pathlib import Path
import json
import sys
import duckdb
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from pyproj import Transformer
from mois_global_2020_GSHHG_scope_20261107 import prepare_land, assign_points


def execute(cache, lights, coast, out):
    out.mkdir(parents=True, exist_ok=True)
    folder = cache / "snapshot"
    receipts = sorted(folder.glob("*.json"))
    assert [p.stem for p in receipts] == [f"{i:06}" for i in range(1, 10347)]
    counts = np.array([json.loads(p.read_text()) for p in receipts])
    assert sum(counts[:, 2] > 0) == len(list(folder.glob("*.parquet")))
    db = duckdb.connect(str(cache / "annual.duckdb"), config={"threads":3,"memory_limit":"1600MB","preserve_insertion_order":"false"})
    db.execute(f"create or replace view obs as from read_parquet('{folder.resolve().as_posix()}/*.parquet')")
    keys = "datasetkey,basisofrecord,decimallatitude,decimallongitude"
    # Source and basis must remain the same throughout a station's history.
    steps = {
        "species_names": "select specieskey as sp,min(species) as species from obs group by 1",
        "early_sites": f"select distinct {keys} from obs where year<=2004 and not(month=1 and day=1)",
        "potential":f"""select row_number() over() as site,* from
            (select {keys} from obs semi join early_sites using({keys}) where not(month=1 and day=1)
             group by all having max(year)-min(year)>=20)""",
        "detections": f"""select distinct site,year::smallint as y,month::tinyint as m,day::tinyint as d,specieskey as sp
            from obs join potential using({keys}) where not(month=1 and day=1)""",
        "survey_days":"select site,y,m,d,count(*) as richness from detections group by all having count(*)>=5",
        "stations":"""select site,min(y) as first_year,max(y) as last_year,count(distinct y) as years,count(*) as dates
            from survey_days group by site having max(y)-min(y)>=20""",
        "survey":"select * from survey_days semi join stations using(site)",
        "annual_month":"select site,y,m,count(*) as visits, avg(richness) as richness from survey group by all",
        "hits":"""select site,y,m,sp,count(*) as hits from detections semi join survey using(site,y,m,d) group by all""",
    }
    for name, query in steps.items():
        if not db.execute("select count(*) from information_schema.tables where table_name=?",[name]).fetchone()[0]:
            db.execute("create table "+name+" as "+query)
        print(name,db.sql("select count(*) from "+name).fetchone()[0],flush=True)
    sites=db.sql("""select site,decimallongitude as longitude,decimallatitude as latitude
                    from potential semi join stations using(site)""").df()
    xy=sites[["longitude","latitude"]].drop_duplicates().reset_index(drop=True)
    xy["unit_id"]=xy.index
    land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
    xy=xy.merge(assign_points(xy,land,gpd.read_file(coast/"GSHHS_h_L2.shp"))[["unit_id","land_id","land_class"]],validate="1:1")
    e,n=Transformer.from_crs(4326,6933,always_xy=True).transform(xy.longitude,xy.latitude)
    xy["block"]=(np.floor(e/1e5)*1000+np.floor(n/1e5)).astype(int)
    geo=sites.merge(xy,validate="m:1")
    db.register("geo_frame",geo)
    db.execute("create or replace table geo as from geo_frame")
    sy=db.sql("select distinct unit_id,y from annual_month join geo using(site)").df()
    samples=[]
    for year,g in sy.groupby("y"):
        pts=g.merge(xy,validate="m:1")
        path=next(lights.glob(f"Harmonized_DN_NTL_{year}_*.tif"))
        with rasterio.open(path) as src:
            pts=pts[pts.latitude.gt(src.bounds.bottom)&pts.latitude.le(src.bounds.top)&pts.longitude.ge(src.bounds.left)&pts.longitude.lt(src.bounds.right)].copy()
            if pts.empty:continue
            pts["dn"]=[float(v[0]) if not np.ma.getmaskarray(v)[0] else np.nan for v in src.sample(pts[["longitude","latitude"]].values,masked=True)]
            samples.append(pts[["unit_id","y","dn"]])
    light=pd.concat(samples).query("0<=dn<=63")
    db.register("light_frame",light)
    db.execute("create or replace table light as from light_frame")
    db.execute("""create or replace table panel as select a.*,g.unit_id,g.land_class,g.land_id,g.block,l.dn
              from annual_month a join geo g using(site) join light l using(unit_id,y)
              where land_class in ('mainland','island')""")
    db.execute("""create or replace table eligible as select * from panel semi join
               (select site from panel group by site having max(y)-min(y)>=20) using(site)""")
    stats=dict(snapshot_records=int(counts[:,0].sum()),vascular=int(counts[:,1].sum()),eligible_records=int(counts[:,2].sum()))
    for k in [*steps,"light","panel","eligible"]:
        stats[k]=db.sql("select count(*) from "+k).fetchone()[0]
    (out/"coverage.json").write_text(json.dumps(stats,indent=2))
    db.sql("""select land_class,count(distinct site) source_stations,count(distinct unit_id) coordinates,
            count(distinct land_id) landmasses,count(distinct block) blocks,count(distinct (site,y)) station_years,
            count(distinct site) filter(where years>=20) stations_with_20_sampled_years
            from eligible join stations using(site) group by land_class""").df().to_csv(out/"geography.csv",index=False)
    db.sql("""select y,count(distinct site) as sites,count(distinct unit_id) coordinates,sum(visits) survey_days
             from eligible group by y order by y""").df().to_csv(out/"year_coverage.csv",index=False)
    db.execute("checkpoint")


if __name__=="__main__":
    execute(*map(Path,sys.argv[1:5]))
