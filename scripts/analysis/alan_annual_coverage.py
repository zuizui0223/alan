"""Coverage of accepted source histories, including years after the focal species' last detection."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import sys
import duckdb
import geopandas as gpd
import pandas as pd
import requests

cache,out=map(Path,sys.argv[1:3])
db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
st=db.sql("""select g.*,p.datasetkey,p.basisofrecord,min(y) as y0,max(y) as y1,count(distinct y) as years,
             sum(visits) as visits from year_surveys a join geo g using(site) join potential p using(site) group by all""").df()
c=gpd.read_file(cache/"countries.zip")
# Natural Earth is a descriptive overlay; retain unresolved small islands/coasts as unknown.
pts=gpd.GeoDataFrame(st,geometry=gpd.points_from_xy(st.longitude,st.latitude),crs=4326)
world=gpd.sjoin(pts,c[["ADMIN","CONTINENT","geometry"]],how="left",predicate="within")
assert not world.site.duplicated().any()
world["ADMIN"]=world.ADMIN.fillna("unassigned");world["CONTINENT"]=world.CONTINENT.fillna("unassigned")
world.groupby(["CONTINENT","ADMIN","land_class"]).agg(source_stations=("site","size"),coordinates=("unit_id","nunique"),survey_years=("years","sum")).to_csv(out/"countries.csv")
st.groupby("years").agg(stations=("site","size"),coordinates=("unit_id","nunique")).to_csv(out/"history.csv")
st.groupby("basisofrecord").agg(stations=("site","size"),visits=("visits","sum")).to_csv(out/"record_basis.csv")
src=st.groupby("datasetkey",as_index=False).agg(stations=("site","size"),visits=("visits","sum"))
def describe(key):
    try:
        r=requests.get("https://api.gbif.org/v1/dataset/"+key,timeout=20)
        r.raise_for_status();data=r.json()
        return [data["title"],data.get("doi"),data.get("license")]
    except requests.RequestException:return ["","",""]
with ThreadPoolExecutor(6) as pool:details=list(pool.map(describe,src.datasetkey))
src[["title","doi","license"]]=details
src.to_csv(out/"sources.csv",index=False)
s=db.sql("""with x as (select r.sp,r.site,s.y,s.visits,coalesce(h.hits,0) as hits,l.last_seen
 from year_risk r join year_surveys s using(site) join last_seen l using(site,sp)
 left join year_hits h using(site,y,sp))
 select count(*) as species_site_years,count(*) filter(where hits=0) as surveyed_zero,
 count(*) filter(where y>last_seen) as after_last, count(distinct (site,sp)) filter(where y>last_seen) as pairs_with_later_survey
 from x""").df()
s.to_csv(out/"surveyed_nondetections.csv",index=False)
print(s.to_string(index=False));print("datasets",len(src),"assigned countries",world.loc[world.ADMIN.ne("unassigned"),"ADMIN"].nunique())
