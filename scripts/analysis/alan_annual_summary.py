"""Publish a single per-candidate matrix and verify source-to-inference counts."""
from pathlib import Path
import sys
import duckdb
import numpy as np
import pandas as pd

cache,out=map(Path,sys.argv[1:3])
cohorts=["year_span20","year_observed20","screen_span20","screen_observed20","long_term_screen"]
matrix=None
counts=[]
for name in cohorts:
    d=pd.read_csv(out/(name+".csv.gz"),low_memory=False)
    assert d.sp.is_unique
    target=d.q.lt(.05)
    counts.append(dict(analysis=name,roster=len(d),estimated=d.p.count(),negative=(target&d.slope.lt(0)).sum(),
                       positive=(target&d.slope.gt(0)).sum(),minimum_q=d.q.min()))
    numeric=d[["sp","species","coordinates","slope","lo","hi","q"]].rename(columns={k:name+"_"+k for k in ["coordinates","slope","lo","hi","q"]})
    matrix=numeric if matrix is None else matrix.merge(numeric,on=["sp","species"],validate="1:1")
pd.DataFrame(counts).to_csv(out/"model_counts.csv",index=False)
targets=matrix.filter(regex="_q$").lt(.05).any(axis=1)
matrix=matrix[targets].merge(pd.read_csv(out/"candidate_flowering.csv")[["sp","field_years","island_coords","mainland_coords","exact_survey_years"]],validate="1:1")
matrix.to_csv(out/"candidate_summary.csv",index=False)
db=duckdb.connect(str(cache/"annual.duckdb"))
assert db.sql("select sum(visits) from year_surveys").fetchone()==db.sql("select count(*) from survey semi join year_surveys using(site,y)").fetchone()
assert not db.sql("select count(*) from year_hits h join year_surveys s using(site,y) where h.hits>s.visits").fetchone()[0]
assert db.sql("select count(distinct y) from year_surveys").fetchone()[0]==33
assert not db.sql("select count(*) from (select site from year_surveys group by site having max(y)-min(y)<20)").fetchone()[0]
panel=pd.read_parquet(cache/"sweden_matched.parquet").dropna(subset=["gap","dn"])
db.register("swedish_matched",panel)
db.execute("""create or replace table fine_field as select p.*,exists (
     select 1 from geo g join year_surveys s using(site) join year_risk r using(site) join species_names n using(sp)
     where g.latitude=p.latitude and g.longitude=p.longitude and s.y=p.year and n.species=p.name) exact_survey_year
     from swedish_matched p""")
side=panel.groupby(["name","land_class"],as_index=False).agg(coordinates=("unit_id","nunique"),landmasses=("land_id","nunique"),
        station_years=("year","size"),mean_gap=("gap","mean"))
variation=panel.groupby(["name","land_class","unit_id"]).dn.nunique().gt(1).groupby(level=[0,1]).sum().rename("coordinates_with_light_change")
side.merge(variation,on=["name","land_class"]).to_csv(out/"swedish_modis_replication.csv",index=False)
print("Verified species candidate matrix and annual survey denominators",len(matrix))
