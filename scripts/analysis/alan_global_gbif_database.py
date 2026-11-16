"""Build the private analytic database and the source/effort denominators."""
from pathlib import Path
import sys,json
import duckdb
import numpy as np
import pandas as pd

cache,out=map(Path,sys.argv[1:3])
out.mkdir(exist_ok=True,parents=True)
files=sorted((cache/"snapshot").glob("*.json"))
assert len(files)==10346 and [p.stem for p in files]==[f"{i:06}" for i in range(1,10347)]
counts=np.array([json.loads(p.read_text()) for p in files])
assert np.sum(counts[:,2]>0)==len(list((cache/"snapshot").glob("*.parquet")))
db=duckdb.connect(str(cache/"global_gbif.duckdb"),config={"threads":3,"memory_limit":"1500MB","preserve_insertion_order":"false"})
db.execute(f"create or replace view occurrences as select * from read_parquet('{(cache/'snapshot'/'*.parquet').resolve().as_posix()}')")
steps={
"names":"""select specieskey, min(species) as species,count(distinct species) as spellings from occurrences group by specieskey""",
"early_places":"""select distinct datasetkey,basisofrecord,decimallatitude,decimallongitude from occurrences
 where year<=1999 and not(month=1 and day=1)""",
"repeat_places":"""select row_number() over () as site, * from (select distinct datasetkey,basisofrecord,decimallatitude,decimallongitude from occurrences
 semi join early_places using(datasetkey,basisofrecord,decimallatitude,decimallongitude)
 where year>=2020 and not(month=1 and day=1))""",
"events":"""select distinct site, year::smallint as y,month::tinyint as m,day::tinyint as d,specieskey as sp from occurrences
 join repeat_places using(datasetkey,basisofrecord,decimallatitude,decimallongitude) where not(month=1 and day=1)""",
"repeat_taxa":"""select site,sp from events group by all having min(y)<=1999 and max(y)>=2020""",
"visits":"""select site,y,m,d,count(*)::integer as richness from events group by all having count(*)>=5""",
"strata":"""select site,m from visits group by all having count(*) filter (where y<=1999)>=2 and count(*) filter (where y>=2020)>=2""",
"eligible_visits":"""select v.* from visits v semi join strata using(site,m)""",
"hits":"""select e.site,e.sp,e.m,(e.y>=2020)::int as late,count(*)::integer as n from events e semi join eligible_visits v using (site,y,m,d)
 semi join repeat_taxa using(site,sp) group by all""",
"repeated":"""select site,sp from hits group by all having count(distinct late)=2""",
"denominators":"""select site,m,(y>=2020)::int as late,count(*) as visits, avg(richness) as richness from eligible_visits group by all""",
"pairs":"""with rates as (select h.site,h.sp,h.m,h.late,100.0*h.n/d.visits as rate from hits h join denominators d using(site,m,late)
 semi join repeated using(site,sp)), n as (select site,count(*) as nmonth from strata group by site)
 select r.site,r.sp,sum(case when late=0 then rate else 0 end)/nmonth as early,
 sum(case when late=1 then rate else 0 end)/nmonth as recent from rates r join n using(site) group by r.site,r.sp,nmonth""",
}
for key,sql in steps.items():
    if db.execute("select count(*) from information_schema.tables where table_name=?",[key]).fetchone()[0]:continue
    db.execute("create table "+key+" as "+sql)
    print(key,db.execute("select count(*) from "+key).fetchone()[0],flush=True)
taxa=db.execute("""select a.sp,n.species,a.raw_sites,coalesce(b.sites,0) as survey_sites
 from (select sp,count(*) as raw_sites from repeat_taxa group by sp) a join names n on n.specieskey=a.sp
 left join (select sp,count(*) as sites from repeated group by sp) b using(sp)""").df()
taxa.to_csv(out/"repeated_taxa.csv",index=False)
totals=dict(snapshot_objects=len(files),snapshot_rows=int(counts[:,0].sum()),vascular_rows=int(counts[:,1].sum()),eligible_rows=int(counts[:,2].sum()))
for k in steps:totals[k]=int(db.execute("select count(*) from "+k).fetchone()[0])
totals["species_in_raw_repeats"]=len(taxa)
totals["species_in_survey_repeats"]=int(taxa.survey_sites.gt(0).sum())
totals["distinct_sampled_coordinates"]=db.execute("select count(*) from (select distinct decimallatitude,decimallongitude from repeat_places semi join repeated using(site))").fetchone()[0]
(out/"coverage.json").write_text(json.dumps(totals,indent=2))
db.execute("checkpoint")
