"""Compare annual trajectories between sites with stronger and weaker long-term brightening."""
from pathlib import Path
import sys
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
from alan_global_gbif_screen import model


def main(cache,out):
    db=duckdb.connect(str(cache/"annual.duckdb"),config={"threads":3,"memory_limit":"1500MB"})
    db.execute("""create or replace table trend_months as select site,m
        from supported group by all having count(distinct y)>=5 and max(y)-min(y)>=20""")
    db.execute("""create or replace table trajectories as
      with slopes as (
       select r.sp,s.site,s.m,s.unit_id,s.block,s.land_class,
         10*regr_slope(dn,y) as light_trend,
         10*regr_slope(100.0*coalesce(h.hits,0)/visits,y) as rate_trend,
         10*regr_slope(ln(richness),y) as effort_trend,arg_min(dn,y) as light0
       from risk r join supported s using(site,m) semi join trend_months using(site,m)
        left join hits h using(site,m,y,sp) group by all),
      stations as (select sp,site,unit_id,block,land_class,avg(light_trend) as l,avg(rate_trend) as r,
        avg(effort_trend) as e,avg(light0) as l0 from slopes group by all)
      select sp,unit_id,block,land_class,avg(l) as delta,avg(r) as change,avg(e) as rich_change,avg(l0) as light0
      from stations group by all""")
    data=db.sql("select * from trajectories").df()
    results=[]
    for sp,g in data.groupby("sp"):
        row=dict(sp=sp,coordinates=len(g),blocks=g.block.nunique(),brightening=int(g.delta.gt(0).sum()),
                 brightening_mean_change_pp_per_decade=g.loc[g.delta.gt(0),"change"].mean())
        fit=model(g)
        if fit:
            b,se,groups=fit[:3]
            ci=t.ppf(.975,groups-1)*se
            row.update(slope=b,lo=b-ci,hi=b+ci,p=2*t.sf(abs(b/se),groups-1))
        results.append(row)
    taxa=db.sql("from species_names").df()
    d=taxa.merge(pd.DataFrame(results),how="left",validate="1:1")
    d["q"]=np.nan
    if "p" in d:
        ix=d.p.notna()
        if ix.any():d.loc[ix,"q"]=multipletests(d.loc[ix,"p"],method="fdr_bh")[1]
    d.to_csv(out/"long_term_screen.csv.gz",index=False,compression="gzip")
    print("site trend candidates",len(d.query("q<.05")),flush=True)


if __name__=="__main__":main(*map(Path,sys.argv[1:3]))
