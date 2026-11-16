"""Compare within-site, within-month annual recording rates across local light conditions."""
import sys
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.special import betainc

def run(cache,out):
    out.mkdir(parents=True,exist_ok=True)
    db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True,config={"threads":3})
    geo=db.sql("select * from geo").df()
    x,y=Transformer.from_crs(4326,6933,always_xy=True).transform(geo.longitude,geo.latitude)
    geo["cell"]=(np.floor(x/10000).astype(int)*10000+np.floor(y/10000).astype(int))
    mainland="M_"+np.floor((geo.longitude+180)/5).astype(int).astype(str)+"_"+np.floor((geo.latitude+90)/5).astype(int).astype(str)
    geo["region"]=np.where(geo.land_class.eq("island"),"I_"+geo.land_id.fillna("").astype(str),mainland)
    db.register("regions",geo[["site","unit_id","region","cell","land_class"]])
    query="""
      with strata as (select site,m,avg(dn) as mid from eligible group by all having max(dn)-min(dn)>=1),
      obs as (select r.sp,s.site,s.m,s.y,s.visits,s.dn,s.dn>k.mid as bright,
                     100.0*coalesce(h.hits,0)/s.visits as rate
          from risk r join eligible s using(site,m) join strata k using(site,m)
          left join hits h using(site,m,y,sp)),
      contrasts as (select sp,site,m,avg(rate) filter(where bright)-avg(rate) filter(where not bright) as delta,
          avg(dn) filter(where bright)-avg(dn) filter(where not bright) as delta_dn,
          count(*) as survey_month_years from obs group by all),
      sources as (select sp,site,avg(delta) delta,avg(delta_dn) delta_dn,sum(survey_month_years) as effort from contrasts group by all),
      points as (select sp,region,cell,unit_id,land_class,avg(delta) delta,avg(delta_dn) delta_dn,sum(effort) effort
                 from sources join regions using(site) group by all),
      cells as (select sp,region,cell,land_class,avg(delta) delta,avg(delta_dn) delta_dn,count(*) as coords,sum(effort) effort
                from points group by all)
      select sp,region,land_class,count(*) as blocks,sum(coords) as coordinates,sum(effort) as survey_month_years,
          count(*) filter(where delta < -1e-8) neg,count(*) filter(where delta>1e-8) pos,
          avg(delta) as mean_delta_pp,median(delta) as median_delta_pp,median(delta_dn) as median_light_contrast_dn
      from cells group by all
    """
    d=db.sql(query).df().merge(db.sql("from species_names").df(),validate="m:1")
    d["ties"]=d.blocks-d.neg-d.pos
    d["prob_negative"]=1-betainc(d.neg+2,d.pos+2,.5)
    d["label"]="limited"
    enough=d.blocks.ge(3)
    d.loc[enough,"label"]="mixed_or_weak"
    for key,sign,prob in [("neg",-1,d.prob_negative),("pos",1,1-d.prob_negative)]:
        flag=enough & prob.ge(.9) & d[key].gt(d.blocks/2) & (sign*d.mean_delta_pp>0)
        d.loc[flag,"label"]=key
    d.sort_values(["species","region"]).to_csv(out/"all_regions.csv.gz",index=False)
    by=d.pivot_table(index=["sp","species"],columns="label",values="region",aggfunc="count",fill_value=0)
    by.columns.name=None
    by.to_csv(out/"species.csv")
    print(d.groupby(["land_class","label"]).agg(pairs=("sp","size"),species=("sp","nunique"),regions=("region","nunique")).to_string())

if __name__=="__main__":run(*map(Path,sys.argv[1:3]))
