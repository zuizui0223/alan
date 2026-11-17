"""Paired plant recording trends by baseline proximity to open water."""
from alan_wet_shore import cache, out, frozen, repo
import duckdb
import pandas as pd
import numpy as np
from scipy.stats import t

g=frozen("coordinate_trends.csv.gz").merge(pd.read_csv(out/"distance.csv").dropna(subset="distance_m"),on="unit_id",validate="1:1")
g["near"]=1-g.distance_m/1000
g["light10"]=g.light/10
g["near_light"]=g.near*g.light10
db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
wet=pd.read_csv(repo/"results/all_species_forest/wet_partner_composition/species.csv")
water=pd.read_parquet(cache/"glad_water/survey_water.parquet").query("status=='ok'")
db.execute("""create temp table sampled as select s.* from year_surveys s join water using(unit_id,y)
    semi join g using(unit_id) qualify count(*) over(partition by site)>=5 and
    max(y) over(partition by site)-min(y) over(partition by site)>=20""")
db.execute("""create temp table roster as select distinct s.site,s.unit_id,h.sp,w.candidate
 from sampled s join year_hits h using(site,y) join wet w using(sp)
 qualify bool_or(candidate) over(partition by s.site) and bool_or(not candidate) over(partition by s.site)""")
coverage=db.sql("select sp,candidate,count(distinct unit_id) as coordinates from roster group by all").df()
coverage.to_csv(out/"species_coverage.csv",index=False)
wells=pd.read_csv(out/"groundwater_years.csv.gz")
wg=db.sql("""with dates as (select distinct unit_id,y from sampled semi join roster using(site))
 select unit_id,first(station) as station,count(*) as matched_years,min(y) as first_y,max(y) as last_y,
  10*regr_slope(level,y) as level_trend from wells join dates using(unit_id,y) group by unit_id""").df()
wg.to_csv(out/"well_support.csv",index=False)
full=g.merge(wg.query("matched_years>=5 and last_y-first_y>=20"),on="unit_id")
env=["light0","temp0","rain0","temp_trend","rain_trend","built0","built","effort","days","spring","summer","autumn","start","span","years","island","nc","nr"]
hydro=["area0","freq0","area_trend","freq_trend"]
rows=[]
def fit(frame,label,adjust):
    predictors=["near_light","light10","near",*env,*adjust]
    outcomes=["candidate_rate","reference_rate","rate_difference"]
    for_slope=frame.groupby("block").filter(lambda z:len(z)>1 and z.light.max()-z.light.min()>1e-8)
    row=dict(model=label,available=len(frame),used=len(for_slope),blocks=for_slope.block.nunique(),brightened=(for_slope.light>1e-8).sum())
    if not len(for_slope):return [row]
    dm=for_slope[predictors+outcomes]-for_slope.groupby("block")[predictors+outcomes].transform("mean")
    x=dm[predictors].to_numpy();y=dm[outcomes].to_numpy()
    rank=np.linalg.matrix_rank(x)
    df=len(x)-row["blocks"]-rank
    if df<1 or row["blocks"]<3 or np.linalg.matrix_rank(x[:,1:])==rank:return [row]
    inv=np.linalg.pinv(x); beta=inv@y
    resid=y-x@beta
    results=[]
    for term,weights in [("near_minus_far",[1,0]),("far_light",[0,1]),("near_light",[1,1])]:
        a=np.zeros(len(predictors));a[:2]=weights
        if not np.allclose(a,a@inv@x):continue
        v=a@inv;b=a@beta
        scores=pd.DataFrame(v[:,None]*resid).groupby(for_slope.block.values).sum()
        share=pd.Series(v**2).groupby(for_slope.block.values).sum();share/=share.sum()
        k=(share>1e-9).sum()
        if k<2:continue
        half=t.ppf(.975,k-1)*np.sqrt((scores**2).sum().values*k/(k-1)*(len(x)-1)/df)
        for i,o in enumerate(outcomes):
            results.append(dict(**row,term=term,outcome=o,estimate=b[i],low=b[i]-half[i],high=b[i]+half[i],effective_blocks=1/(share@share)))
    return results
for label,frame,adj in [("base",g,[]),("surface",g,hydro),("surface_brightened",g[g.light>1e-8],hydro),("well_subset_surface",full,hydro),("well_level",full,hydro+["level_trend"])]:
    rows.extend(fit(frame,label,adj))
pd.DataFrame(rows).to_csv(out/"effects.csv",index=False)
g["distance_band"]=pd.cut(g.distance_m,[-1,100,500,999.999999,1000],labels=["within_100m","100_500m","500_1000m","capped_1km"])
g.groupby("distance_band",observed=True).agg(sites=("unit_id","size"),brightened=("light",lambda x:(x>1e-8).sum())).to_csv(out/"distance_support.csv")
print(pd.DataFrame(rows).query("term=='near_minus_far' or term!=term").round(3).to_string(index=False))
print("wells",len(full),full.station.nunique(),full.block.nunique(),full.island.value_counts().to_dict(),(full.light>1e-8).sum())
print(coverage.groupby("candidate").size())
