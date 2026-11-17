"""Compare seven frozen plants with the other insect-pollinated wet plants surveyed at each site."""
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb
import xarray as xr
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

repo=Path(__file__).resolve().parents[2]
root=repo/"results/all_species_forest"
out=root/"seven_moth"
cache=repo.parent/"gbif_annual"
species=pd.read_csv(root/"wet_reproduction/species.csv")
species["focus"]=species.candidate & species.documented_transport
assert species.focus.sum()==7
species=species[species["Pollination syndrome"].str.contains("insect-pollination",na=False)].copy()
assert species.focus.sum()==7 and len(species)==140
db=duckdb.connect(str(cache/"annual.duckdb"),read_only=True)
# Retain a single original survey series per coordinate to keep exposure and survey dates identical.
points=db.sql("""select unit_id,first(latitude) latitude,first(longitude) longitude
  from geo group by unit_id having count(distinct site)=1""").df()
dates=db.sql("select distinct unit_id,y from year_surveys semi join points using(unit_id)").df().merge(points)
for var in ["air","prate"]:
    with xr.open_dataset(next((cache/"ncep").glob(var+"*.nc"))) as ds:
        a=ds[var].sel(time=slice("1992","2024")).groupby("time.year").mean().load()
        a=xr.concat([a,a.isel(lon=[0]).assign_coords(lon=[360.])],dim="lon")
        dates[var]=a.sel(year=xr.DataArray(dates.y.values,dims="point"),
                         lat=xr.DataArray(dates.latitude.values,dims="point"),
                         lon=xr.DataArray(dates.longitude.values%360,dims="point"),method="nearest").values
built=pd.read_parquet(cache/"all_species_rf_sites.parquet",columns=["unit_id","built0","built"]).drop_duplicates()
built=built[built.unit_id.isin(points.unit_id)]
assert (built.groupby("unit_id").max()-built.groupby("unit_id").min()).abs().max().max()<1e-8
built=built.groupby("unit_id",as_index=False).mean()
dates=dates.merge(built,on="unit_id",validate="m:1")
sat=pd.read_parquet(cache/"glad_water/survey_water.parquet").query("status=='ok'")
wells=pd.read_csv(out/"groundwater_years.csv.gz")
controls=["temp","rain","built","temp0","rain0","built0","light0","effort","days","start","span","years","spring","summer","autumn"]
outputs=[]
inter=[]
equal={}
well_audit=[]
detail=[]
def fit(frame,extra):
    g=frame.groupby("block").filter(lambda z:len(z)>1 and z.light.max()-z.light.min()>1e-8).copy()
    row=dict(sites=len(frame),used=len(g),blocks=g.block.nunique(),brightening=int(g.light.gt(0).sum()))
    if g.empty:return row,None
    bright=g[g.light.gt(0)]
    row.update(focal_bright=bright.focal.mean(),reference_bright=bright.reference.mean(),bright_diff=bright.diff_rate.mean())
    if row["blocks"]<3:return row,None
    cols=["diff_rate","light",*controls,*extra]
    X=g[cols]-g.groupby("block")[cols].transform("mean")
    y=X.pop("diff_rate").to_numpy()
    X["light"]/=10
    X=X.loc[:,X.std().gt(1e-9)]
    if "light" not in X:return row,None
    mat=X.to_numpy()
    rank=np.linalg.matrix_rank(mat)
    df=len(g)-row["blocks"]-rank
    if df<1 or np.linalg.matrix_rank(X.drop(columns="light"))==rank:return row,None
    # Cluster influence on the light coefficient, allowing shared blocks across land classes.
    pinv=np.linalg.pinv(mat)
    beta=pinv@y
    influence=pd.Series(pinv[X.columns.get_loc("light")]*(y-mat@beta),index=g.index).groupby(g.block).sum()
    influence*=np.sqrt(row["blocks"]/(row["blocks"]-1)*(len(g)-1)/df)
    b=beta[X.columns.get_loc("light")]
    se=np.linalg.norm(influence)
    half=t.ppf(.975,row["blocks"]-1)*se
    row.update(estimate=b,low=b-half,high=b+half,p=2*t.sf(abs(b/se),row["blocks"]-1))
    return row,influence
for mode in ["long","satellite","groundwater"]:
    env=dates.copy()
    extra=[]
    if mode=="satellite":
        env=env.merge(sat,on=["unit_id","y"],validate="1:1")
        extra=["area0","freq0","area_trend","freq_trend"]
    if mode=="groundwater":
        env=env.merge(wells[["unit_id","y","level"]],on=["unit_id","y"],validate="1:1")
        extra=["level_trend"]
    db.execute("""create or replace temp table sampled as select s.*,e.* exclude(unit_id,y) from year_surveys s
        join env e using(unit_id,y) qualify count(*) over(partition by site)>=5 and
        max(y) over(partition by site)-min(y) over(partition by site)>=20""")
    if db.sql("select count(*) from sampled").fetchone()[0]==0:continue
    db.execute("""create or replace temp table roster as select distinct h.site,sp,species,focus from year_hits h
         join species using(sp) semi join sampled using(site,y)""")
    db.execute("""create or replace temp table rates as select site,sp,species,focus,
          10*regr_slope(100.*coalesce(hits,0)/visits,y) rate from sampled join roster using(site)
          left join year_hits using(site,sp,y) group by all""")
    hydro=",arg_min(water_area,y) area0,arg_min(water_freq,y) freq0,10*regr_slope(water_area,y) area_trend,10*regr_slope(water_freq,y) freq_trend" if mode=="satellite" else ",10*regr_slope(level,y) level_trend" if mode=="groundwater" else ""
    meta=db.sql("""select site,unit_id,block,land_class,10*regr_slope(dn,y) light,arg_min(dn,y) light0,
       10*regr_slope(air,y) as temp,864000*regr_slope(prate,y) rain,avg(air) temp0,86400*avg(prate) rain0,
       avg(built) built,avg(built0) built0,10*regr_slope(ln(richness),y) effort,10*regr_slope(ln(visits),y) as days,
       min(y) as start,max(y)-min(y) span,count(*) as years,
       10*regr_slope(month_fraction_3+month_fraction_4+month_fraction_5,y) spring,
       10*regr_slope(month_fraction_6+month_fraction_7+month_fraction_8,y) summer,
       10*regr_slope(month_fraction_9+month_fraction_10+month_fraction_11,y) autumn"""+hydro+" from sampled group by all").df()
    r=db.table("rates").df()
    ref=r[~r.focus].groupby("site").rate.mean().rename("reference")
    focal=r[r.focus]
    for label in ["seven_site_mean",*species.loc[species.focus,"species"]]:
        f=focal if label=="seven_site_mean" else focal[focal.species.eq(label)]
        frame=f.groupby("site").rate.mean().rename("focal").to_frame().join(ref,how="inner").join(meta.set_index("site")).reset_index(drop=True)
        frame["diff_rate"]=frame.focal-frame.reference
        if label=="seven_site_mean":
            detail.append(frame.assign(panel=mode))
        for adjustment,variables in ([("base",[]),("water",extra)] if extra else [("base",[])]):
            effects={}
            for land,g in [("all",frame),("island",frame[frame.land_class.eq("island")]),("mainland",frame[frame.land_class.eq("mainland")])]:
                row,inf=fit(g,variables)
                outputs.append(dict(panel=mode,adjustment=adjustment,plant=label,land=land,**row))
                if inf is not None:
                    effects[land]=(row,inf)
                    if label!="seven_site_mean":
                        equal.setdefault((mode,adjustment,land),[]).append((row,inf))
                if label=="seven_site_mean" and mode=="groundwater" and land=="all" and variables:
                    for block in g.block.unique():
                        loo,_=fit(g[g.block.ne(block)],variables)
                        well_audit.append(dict(omitted_block=block,**loo))
            if "island" in effects and "mainland" in effects:
                a,ia=effects["island"];b,ib=effects["mainland"]
                diff=a["estimate"]-b["estimate"]
                err=np.linalg.norm(ia.subtract(ib,fill_value=0))
                half=t.ppf(.975,min(a["blocks"],b["blocks"])-1)*err
                inter.append(dict(panel=mode,adjustment=adjustment,plant=label,difference=diff,low=diff-half,high=diff+half))
equal_summary=[]
for (mode,adjustment,land),items in equal.items():
    if len(items)!=7:continue
    beta=np.mean([v["estimate"] for v,_ in items])
    scores=pd.concat([x for _,x in items],axis=1).fillna(0).mean(axis=1)
    half=t.ppf(.975,min(v["blocks"] for v,_ in items)-1)*np.linalg.norm(scores)
    equal_summary.append(dict(panel=mode,adjustment=adjustment,land=land,estimate=beta,low=beta-half,high=beta+half))
pd.DataFrame(equal_summary).to_csv(out/"equal_species.csv",index=False)
pd.DataFrame(well_audit).to_csv(out/"well_influence.csv",index=False)
tab=pd.DataFrame(outputs)
for _,g in tab[tab.plant.ne("seven_site_mean") & tab.p.notna()].groupby(["panel","adjustment","land"]):
    tab.loc[g.index,"holm_seven"]=multipletests(list(g.p)+[1]*(7-len(g)),method="holm")[1][:len(g)]
tab.to_csv(out/"site_contrasts.csv",index=False)
pd.DataFrame(inter).to_csv(out/"island_minus_mainland.csv",index=False)
pd.concat(detail,ignore_index=True).to_csv(out/"matched_coordinates.csv.gz",index=False)
