"""Spatial holdout ablation of four environmental changes for long-term recording-rate slopes."""
import sys
from pathlib import Path
from zipfile import ZipFile

import duckdb
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from joblib import Parallel, delayed
from pyproj import Transformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold


def environmental_months(dates,cache,built_end=2020):
    sites=dates.drop_duplicates("unit_id").set_index("unit_id")
    for var in ["air","prate"]:
        with xr.open_dataset(next(cache.glob(var+"*.nc"))) as dataset:
            arr=dataset[var].sel(time=slice("1992","2024")).load()
            # Nearest grid; cyclic longitude uses the equivalent first grid column at 360 degrees.
            arr=xr.concat([arr,arr.isel(lon=[0]).assign_coords(lon=[360.])],dim="lon")
            sub=arr.sel(lat=xr.DataArray(sites.latitude.values,dims="point"),
                        lon=xr.DataArray(sites.longitude.values%360,dims="point"),method="nearest")
            t=pd.DatetimeIndex(sub.time.values)
            lookup=pd.DataFrame(sub.transpose("point","time").values,index=sites.index,
                                columns=t.year*100+t.month)
            month_index=lookup.columns.get_indexer(dates.y.astype("int64")*100+dates.m)
            point_index=lookup.index.get_indexer(dates.unit_id)
            assert min(month_index)>=0 and min(point_index)>=0 and len(set(month_index))>200
            dates[var]=lookup.to_numpy()[point_index,month_index]
    for yr in [1990,built_end]:
        path=next(cache.glob(f"*E{yr}*.zip"))
        with ZipFile(path) as z:
            member=next(n for n in z.namelist() if n.lower().endswith(".tif"))
            # Extract once to avoid repeated decompression for scattered point queries.
            target=cache/Path(member).name
            if not target.exists():
                import shutil
                with z.open(member) as inp,target.open("wb") as out:
                    shutil.copyfileobj(inp,out)
        with rasterio.open(target) as src:
            x,y=Transformer.from_crs(4326,src.crs,always_xy=True).transform(sites.longitude,sites.latitude)
            sites[f"built{yr}"]=[float(v[0]) if not np.ma.is_masked(v[0]) else np.nan for v in src.sample(zip(x,y),masked=True)]
    sites["built_change"]=(sites[f"built{built_end}"]-sites.built1990)/((built_end-1990)/10)/1e4
    dates["built0"]=dates.unit_id.map(sites.built1990)/1e4
    dates["built_change"]=dates.unit_id.map(sites.built_change)
    return dates


def exposures(db,cache):
    dates=db.sql("""select s.site,s.m,s.y,s.unit_id,g.longitude,g.latitude
          from supported s join trend_months using(site,m) join geo g on g.site=s.site""").df()
    dates=environmental_months(dates,cache)
    db.register("months_env",dates)
    return db.sql("""with mm as (select site,m,unit_id,10*regr_slope(air,y) as temp,
         864000*regr_slope(prate,y) as rain, avg(air) as temp0,86400*avg(prate) as rain0,
         min(y) as start, max(y)-min(y) as span,count(*) as years,
         avg(built0) as built0,avg(built_change) as built_change from months_env group by all),
         per_source as (select sp,site,unit_id,avg(temp) as temp,avg(rain) rain,avg(temp0) temp0,
           avg(rain0) rain0,avg(start) as first_year,avg(span) as span,avg(years) as years,
           avg(built0) built0,avg(built_change) built_change from mm join risk using(site,m) group by all),
         per_point as (select sp,unit_id,avg(temp) as temp, avg(rain) rain, avg(temp0) temp0, avg(rain0) rain0,
           avg(first_year) first_year,avg(span) span,avg(years) as years,avg(built0) built0,
           avg(built_change) built_change from per_source group by all)
         select t.*,p.* exclude(sp,unit_id) from trajectories t join per_point p using(sp,unit_id)""").df()


def screen(sp,g,repeat=False):
    changes=["delta","temp","rain","built_change"]
    base=["light0","temp0","rain0","built0","rich_change","first_year","span","years","island"]
    g=g.copy().sort_values("unit_id")
    g["island"]=g.land_class.eq("island").astype(int)
    g=g.dropna(subset=base+changes+["change"])
    row=dict(sp=sp,sites=len(g),blocks=g.block.nunique())
    if row["blocks"]<10 or len(g)<20:
        return row
    cols=base+changes
    X=g[cols].to_numpy()
    y=g.change.to_numpy()
    groups=g.block.to_numpy()
    # Each held-out block contributes equally when evaluating all folds.
    _,idx,counts=np.unique(groups,return_inverse=True,return_counts=True)
    weights=1/counts[idx]
    errors={k:np.zeros(len(g)) for k in ["null","full",*changes]}
    folds=np.zeros(len(g),int)
    split=GroupKFold(n_splits=5,shuffle=True,random_state=73) if repeat else GroupKFold(n_splits=5)
    changes_in_prediction=[]
    for k,(train,test) in enumerate(split.split(X,y,groups)):
        folds[test]=k
        errors["null"][test]=(y[test]-y[train].mean())**2
        for drop in ["full",*changes]:
            mask=np.array([c!=drop for c in cols])
            rf=RandomForestRegressor(n_estimators=96 if repeat else 48,min_samples_leaf=5,max_features=1.,random_state=(710 if repeat else 170)+k)
            rf.fit(X[train][:,mask],y[train])
            errors[drop][test]=(y[test]-rf.predict(X[test][:,mask]))**2
            if repeat and drop=="full":
                # Local accumulated-effect increments. Only evaluate test sites inside training support.
                j=cols.index("delta")
                edges=np.unique(np.quantile(X[train,j],np.linspace(0,1,6)))
                for a,b in zip(edges[:-1],edges[1:]):
                    within=test[(X[test,j]>=a)&((X[test,j]<b)|((b==edges[-1])&(X[test,j]==b)))]
                    if len(within):
                        left=X[within].copy()
                        right=left.copy()
                        left[:,j]=a
                        right[:,j]=b
                        changes_in_prediction.extend(zip((rf.predict(right)-rf.predict(left))/(b-a),weights[within]))
    mse={key:np.average(val,weights=weights) for key,val in errors.items()}
    row["skill"]=1-mse["full"]/mse["null"]
    for v in changes:
        row["gain_"+v]=(mse[v]-mse["full"])/mse["null"]
    row["light_winning_folds"]=sum(np.average((errors["delta"]-errors["full"])[folds==k],weights=weights[folds==k])>0 for k in range(5))
    row["top"]=max(changes,key=lambda v:row["gain_"+v])
    if changes_in_prediction:
        local,weight=zip(*changes_in_prediction)
        row["local_rf_direction"]=np.average(local,weights=weight)
    return row


def main(work,results,repeat_only=False):
    parent=results/"long_term_land_contrasts"
    out=parent/"environment_forest"
    out.mkdir(exist_ok=True)
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    env=exposures(db,work/"gbif_annual/ncep")
    known=pd.read_csv(parent/"all_species.csv.gz",dtype={"sp":str}).dropna(subset=["slope"])
    if repeat_only:
        first=pd.read_csv(out/"species.csv",dtype={"sp":str})
        known=known[known.sp.isin(first.loc[first.light_candidate,"sp"])]
    data=env[env.sp.isin(known.sp)]
    rows=Parallel(n_jobs=6,verbose=1)(delayed(screen)(sp,g,repeat_only) for sp,g in data.groupby("sp"))
    res=known[["sp","species","slope","q","brightening_mean_change_pp_per_decade"]].merge(pd.DataFrame(rows),validate="1:1")
    assert len(res)==len(known)
    res["light_candidate"]=(res.top.eq("delta") & res.skill.gt(0) & res.gain_delta.gt(0) &
                            res.light_winning_folds.ge(4) & res.slope.lt(0) & res.brightening_mean_change_pp_per_decade.lt(0))
    res.sort_values("gain_delta",ascending=False).to_csv(out/("repartition.csv" if repeat_only else "species.csv"),index=False)


if __name__=="__main__":
    main(*map(Path,sys.argv[1:3]),repeat_only="--repeat" in sys.argv)
