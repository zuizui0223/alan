"""All-species annual recording trends; environmental predictor ranking in held-out blocks."""
from pathlib import Path
import sys
import duckdb
import numpy as np
import pandas as pd
import rasterio
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from alan_long_term_environment_forest import environmental_months

CHANGES=["light","temperature","precipitation","built"]
BASE=["light0","temp_mean","rain_mean","built0","effort","log_days_trend","start_year","span","years","island",
      "spring_shift","summer_shift","autumn_shift"]


def annual_data(db,work,green=False):
    if green:
        points=db.sql("select distinct unit_id,longitude,latitude from geo where latitude>=30").df()
        yy=[]
        for path in sorted((work/"gbif_macro/pheno").glob("GIMMS3g_v1_SOS_*_Method_4_30_90.tif")):
            year=int(path.name.split("_")[3])
            with rasterio.open(path) as ds:
                rr,cc=rasterio.transform.rowcol(ds.transform,points.longitude,points.latitude)
                inside=(rr>=0)&(cc>=0)&(rr<ds.height)&(cc<ds.width)
                vals=np.full(len(points),np.nan)
                vals[inside]=ds.read(1)[rr[inside],cc[inside]]
            ok=(vals>=1)&(vals<=366)
            yy.append(points.loc[ok,["unit_id"]].assign(y=year,green=vals[ok]))
        onset=pd.concat(yy)
        assert set(onset.y)==set(range(1992,2015)) and not onset.duplicated(["unit_id","y"]).any()
        db.register("onset",onset)
        db.execute("create temp view panel_y as select * from year_surveys join onset using(unit_id,y)")
    else:
        db.execute("create temp view panel_y as select * from year_surveys")
    db.execute("""create temp table stations5 as select site from panel_y group by site
                  having count(distinct y)>=5 and max(y)-min(y)>=20""")
    xy=db.sql("""select distinct unit_id,longitude,latitude from geo semi join stations5 using(site)""").df()
    yr=db.sql("""select distinct unit_id,y from panel_y semi join stations5 using(site)""").df()
    dates=yr.merge(xy,on="unit_id",validate="m:1").merge(pd.DataFrame({"m":range(1,13)}),how="cross")
    e=environmental_months(dates,work/"gbif_annual/ncep",2015 if green else 2020)
    env=e.groupby(["unit_id","y"],as_index=False).agg(air=("air","mean"),prate=("prate","mean"),built0=("built0","first"),built=("built_change","first"))
    db.register("environment",env)
    # Shared survey effort is retained even for zero-detection years after a species' last record.
    db.execute("""create temp table source_env as select site,unit_id,block,land_class,
      10*regr_slope(dn,y) light, arg_min(dn,y) light0, 10*regr_slope(air,y) temperature,
      864000*regr_slope(prate,y) precipitation,avg(air) temp_mean,86400*avg(prate) rain_mean,
      avg(built0) built0,avg(built) built,10*regr_slope(ln(richness),y) effort,
      10*regr_slope(ln(visits),y) log_days_trend,min(y) start_year,max(y)-min(y) span,
      count(*) as years,10*regr_slope(month_fraction_3+month_fraction_4+month_fraction_5,y) spring_shift,
      10*regr_slope(month_fraction_6+month_fraction_7+month_fraction_8,y) summer_shift,
      10*regr_slope(month_fraction_9+month_fraction_10+month_fraction_11,y) autumn_shift
      """+(",10*regr_slope(green,y) as greenup,avg(green) as green_mean " if green else "")+
      "from panel_y join stations5 using(site) join environment using(unit_id,y) group by all")
    extra=["greenup","green_mean"] if green else []
    cols=", ".join(f"avg({c}) as {c}" for c in CHANGES+BASE+extra if c!="island")
    db.execute("""create temp table at_risk as select distinct site,sp from year_hits
                  semi join panel_y using(site,y) semi join stations5 using(site)""")
    return db.sql("""with response as (select sp,s.site,10*regr_slope(100.*coalesce(hits,0)/visits,y) as change
        from panel_y s join stations5 using(site) join at_risk using(site)
        left join year_hits using(site,y,sp) group by all)
        select sp,unit_id,block,(land_class='island')::int as island,avg(change) as change,"""+cols+
        " from response join source_env using(site) group by all").df()


def evaluate(sp,g,confirm=False):
    g=g.sort_values("unit_id")
    row=dict(sp=sp,trajectory_sites=len(g))
    changes=CHANGES+(["greenup"] if "greenup" in g else [])
    cols=BASE+(["green_mean"] if "greenup" in g else [])+changes
    g=g.dropna(subset=cols+["change"])
    groups=g.block.to_numpy()
    n=len(g)
    ng=g.block.nunique()
    row.update(sites=n,blocks=ng)
    if ng<2:
        return dict(row,status="insufficient_spatial_blocks")
    X=g[cols].to_numpy(dtype=float)
    y=g.change.to_numpy()
    k=min(5,ng)
    split=GroupKFold(n_splits=k,shuffle=confirm,random_state=391 if confirm else None)
    folds=list(split.split(X,y,groups))
    if any(len(tr)<2 for tr,_ in folds) or np.ptp(y)<1e-9:
        return dict(row,status="insufficient_train_or_response_variation")
    _,inverse,counts=np.unique(groups,return_inverse=True,return_counts=True)
    w=1/counts[inverse]
    mse=np.zeros((k,2+len(changes)))
    local_diff=0.
    local_mass=0.
    for j,(tr,te) in enumerate(folds):
        mse[j,0]=np.dot(w[te],(y[te]-np.mean(y[tr]))**2)
        for h,drop in enumerate([None,*changes],1):
            keep=np.array([c!=drop for c in cols])
            model=RandomForestRegressor(n_estimators=64 if confirm else 24,
                        min_samples_leaf=min(5,max(1,len(tr)//12)),random_state=(1909 if confirm else 293)+j)
            model.fit(X[tr][:,keep],y[tr])
            err=y[te]-model.predict(X[te][:,keep])
            mse[j,h]=np.dot(w[te],err**2)
            if not drop:
                i=cols.index("light")
                edges=np.unique(np.quantile(X[tr,i],np.linspace(0,1,6)))
                for a,b in zip(edges[:-1],edges[1:]):
                    ids=te[(X[te,i]>=a)&((X[te,i]<b)|((b==edges[-1])&(X[te,i]==b)))]
                    if not len(ids):continue
                    left=X[ids].copy()
                    right=left.copy()
                    left[:,i]=a
                    right[:,i]=b
                    effects=(model.predict(right)-model.predict(left))/(b-a)
                    local_diff+=np.dot(w[ids],effects)
                    local_mass+=w[ids].sum()
    total=mse.sum(axis=0)
    if total[0]<=1e-10:return dict(row,status="constant_test_baseline")
    gains=(total[2:]-total[1])/total[0]
    row.update(zip(["gain_"+v for v in changes],gains))
    row["top"]=changes[gains.argmax()] if (gains==gains.max()).sum()==1 else "tie"
    row.update(status="evaluated",folds=k,skill=1-total[1]/total[0],
               light_winning_folds=int((mse[:,2]>mse[:,1]).sum()),
               light_response=local_diff/local_mass if local_mass else np.nan,
               observed_at_brightening=g.loc[g.light>1.e-8,"change"].mean())
    return row


def main(work,out,green=False):
    out.mkdir(exist_ok=True,parents=True)
    db=duckdb.connect(str(work/"gbif_annual/annual.duckdb"),read_only=True)
    taxa=db.sql("from species_names").df()
    data=annual_data(db,work,green)
    # Cache inputs only in the acquisition workspace.
    data.to_parquet(work/"gbif_annual"/("green_rf_sites.parquet" if green else "all_species_rf_sites.parquet"),index=False)
    for confirm in [False,True]:
        if confirm:
            ids=ledger.loc[ledger.status.eq("evaluated") & ledger.top.eq("light") &
                           ledger.skill.gt(0) & ledger.gain_light.gt(0),"sp"]
            cohort=data[data.sp.isin(ids)]
        else:cohort=data
        fits=Parallel(n_jobs=4 if green else 8,verbose=1)(delayed(evaluate)(sp,g,confirm) for sp,g in cohort.groupby("sp"))
        res=pd.DataFrame(fits)
        if not confirm:
            ledger=taxa.merge(res,how="left",validate="1:1")
            ledger["status"]=ledger.status.fillna("no_qualifying_long_series")
            ledger.to_csv(out/"all_species.csv.gz",index=False)
        else:
            keep=["sp","gain_light","light_response","skill","top","light_winning_folds","folds"]
            r=ledger[ledger.sp.isin(ids)].merge(res[keep],on="sp",suffixes=("","_repeat"),validate="1:1")
            r["stable_light"]=(r.top_repeat.eq("light") & r.skill_repeat.gt(0) & r.gain_light_repeat.gt(0)
                & r.light_winning_folds.ge(4) & r.light_winning_folds_repeat.ge(4) & r.folds.eq(5)
                & (r.light_response*r.light_response_repeat).gt(0))
            r["direction"]=np.select([r.light_response.gt(0),r.light_response.lt(0)],
                                     ["positive","negative"],default="undetermined")
            r["net_agrees"]=(r.observed_at_brightening*r.light_response).gt(0)
            r.sort_values("gain_light_repeat",ascending=False).to_csv(out/"repartition.csv",index=False)
            r[r.stable_light].sort_values("gain_light_repeat",ascending=False).to_csv(out/"candidates.csv",index=False)
            assert len(ledger)==len(taxa) and set(r.sp)==set(ids)
            print(ledger.status.value_counts(),r[r.stable_light].groupby(["direction","net_agrees"]).size())


if __name__=="__main__":
    main(*map(Path,sys.argv[1:3]),green="--green" in sys.argv)
