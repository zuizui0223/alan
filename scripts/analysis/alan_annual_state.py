"""Relate annual wetland state to recording rates of the two fixed wet-plant groups."""
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
import duckdb
import xarray as xr
from scipy.stats import t
from statsmodels.stats.multitest import multipletests
from alan_wet_shore import repo,cache

out=repo/"results/all_species_forest/wet_annual"
g=pd.read_csv(out/"habitats.csv").merge(pd.read_csv(out/"coastal_screen.csv"),validate="1:1")
g["inland"]=g.coastal_classes_in_2km_square.isna()
annual=pd.read_csv(out/"water.csv").merge(g,on="unit_id",validate="m:1")
assert annual.status.eq("ok").all()
for y,idx in annual.groupby("y").groups.items():
    path=next((repo.parent/"phenology_run/decades").glob(f"Harmonized_DN_NTL_{y}_*.tif"))
    with rasterio.open(path) as im:
        assert im.crs.to_epsg()==4326
        annual.loc[idx,"light"]=[float(v[0]) for v in im.sample(annual.loc[idx,["longitude","latitude"]].to_numpy())]
for var in ["air","prate"]:
    with xr.open_dataset(next((cache/"ncep").glob(var+"*.nc"))) as ds:
        a=ds[var].sel(time=slice("1999","2023")).groupby("time.year").mean().load()
        a=xr.concat([a,a.isel(lon=[0]).assign_coords(lon=[360.])],dim="lon")
        annual[var]=a.sel(year=xr.DataArray(annual.y.values,dims="point"),lat=xr.DataArray(annual.latitude.values,dims="point"),lon=xr.DataArray(annual.longitude.values%360,dims="point"),method="nearest").values
annual["rain"]=annual.prate*86400
veg=pd.read_csv(out/"vegetation.csv.gz").dropna(subset=["ndvi"])
veg["date"]=pd.to_datetime(veg.y.astype(str)+"-01-01")+pd.to_timedelta(veg.day-1,unit="D")
veg=veg[veg.date.dt.month.between(4,9)].copy()
veg["season"]=((veg.date.dt.month-4)//2).astype(int)
seasons=veg.groupby(["unit_id","y","season"]).ndvi.median().unstack()
ndvi=seasons.mean(axis=1).where(seasons.notna().sum(axis=1).eq(3))
annual=annual.merge(ndvi.rename("ndvi"),on=["unit_id","y"],how="left",validate="1:1")
monthly=pd.read_csv(out/"monthly_water.csv.gz")
monthly["quarter"]=(monthly.m-1)//3
def summarize(f):
    a=f.dropna(subset=["gain","loss"])
    eligible=len(a)>=6 and a.quarter.nunique()==4
    v=a.groupby("quarter")[["gain","loss"]].mean().mean()
    return pd.Series(dict(month_pairs=len(a),quarters=a.quarter.nunique(),rewet=v.gain if eligible else np.nan,dry=v.loss if eligible else np.nan))
trans=monthly.groupby(["unit_id","y"]).apply(summarize,include_groups=False)
annual=annual.merge(trans,on=["unit_id","y"],validate="1:1").sort_values(["unit_id","y"])
for var in ["light","water_freq","water_area","ndvi"]:
    annual["change_"+var]=annual.groupby("unit_id")[var].diff()
    annual["abs_"+var]=annual["change_"+var].abs()
annual["turnover"]=annual.rewet+annual.dry
summer=monthly[monthly.m.between(6,8)].groupby(["unit_id","y"])[["gain","loss"]]
summer_means=summer.mean().where(summer.count()>=2).rename(columns={"gain":"summer_rewet","loss":"summer_dry"})
annual=annual.merge(summer_means,on=["unit_id","y"],validate="1:1")
annual.to_csv(out/"annual_environment.csv.gz",index=False)
survey=pd.read_csv(out/"survey_series.csv.gz").query("ever_both")
assert not survey.duplicated(["unit_id","y"]).any()
panel=survey.merge(annual,on=["unit_id","y"],validate="1:1").query("inland").copy()
np.testing.assert_allclose(panel.dn,panel.light)
panel["difference"]=panel.candidate_rate-panel.reference_rate
panel["logvisits"]=np.log(panel.visits)
panel.to_csv(out/"annual_survey.csv.gz",index=False)
def slope(x,y):
    ok=np.isfinite(x)&np.isfinite(y)
    return 10*np.polyfit(x[ok],y[ok],1)[0] if ok.sum()>=5 and np.ptp(x[ok])>=15 else np.nan
rows=[]
for site,d in panel.groupby("unit_id"):
    f=annual[annual.unit_id.eq(site)]
    row=dict(unit_id=site,habitat=d.habitat.iloc[0],block=d.block.iloc[0],survey_years=len(d),light_initial=d.sort_values("y").dn.iloc[0])
    for col in ["candidate_rate","reference_rate","difference"]:
        row[col]=slope(d.y,d[col])
    row["light"]=slope(f.y,f.light)
    row["survey_light"]=slope(d.y,d.light)
    for col in ["water_freq","ndvi"]:
        row["drift_"+col]=abs(slope(f.y,f[col]))
    for col in ["turnover","rewet","dry","abs_water_freq","abs_ndvi"]:
        good=f.dropna(subset=[col])
        row[col]=good[col].mean() if len(good)>=5 and good.y.max()-good.y.min()>=15 else np.nan
    row["waterfreq_mean"]=f.water_freq.mean()
    rows.append(row)
sites=pd.DataFrame(rows)
sites.to_csv(out/"site_trends.csv",index=False)
def model(frame,test,controls,scope,scale,fe=()):
    f=frame.dropna(subset=[test,*controls,*fe,"difference"]).copy()
    X=pd.DataFrame({"intercept":np.ones(len(f)),"target":f[test]/scale},index=f.index)
    X=pd.concat([X,f[controls],*[pd.get_dummies(f[k],prefix=k,drop_first=True,dtype=float) for k in fe]],axis=1)
    x=X.to_numpy(float);ys=f[["candidate_rate","reference_rate","difference"]].to_numpy()
    if len(x)<4 or f.block.nunique()<3:return []
    rank=np.linalg.matrix_rank(x)
    if len(f)<=rank+1 or np.linalg.matrix_rank(np.delete(x,1,axis=1))==rank:return []
    pinv=np.linalg.pinv(x);beta=pinv@ys
    score=pd.DataFrame(pinv[1,:,None]*(ys-x@beta)).groupby(f.block.to_numpy()).sum()
    k=len(score);se=np.sqrt((score**2).sum()*k/(k-1)*(len(f)-1)/(len(f)-rank)).to_numpy()
    half=t.ppf(.975,k-1)*se
    loo=[]
    for block in f.block.unique():
        mask=f.block.ne(block)
        a=x[mask]
        if np.linalg.matrix_rank(np.delete(a,1,axis=1))<np.linalg.matrix_rank(a):loo.append((np.linalg.pinv(a)@ys[mask])[1])
    return [dict(scope=scope,state=test,units=scale,outcome=outcome,rows=len(f),sites=f.unit_id.nunique(),blocks=k,estimate=beta[1,j],low=beta[1,j]-half[j],high=beta[1,j]+half[j],p=2*t.sf(abs(beta[1,j]/se[j]),k-1),loo_min=np.min(loo,axis=0)[j],loo_max=np.max(loo,axis=0)[j]) for j,outcome in enumerate(["candidate_rate","reference_rate","difference"])]
res=[]
sites["bog"]=(sites.habitat==412).astype(int)
for state,scale in [("drift_water_freq",1),("drift_ndvi",.01),("abs_water_freq",1),("turnover",1)]:
    res+=model(sites,state,["light","bog"],"long_term",scale)
# Prior year completed before any botanical survey in the outcome year.
for state,scale in [("abs_water_freq",1),("turnover",1),("rewet",1),("dry",1),("summer_dry",1),("summer_rewet",1),("abs_ndvi",.1)]:
    previous=annual[["unit_id","y",state]].copy();previous.y+=1
    f=panel.drop(columns=state).merge(previous,on=["unit_id","y"],validate="1:1")
    res+=model(f,state,["light","change_light","air","rain","logvisits","summer"],"previous_year",scale,fe=["unit_id","y"])
# Restrict to consecutive plant survey years and estimate actual year-on-year changes.
pairs=panel.sort_values(["unit_id","y"]).copy()
group=pairs.groupby("unit_id")
one_year=group.y.diff().eq(1)
for col in ["candidate_rate","reference_rate","difference","logvisits","summer","air","rain"]:
    pairs[col]=group[col].diff()
pairs=pairs[one_year].copy()
for state,scale in [("abs_water_freq",1),("abs_ndvi",.1),("turnover",1)]:
    res+=model(pairs,state,["change_light","logvisits","summer","air","rain"],"consecutive_surveys",scale)
pairs.to_csv(out/"consecutive_surveys.csv.gz",index=False)
res=pd.DataFrame(res)
for _,idx in res.groupby(["scope","outcome"]).groups.items():res.loc[idx,"q"]=multipletests(res.loc[idx,"p"],method="fdr_bh")[1]
res.to_csv(out/"contrasts.csv",index=False)
print(res[res.outcome.eq("difference")].round(3).to_string(index=False))
