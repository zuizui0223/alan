"""Align Swiss observed flowering onset with annual satellite green-up and night light."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

def demean(array,groups):
    for codes in groups:
        counts=np.bincount(codes)
        for v in array.T:
            v-= (np.bincount(codes,weights=v)/counts)[codes]

def estimate(g):
    if g.station.nunique()<10 or g.block.nunique()<10:return []
    vals=g[["dn","flower","sos","signed_gap","abs_gap"]].to_numpy(float)
    vals[:,0]/=10
    groups=[pd.factorize(g[k])[0] for k in ["station","year"]]
    residual=vals.copy()
    for _ in range(200):
        prev=residual.copy()
        demean(residual,groups)
        if np.max(np.abs(residual-prev))<1e-9:break
    x=residual[:,0]
    ss=x@x
    clusters=pd.factorize(g.block)[0]
    G=len(set(clusters))
    if G<10 or g.station.nunique()<10 or ss<.01:return []
    result=[]
    for j,endpoint in enumerate(["flower","sos","signed_gap","abs_gap"],1):
        y=residual[:,j]
        beta=x@y/ss
        u=y-beta*x
        me=np.bincount(clusters,weights=x*u)
        df=len(g)-g.station.nunique()-g.year.nunique()
        se=np.sqrt(me@me/ss**2 * G/(G-1)*(len(g)-1)/df)
        half=t.ppf(.975,G-1)*se
        result.append(dict(outcome=endpoint,n=len(g),stations=g.station.nunique(),years=g.year.nunique(),blocks=G,
                           days_per_10dn=beta,lo=beta-half,hi=beta+half,p=2*t.sf(abs(beta/se),G-1)))
    return result

def summarize(panel):
    estimates=[]
    for name,g in panel.groupby("species"):
        g=g.groupby("station").filter(lambda h:len(h)>=3)
        for row in estimate(g):estimates.append(dict(species=name,**row))
    e=pd.DataFrame(estimates)
    for _,g in e.groupby("outcome"):
        e.loc[g.index,"q"]=multipletests(g.p,method="fdr_bh")[1]
    return e

def panel_at_stations(cache,lights,description="start of flowering"):
    stations=pd.read_csv(cache/"stations.csv",sep=";",encoding="latin1").rename(columns={
        "station_abbr":"station","station_coordinates_wgs84_lat":"lat","station_coordinates_wgs84_lon":"lon"})
    params=pd.read_csv(cache/"parameters.csv",sep=";",encoding="latin1")
    onset=params[params.parameter_description_en.str.contains(description,case=False,na=False,regex=False)]
    species=dict(zip(onset.parameter_shortname,onset.parameter_description_en.str.split(" / ").str[-1]))
    flowers=[]
    for file in cache.glob("*_y.csv"):
        f=pd.read_csv(file,sep=";",encoding="latin1")
        f=f.rename(columns={"station_abbr":"station"})
        f["year"]=pd.to_datetime(f.reference_timestamp,format="%d.%m.%Y %H:%M").dt.year
        f=f[f.year.between(1992,2014)]
        l=f.melt(id_vars=["station","year"],value_vars=list(species),var_name="param",value_name="value").dropna()
        date=pd.to_datetime(l.value.astype(int).astype(str),format="%Y%m%d")
        l["flower"]=(date-pd.to_datetime(l.year.astype(str)+"-01-01")).dt.days+1
        l["species"]=l.param.map(species)
        flowers.append(l[["station","year","species","flower"]])
    flowering=pd.concat(flowers)
    assert not flowering.duplicated(["station","year","species"]).any()
    panel=[]
    for year,g in flowering.merge(stations[["station","lat","lon"]],validate="m:1").groupby("year"):
        xy=g[["lon","lat"]].drop_duplicates().copy()
        for path,col in [(cache/f"GIMMS3g_v1_SOS_{year}_Method_4_30_90.tif","sos"),
                         (next(lights.glob(f"*{year}_{'calDMSP' if year<=2013 else 'simVIIRS'}.tif")),"dn")]:
            with rasterio.open(path) as s:
                xy[col]=[float(v[0]) if not np.ma.getmaskarray(v)[0] else np.nan for v in s.sample(xy[["lon","lat"]].to_numpy(),masked=True)]
        panel.append(g.merge(xy,validate="m:1"))
    d=pd.concat(panel).dropna()
    d=d[d.sos.between(1,240)&d.dn.between(0,63)].copy()
    d["signed_gap"]=d.flower-d.sos
    d["abs_gap"]=d.signed_gap.abs()
    e,n=Transformer.from_crs(4326,2056,always_xy=True).transform(d.lon,d.lat)
    d["block"]=(e//25000).astype(int)*100+(n//25000).astype(int)
    d["gx"]=(e//5000).astype(int)
    d["gy"]=(n//5000).astype(int)
    return d

if __name__=="__main__":
    cache,lights,out=map(Path,sys.argv[1:4])
    out.mkdir(exist_ok=True,parents=True)
    d=panel_at_stations(cache,lights)
    d.to_parquet(cache/"aligned.parquet",index=False)
    r=summarize(d)
    r.to_csv(out/"flowering_models.csv",index=False)
    print(len(d),d.station.nunique(),d.species.nunique(),len(d.drop_duplicates(["station","year"])))
    print(r[r.outcome.eq("abs_gap")][["species","days_per_10dn","lo","hi","q"]].to_string(index=False))

    full=panel_at_stations(cache,lights,"flowering (50%)")
    full.to_parquet(cache/"full_flower_aligned.parquet",index=False)
    r=summarize(full)
    r.to_csv(out/"full_flower_models.csv",index=False)
    print("50% flowering",len(full),full.station.nunique(),full.year.min(),full.year.max(),full.species.nunique(),
          r[r.outcome.eq("abs_gap")].q.min())

