from alan_extended11_ndvi import *
from alan_extended11_moths import links
from scipy.stats import theilslopes
import statsmodels.formula.api as smf

def greens(panel):
    series={}
    for st,i in jobs(panel):
        token=i.split(".")[1];year=int(token[1:5]); day=int(token[5:])
        path=C/f"px_{st}_{i}.json"
        if not path.exists():path=F/"long_ndvi"/f"ndvi_{st}_{i}.json"
        if not path.exists():raise RuntimeError(f"missing {st} {i}")
        data=json.loads(path.read_bytes())
        raw,qa,doy=data["values"]
        if qa==0 and raw is not None and -2000<=raw<=10000 and doy is not None and 1<=doy<=366:
            series.setdefault((st,year),{})[int(doy)]=raw/10000
    green={}
    for key,v in series.items():
        w=[x for d,x in v.items() if d<=80];s=[x for d,x in v.items() if 160<=d<=240]
        if min(len(w),len(s))<2 or np.median(s)<=np.median(w):continue
        half=(np.median(w)+np.median(s))/2
        days=sorted(v)
        for a,b in zip(days,days[1:]):
            if 45<=a<b<=180 and b-a<=32 and v[a]<half<=v[b]:
                green[key]=a+(b-a)*(half-v[a])/(v[b]-v[a]);break
    return green

def lights():
    vals={}
    for y in YEARS:
        for k in range(0,len(coords),80):
            group=coords.iloc[k:k+80]
            rows=group.rename(columns={"longitude":"query_longitude","latitude":"query_latitude"}).to_dict("records")
            for m in MONTHS:
                d=json.loads((C/f"lights_{y}_{m}_{k}.json").read_bytes())
                if "error" in d:continue
                for st,item in zip(group.index,extract_samples("viirs",m,rows,d,y)):
                    if item["value"] is not None:vals.setdefault((st,y),[]).append(item["value"])
    return {k:np.median(v) for k,v in vals.items() if len(v)>=3}

def fly():
    days={}
    for t in targets.itertuples():
        names=links[t.plant]
        for f in C.glob(f"gbif_{t.station}_*.json"):
            for rec in json.loads(f.read_bytes())["results"]:
                if rec.get("species") not in names or rec.get("lifeStage")!="Adult" or not rec.get("day"):continue
                k=t.station,t.plant,rec["year"]
                days.setdefault(k,set()).add(pd.Timestamp(rec["year"],rec["month"],rec["day"]).dayofyear)
    return days

def report():
    p=pd.read_parquet(C/"panel.parquet")
    green,light,flight=greens(p),lights(),fly()
    p["green"]=[green.get((r.station,r.year)) for r in p.itertuples()]
    p["radiance"]=[light.get((r.station,r.year)) for r in p.itertuples()]
    p["gap"]=p.onset-p.green
    p["loglight"]=np.log1p(p.radiance.clip(lower=0))
    p["unit"]=p.station.astype(str)+"|"+p.plant
    p["flight_days"]=[len(flight.get((r.station,r.plant,r.year),set())) for r in p.itertuples()]
    p["overlap_days"]=[len(flight.get((r.station,r.plant,r.year),set()) & {int(d) for d in r.open_days.split("|") if d}) for r in p.itertuples()]
    event_rows=[]
    for t in targets.itertuples():
        for f in C.glob(f"gbif_{t.station}_*.json"):
            for rec in json.loads(f.read_bytes())["results"]:
                if rec.get("species") in links[t.plant] and rec.get("lifeStage")=="Adult" and rec.get("day"):
                    event_rows.append(dict(station=t.station,plant=t.plant,land=t.land_class,moth=rec["species"],date=rec["eventDate"][:10],gbif_key=rec["key"]))
    pd.DataFrame(event_rows).drop_duplicates().to_csv(O/"moth_dates.csv",index=False)
    p["nearest_days"]=[min((abs(d-int(e)) for d in flight.get((r.station,r.plant,r.year),set())
                for e in r.open_days.split("|") if e),default=None) for r in p.itertuples()]
    a=p.dropna(subset=["gap","loglight"]).copy()
    a=a[a.groupby("unit").year.transform("size").ge(2)].copy()
    a["island"]=a.land.eq("island").astype(int)
    a["islandlight"]=a.island*a.loglight
    trend=[]
    for u,g in a.groupby("unit"):
        if len(g)>=3:
            trend.append(dict(unit=u,plant=g.plant.iloc[0],land=g.land.iloc[0],station=g.station.iloc[0],years=len(g),
                light_per_year=theilslopes(g.loglight,g.year).slope,gap_per_year=theilslopes(g.gap,g.year).slope))
    trends=pd.DataFrame(trend)
    trendfits=[]
    if len(trends):
        for side,g in trends.groupby("land"):
            if g.station.nunique()>=3 and g.light_per_year.nunique()>=2:
                fit=smf.ols("gap_per_year ~ light_per_year",g).fit(cov_type="cluster",cov_kwds={"groups":g.station},use_t=True)
                lo,hi=fit.conf_int().loc["light_per_year"]
                trendfits.append(dict(land=side,stations=g.station.nunique(),species_stations=len(g),
                     slope=fit.params.light_per_year,low=lo,high=hi))
        pd.DataFrame(trendfits).to_csv(O/"trend_association.csv",index=False)
        print("trend association",trendfits)
    contrasts=[]
    if a.island.nunique()==2:
        from numpy.linalg import matrix_rank
        # Explicit formula environment avoids module cache-path C shadowing Patsy's C().
        mod=smf.ols("gap ~ loglight + islandlight + C(unit) + C(year)",data=a,eval_env=-1)
        assert matrix_rank(mod.exog)==mod.exog.shape[1]
        fit=mod.fit(cov_type="cluster",cov_kwds={"groups":a.station},use_t=True)
        for label,f in [("mainland","loglight"),("island","loglight + islandlight"),("island_minus_mainland","islandlight")]:
            test=fit.t_test(f+"=0")
            effect=test.effect.item();lo,hi=test.conf_int()[0]
            contrasts.append(dict(contrast=label,days_per_log1p=effect,low=lo,high=hi))
    complete=p.dropna(subset=["gap","radiance"])
    bilateral=complete.groupby(["plant","region","year"]).land.nunique().eq(2)
    bilateral[bilateral].to_csv(O/"bilateral_years.csv")
    print("complete",len(complete),"bilateral",bilateral.sum(),bilateral[bilateral].index.get_level_values("plant").unique().tolist())
    p.to_csv(O/"panel.csv",index=False)
    trends.to_csv(O/"longitudinal_slopes.csv",index=False)
    # Brightening expressed as a twofold increase in (1 + radiance).
    for c in contrasts:
        c["days_per_doubling"] = np.log(2)*c["days_per_log1p"]
        c["doubling_low"] = np.log(2)*c["low"]
        c["doubling_high"] = np.log(2)*c["high"]
    pd.DataFrame(contrasts).to_csv(O/"within_site.csv",index=False)
    counts=p.groupby(["plant","land"]).agg(stations=("station","nunique"),onsets=("onset","count"),flower_and_flight=("nearest_days","count"))
    counts["complete_years"]=p.dropna(subset=["gap","radiance"]).groupby(["plant","land"]).size()
    counts["repeat_years"]=a.groupby(["plant","land"]).size()
    counts.fillna(0).astype(int).to_csv(O/"coverage.csv")
    print(counts,flush=True)
    print("within",a.groupby("land").agg(rows=("year","size"),stations=("station","nunique"),plants=("plant","nunique")),contrasts)
    print("trends",trends.groupby("land").size(),"flight",p.nearest_days.eq(0).sum(),p.nearest_days.notna().sum())
    if len(trends):
        print(trends.groupby("land")[["light_per_year","gap_per_year"]].corr(method="spearman"))
if __name__=="__main__":report()
