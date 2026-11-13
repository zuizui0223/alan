from alan_extended11 import *
from alan_long_ndvi_cogs import BANDS
import math

def onsets():
    rows=[]
    for t in targets.itertuples():
        for y in YEARS:
            path=C/f"{t.plant}_{y}.json"
            v=pd.DataFrame(json.loads(path.read_bytes()))
            if len(v):v=v[v.site_id.eq(t.station)&v.phenophase_status.isin([0,1])].copy()
            pos=set();on=[]
            if len(v):
                v["day"]=pd.to_datetime(v.observation_date).dt.dayofyear
                assert v.latitude.sub(t.latitude).abs().max()<1e-5 and v.longitude.sub(t.longitude).abs().max()<1e-5
                pos=set(v.loc[v.phenophase_status.eq(1),"day"])
                for _,g in v.groupby("individual_id"):
                    if t.plant=="Hamamelis virginiana":g=g[g.day>=183]
                    yes=g.loc[g.phenophase_status.eq(1),"day"].min()
                    no=g.loc[g.phenophase_status.eq(0)&g.day.lt(yes),"day"].max()
                    if yes-no<=30:on.append((yes+no)/2)
            rows.append(dict(station=t.station,plant=t.plant,land=t.land_class,region=t.region5,year=y,
                onset=np.median(on) if on else None,open_days="|".join(map(str,sorted(pos))),obs=len(v)))
    return pd.DataFrame(rows)

def catalog(tile,lat,lon,y):
    params={"collections":"modis-13Q1-061","bbox":f"{lon-.00001},{lat-.00001},{lon+.00001},{lat+.00001}",
            "datetime":f"{y}-01-01/{y}-12-31","limit":500,"fields":"id,links"}
    data=get(C/f"catalog_{tile}_{y}.json",API+"/stac/v1/search",params)
    assert not any(l["rel"]=="next" for l in data["links"])
    ids=sorted(r["id"] for r in data["features"] if r["id"].startswith(f"MOD13Q1.A{y}") and r["id"].split(".")[2]==tile)
    return {int(i.split(".")[1][5:]):i for i in ids}

def fetch(job):
    st,i=job; lat,lon=coords.loc[st,["latitude","longitude"]]
    path=C/f"px_{st}_{i}.json"
    src=F/"long_ndvi"/f"ndvi_{st}_{i}.json"
    if src.exists() or path.exists():return
    d=get(path,API+f"/data/v1/item/point/{lon},{lat}",{"collection":"modis-13Q1-061","item":i,"assets":BANDS})
    assert d["band_names"]==[b+"_b1" for b in BANDS]
def jobs(panel):
    output=[]
    tileyears={}
    for st,y in panel.dropna(subset=["onset"])[["station","year"]].drop_duplicates().itertuples(index=False):
        lat,lon=coords.loc[st,["latitude","longitude"]]
        tile=f"h{int((lon*math.cos(math.radians(lat))+180)//10):02d}v{int((90-lat)//10):02d}"
        if (tile,y) not in tileyears:tileyears[tile,y]=catalog(tile,lat,lon,y)
        output.extend((st,i) for d,i in tileyears[tile,y].items() if d<=240)
    return output
if __name__=="__main__":
    panel=onsets()
    panel.to_parquet(C/"panel.parquet")
    tasks=jobs(panel)
    print("satellite work",len(tasks),flush=True)
    batch(fetch,tasks,18)
