"""Join repeated open-flower status, QA-screened Terra NDVI, VIIRS and adult moth dates."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from alan_long_ndvi_download import STATIONS, YEARS
from alan_long_ndvi_cogs import BANDS
from mois_global_multipoint_numeric_20261107 import MONTHS, extract_samples

OUT=Path(__file__).resolve().parents[2]/"results/long_phenology"
LABEL={19586:("Hamamelis virginiana","island"),2:("Hamamelis virginiana","mainland"),
       25009:("Holodiscus discolor","island"),33137:("Holodiscus discolor","mainland")}


def run(folder):
    cache=folder/"long_ndvi"
    positions=pd.read_parquet(folder/"npn/individual_visits.parquet").groupby("station").first().loc[list(STATIONS)]
    coords=positions.rename(columns={"latitude":"query_latitude","longitude":"query_longitude"}).to_dict("records")
    # Hash each input family without publishing raw station observations.
    manifest={}
    for prefix in ["flower","light","moth","modis","ndvi"]:
        paths=sorted(cache.glob(prefix+"_*.json"))
        hashes=[hashlib.sha256(f.read_bytes()).hexdigest() for f in paths]
        manifest[prefix]={"files":len(paths),"ordered_hash_digest":hashlib.sha256("\n".join(hashes).encode()).hexdigest()}
    rows=[]
    series=[]
    for year in YEARS:
        light={s:[] for s in STATIONS}
        for month in MONTHS:
            src=json.loads((cache/f"light_{year}_{month}.json").read_bytes())
            if src.get("error"): continue
            for station, item in zip(STATIONS,extract_samples("viirs",month,coords,src,year)):
                if item["value"] is not None: light[station].append(item["value"])
        for station in STATIONS:
            plant,side=LABEL[station]
            v=pd.DataFrame(json.loads((cache/f"flower_{station}_{year}.json").read_bytes()))
            positive=set()
            onsets=[]
            if len(v):
                assert v.observation_id.is_unique and v.site_id.eq(station).all()
                v["day"]=pd.to_datetime(v.observation_date).dt.dayofyear
                v=v[v.phenophase_status.isin([0,1])]
                positive=set(v.loc[v.phenophase_status.eq(1),"day"])
                for _,g in v.groupby("individual_id"):
                    if STATIONS[station]==76: g=g[g.day>=183]
                    yes=g.loc[g.phenophase_status.eq(1),"day"]
                    if yes.empty:continue
                    start=yes.min()
                    no=g.loc[g.phenophase_status.eq(0)&g.day.lt(start),"day"]
                    if len(no) and start-no.max()<=30:
                        onsets.append((start+no.max())/2)
            nd={}
            for path in cache.glob(f"ndvi_{station}_MOD13Q1.A{year}*.json"):
                nd[int(path.name.split(".")[1][5:])]=json.loads(path.read_bytes())["values"]
            for path in cache.glob(f"modis_{station}_{year}_*.json"):
                by={}
                for rec in json.loads(path.read_bytes())["subset"]:
                    by.setdefault(int(rec["modis_date"][5:]),{})[rec["band"]]=rec["data"][0]
                for k,b in by.items():
                    nd[k]=[b.get(name) for name in BANDS]
            good={int(day):raw/10000 for raw,qa,day in nd.values() if None not in (raw,qa,day) and
                  qa==0 and -2000<=raw<=10000 and 1<=day<=366}
            for comp,(raw,qa,day) in sorted(nd.items()):
                series.append(dict(plant=plant,land=side,year=year,composite_start_doy=comp,
                    acquisition_doy=day,reliability=qa,ndvi=raw/10000 if raw is not None and -2000<=raw<=10000 else None))
            # Bracket the first spring crossing of the local winter-to-summer NDVI midpoint.
            winter=[v for day,v in good.items() if day<=80]
            summer=[v for day,v in good.items() if 160<=day<=240]
            green=np.nan
            if len(winter)>=2 and len(summer)>=2 and np.median(summer)>np.median(winter):
                threshold=(np.median(summer)+np.median(winter))/2
                days=sorted(good)
                for a,b in zip(days,days[1:]):
                    if 45<=a<b<=180 and b-a<=32 and good[a]<threshold<=good[b]:
                        green=a+(b-a)*(threshold-good[a])/(good[b]-good[a])
                        break
            near=[good[min(good,key=lambda x:(abs(x-p),x))] for p in sorted(positive) if good and min(abs(x-p) for x in good)<=16]
            mothdays=set()
            for path in cache.glob(f"moth_{station}_*.json"):
                for m in json.loads(path.read_bytes())["results"]:
                    if m["year"]==year and m.get("day") and m.get("lifeStage")=="Adult":
                        mothdays.add(pd.Timestamp(year,m["month"],m["day"]).dayofyear)
            rows.append(dict(plant=plant,land=side,year=year,visits=len(v),open_days=len(positive),
                first_open=min(positive,default=None),last_open=max(positive,default=None),
                onset_midpoint=np.median(onsets) if onsets else None,onset_individuals=len(onsets),
                spring_ndvi_midpoint=green,flower_minus_green_days=np.median(onsets)-green if onsets else None,
                returned_composites=len(nd),good_ndvi_days=len(good),open_days_with_good_ndvi=len(near),
                ndvi_on_open_dates=np.median(near) if near else None,viirs_months=len(light[station]),
                viirs=np.median(light[station]) if len(light[station])>=3 else None,
                adult_moth_days=len(mothdays),same_day_records=len(mothdays & positive),
                closest_record_gap_days=min((abs(m-p) for m in mothdays for p in positive),default=None)))
    frame=pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    pd.DataFrame(series).to_csv(OUT/"ndvi_16day.csv",index=False)
    (OUT/"source_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    frame.to_csv(OUT/"station_year.csv",index=False)
    print(frame.groupby(["plant","land"]).agg(n_years=("year","size"),visits=("visits","sum"),
          open_years=("first_open","count"),onset_years=("onset_midpoint","count"),satellite=("returned_composites","sum"),
          flowering_ndvi=("ndvi_on_open_dates","count"),light=("viirs","count"),moth_days=("adult_moth_days","sum")).to_string())
    print(frame.dropna(subset=["closest_record_gap_days"])[["year","closest_record_gap_days"]])
    print("triple",frame.dropna(subset=["ndvi_on_open_dates","viirs"]).groupby(["plant","land"]).size())
    print("paired onset",frame.dropna(subset=["onset_midpoint","ndvi_on_open_dates","viirs"]).groupby(["plant","year"]).land.nunique().loc[lambda x:x==2])


if __name__=="__main__":
    run(Path(sys.argv[1]))
