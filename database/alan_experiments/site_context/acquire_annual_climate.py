"""Acquire ERA5 context for verified coordinates and recorded experiment years."""
from pathlib import Path
import calendar
from datetime import datetime, timezone
import hashlib
import json
import math
import numpy as np
import pandas as pd
import requests

OUT = Path(__file__).resolve().parent
CACHE = OUT.parents[3] / "gbif_annual/alan_experiment_sources/annual_climate"
CACHE.mkdir(parents=True, exist_ok=True)
SOURCE = OUT / "primary_coordinates.csv"
URL = "https://archive-api.open-meteo.com/v1/archive"

def distance_km(a,b,c,d):
    a,b,c,d=map(math.radians,[a,b,c,d])
    h=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 6371*2*math.asin(min(1,math.sqrt(h)))

def main():
    sites=pd.read_csv(SOURCE)
    daily, annual, requests_log=[],[],[]
    for _,s in sites.dropna(subset=["lat","lon"]).iterrows():
        years=str(s.experiment_year).split("-")
        first,last=int(years[0]),int(years[-1])
        assert 1940<=first<=last<=2025
        params=dict(latitude=float(s.lat),longitude=float(s.lon),start_date=f"{first}-01-01",end_date=f"{last}-12-31",
                    daily="temperature_2m_mean,precipitation_sum",models="era5",timezone="GMT",elevation="nan",cell_selection="nearest")
        key=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()[:20]
        path=CACHE/(key+".json"); meta=CACHE/(key+"_request.json")
        if not path.exists():
            response=requests.get(URL,params=params,timeout=30)
            response.raise_for_status()
            payload=response.json()
            assert not payload.get("error") and "daily" in payload
            path.write_bytes(response.content)
            meta.write_text(json.dumps(dict(url=response.url,params=params,retrieved_utc=datetime.now(timezone.utc).isoformat()),indent=2),encoding="utf-8")
        j=json.loads(path.read_text(encoding="utf-8"))
        request=json.loads(meta.read_text(encoding="utf-8"))
        assert request["params"]==params
        assert j["daily_units"]["temperature_2m_mean"]=="°C" and j["daily_units"]["precipitation_sum"]=="mm"
        assert j["utc_offset_seconds"]==0
        d=pd.DataFrame(j["daily"])
        expected=pd.date_range(f"{first}-01-01",f"{last}-12-31").strftime("%Y-%m-%d").tolist()
        assert d.time.tolist()==expected
        assert np.isfinite(d[["temperature_2m_mean","precipitation_sum"]].to_numpy()).all()
        assert d.precipitation_sum.ge(0).all()
        common=dict(study_id=s.study_id,programme_id=s.programme_id,site_id=s.site_id,
                    requested_lat=s.lat,requested_lon=s.lon,grid_lat=j["latitude"],grid_lon=j["longitude"],
                    grid_distance_km=distance_km(s.lat,s.lon,j["latitude"],j["longitude"]),model="era5",timezone="GMT")
        d["site_id"]=s.site_id;d["year"]=pd.to_datetime(d.time).dt.year
        daily.append(d)
        for year,g in d.groupby("year"):
            assert len(g)==365+calendar.isleap(year)
            annual.append(dict(**common,year=int(year),days=len(g),temperature_mean_C=float(g.temperature_2m_mean.mean()),precipitation_total_mm=float(g.precipitation_sum.sum()),period="calendar_year_not_exact_experiment_window",status="context_only_not_hydrology_or_final_moderator"))
        requests_log.append(dict(**common,request_url=request["url"],retrieved_utc=request["retrieved_utc"],raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),cache_file=path.name))
        print(s.site_id,first,last,flush=True)
    a=pd.DataFrame(annual);a.to_csv(OUT/"annual_climate.csv",index=False)
    pd.concat(daily,ignore_index=True).to_csv(OUT/"daily_climate.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    pd.DataFrame(requests_log).to_csv(OUT/"climate_requests.csv",index=False)
    summary=dict(site_years=len(a),sites=a.site_id.nunique(),programmes=a.programme_id.nunique(),
                 unique_grid_years=len(a[["grid_lat","grid_lon","year"]].drop_duplicates()),
                 source_coordinates_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                 source="Open-Meteo ERA5 historical API",documentation="https://open-meteo.com/en/docs/historical-weather-api",
                 limitations="Calendar-year context only; not measured plot wetness, inundation or baseline illumination. Shared grid cells are not independent climate observations. No meta-regression fitted.")
    (OUT/"climate_audit.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
