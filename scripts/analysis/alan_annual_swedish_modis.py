"""Fine-scale NDVI match for candidate species observed on both Swedish islands and mainland."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys, math, json, time
import numpy as np
import pandas as pd
import requests
from alan_annual_flowering import sample
import geopandas as gpd
from mois_global_2020_GSHHG_scope_20261107 import prepare_land,assign_points

API="https://planetarycomputer.microsoft.com/api"
BANDS=["250m_16_days_NDVI","250m_16_days_pixel_reliability","250m_16_days_composite_day_of_the_year"]


def main(cache,lights,coast,out):
    folder=cache/"swedish_modis";folder.mkdir(exist_ok=True)
    a=pd.read_csv(out/"candidate_flowering.csv")
    names=set(a[a.priority.eq("FDR")].species)
    cross=pd.read_csv(out/"flowering_name_crosswalk.csv")
    names.update(cross.loc[cross.sp.isin(a.loc[a.priority.eq("FDR"),"sp"]),"field_name"])
    d=pd.read_csv(cache/"sweden"/"plants_calendar_2008-2024.csv").query("quality==1 and phase_id==53")
    d=d[d.scientific_name.isin(names)].rename(columns={"obs_year":"year","scientific_name":"name","obs_doy":"flower"})
    d=d.groupby(["latitude","longitude","name","year"],as_index=False).flower.median()
    xy=d[["longitude","latitude"]].drop_duplicates().sort_values(["longitude","latitude"]).reset_index(drop=True)
    xy["unit_id"]=xy.index
    land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
    xy=xy.merge(assign_points(xy,land,gpd.read_file(coast/"GSHHS_h_L2.shp"))[["unit_id","land_class","land_id"]],validate="1:1")
    d=d.merge(xy,validate="m:1").query("land_class in ['island','mainland']")
    supported=d.groupby("name").land_class.nunique()
    d=d[d.name.isin(supported[supported.eq(2)].index)]
    coords=d[["year","latitude","longitude","unit_id"]].drop_duplicates().reset_index(drop=True)
    for year,g in coords.groupby("year"):
        coords.loc[g.index,"dn"]=sample(g,next(lights.glob(f"Harmonized_DN_NTL_{year}_*.tif")))
    def get(path,url,params):
        if path.exists():return json.loads(path.read_bytes())
        for i in range(4):
            try:
                r=requests.get(url,params=params,timeout=35);r.raise_for_status()
                obj=r.json();path.write_bytes(r.content);return obj
            except requests.RequestException:
                if i==3:raise
                time.sleep(1+i)
    def tile(q):
        return f"h{int((q.longitude*math.cos(math.radians(q.latitude))+180)//10):02d}v{int((90-q.latitude)//10):02d}"
    coords["tile"]=[tile(q) for q in coords.itertuples()]
    def catalog(q):
        lon,lat,y=q.longitude,q.latitude,q.year
        p={"collections":"modis-13Q1-061","bbox":f"{lon-.00001},{lat-.00001},{lon+.00001},{lat+.00001}","datetime":f"{y}-01-01/{y}-12-31","limit":500,"fields":"id,links"}
        j=get(folder/f"catalog_{q.tile}_{y}.json",API+"/stac/v1/search",p)
        assert not any(l["rel"]=="next" for l in j["links"])
        items={int(f["id"].split(".")[1][5:]):f["id"] for f in j["features"] if f["id"].startswith(f"MOD13Q1.A{y}") and f["id"].split(".")[2]==q.tile}
        assert len(items)>=20
        return (q.tile,y),[items[k] for k in sorted(items) if k<=240]
    with ThreadPoolExecutor(6) as pool: catalogs=dict(pool.map(catalog,coords.drop_duplicates(["tile","year"]).itertuples()))
    jobs=[(q,item) for q in coords.itertuples() for item in catalogs[q.tile,q.year]]
    def pixel(job):
        q,item=job
        path=folder/f"px_{q.longitude}_{q.latitude}_{item}.json"
        data=get(path,API+f"/data/v1/item/point/{q.longitude},{q.latitude}",dict(collection="modis-13Q1-061",item=item,assets=BANDS))
        assert data["band_names"]==[band+"_b1" for band in BANDS] and np.allclose(data["coordinates"],[q.longitude,q.latitude])
        value,qa,day=data["values"]
        if qa==0 and value is not None and -2000<=value<=10000 and day is not None and 1<=day<=366:return q.unit_id,q.year,day,value
    print("Swedish matched taxa, station years, pixel requests",list(d.name.unique()),len(coords),len(jobs),flush=True)
    with ThreadPoolExecutor(12) as pool: data=[v for v in pool.map(pixel,jobs) if v is not None]
    green=[]
    for (site,y),g in pd.DataFrame(data,columns=["unit_id","year","day","ndvi"]).groupby(["unit_id","year"]):
        v=g.groupby("day").ndvi.median().sort_index()
        winter=v[v.index<=80];summer=v[(v.index>=160)&(v.index<=240)]
        if min(len(winter),len(summer))<2 or winter.median()>=summer.median():continue
        half=(winter.median()+summer.median())/2
        for a,b in zip(v.index,v.index[1:]):
            if 45<=a<b<=180 and b-a<=32 and v[a]<half<=v[b]:
                green.append([site,y,a+(b-a)*(half-v[a])/(v[b]-v[a])]);break
    coords=coords.merge(pd.DataFrame(green,columns=["unit_id","year","sos"]),how="left",validate="1:1")
    d=d.merge(coords,validate="m:1")
    d["gap"]=(d.flower-d.sos).abs()
    d.to_parquet(cache/"sweden_matched.parquet",index=False)
    d.groupby(["name","land_class"]).agg(onset_station_years=("flower","size"),matched=("gap","count"),mean_abs_days=("gap","mean")).to_csv(out/"swedish_modis_coverage.csv")
    shared=d.dropna(subset=["gap"]).groupby(["name","year","land_class"]).agg(mean_gap=("gap","mean"),stations=("unit_id","nunique")).unstack()
    shared.columns=["_".join(c) for c in shared.columns]
    for col in ["mean_gap_island","mean_gap_mainland"]:
        if col not in shared:shared[col]=np.nan
    shared["island_minus_mainland"]=shared.mean_gap_island-shared.mean_gap_mainland
    shared.to_csv(out/"swedish_modis_year_comparisons.csv")
    print("Common year contrasts",shared.island_minus_mainland.dropna().to_string())

if __name__=="__main__":main(*map(Path,sys.argv[1:5]))
