"""Join measured plant flowering onset to the spring NDVI transition and harmonized DN."""
from pathlib import Path
from zipfile import ZipFile
from alan_annual_jma import read_jma
import sys
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
from pyproj import Transformer
from mois_global_2020_GSHHG_scope_20261107 import prepare_land,assign_points


def observations(old,swiss,cache):
    cols=["source","station","latitude","longitude","species","year","flower"]
    with ZipFile(old/"chronicle-of-nature-calendar.v1.0.5.zip") as archive:
        def read(name):
            return pd.read_csv(archive.open(next(f for f in archive.namelist() if f.endswith("/"+name+".csv"))),low_memory=False)
        d=read("phenology")
        d=d[d.eventtype.eq("onset of blooming")&d.quality.eq("OK")].merge(read("studysites"),validate="m:1")
        d=d.rename(columns={"studysite":"station","dayofyear":"flower","taxon":"species"}).assign(source="Chronicles")
    pep=pd.concat([pd.read_csv(p,usecols=lambda c:c in {"s_id","lon","lat","species","year","Flower_open","DOY_Flower_open"}).rename(columns={"DOY_Flower_open":"Flower_open"}) for p in
                   [old/"final_data_TG_tmp.csv",old/"final_data_SG_tmp.csv"]]).drop_duplicates()
    pep=pep.rename(columns={"s_id":"station","lon":"longitude","lat":"latitude","Flower_open":"flower"}).assign(source="PEP725")
    dwd=pd.read_parquet(cache/"dwd"/"onset.parquet")
    countries=gpd.read_file(cache/"countries.zip")
    germany=countries[countries.ADMIN.eq("Germany")].geometry.unary_union
    pts=pep[["longitude","latitude"]].drop_duplicates()
    pts["german_pep"]=gpd.points_from_xy(pts.longitude,pts.latitude).within(germany)
    pep=pep.merge(pts,validate="m:1")
    pep=pep[~(pep.german_pep & pep.species.isin(dwd.species))]
    ch=pd.read_parquet(swiss/"aligned.parquet").rename(columns={"lon":"longitude","lat":"latitude"}).assign(source="MeteoSwiss")
    usa=[]
    for year in range(2009,2015):
        f=pd.read_json(cache/f"npn_{year}.json")
        f=f[f.numdays_since_prior_no.between(1,30)&f.first_yes_year.eq(year)&f.phenophase_description.str.strip().str.lower().eq("open flowers")].copy()
        f["flower"]=f.first_yes_doy-f.numdays_since_prior_no/2
        f["species"]=f.genus+" "+f.species
        usa.append(f.rename(columns={"site_id":"station"}).assign(year=year,source="USA-NPN")[cols])
    swe_file=cache/"sweden"/"plants_calendar_2008-2024.csv"
    if not swe_file.exists():
        import requests
        swe_file.parent.mkdir(exist_ok=True)
        response=requests.get("https://api.researchdata.se/dataset/2022-86-1/2/file/data/"+swe_file.name,timeout=90)
        response.raise_for_status();swe_file.write_bytes(response.content)
    swe=pd.read_csv(swe_file).query("quality==1 and phase_id==53").rename(columns={
        "scientific_name":"species","obs_year":"year","obs_doy":"flower"})
    swe=swe.assign(source="SWE-NPN",station="coordinate")
    panel=pd.concat([swe[cols],d[cols],pep[cols],ch[cols],read_jma(cache/"jma")[cols],dwd[cols],*usa],ignore_index=True).dropna()
    panel=panel[panel.year.between(1992,2014)&panel.flower.between(1,366)&panel.latitude.ge(30)].copy()
    panel["year"]=panel.year.astype(int)
    # Station identifiers can change or be reused. Fixed sites are geographic.
    return panel.groupby(["source","latitude","longitude","species","year"],as_index=False).flower.median()


def sample(pts,path):
    values=np.full(len(pts),np.nan)
    with rasterio.open(path) as s:
        rr,cc=map(np.asarray,rasterio.transform.rowcol(s.transform,pts.longitude,pts.latitude))
        ok=(rr>=0)&(cc>=0)&(rr<s.height)&(cc<s.width)
        ix=np.flatnonzero(ok)
        h,w=s.block_shapes[0]
        for (i,j),sub in pd.Series(ix).groupby([rr[ix]//h,cc[ix]//w]):
            win=s.block_window(1,int(i),int(j))
            ids=sub.values
            v=s.read(1,window=win,masked=True)[rr[ids]-int(win.row_off),cc[ids]-int(win.col_off)]
            values[ids]=v.astype(float).filled(np.nan)
    return values


def main(cache,old,swiss,lights,coast):
    d=observations(old,swiss,cache)
    xy=d[["longitude","latitude"]].drop_duplicates().reset_index(drop=True)
    xy["unit_id"]=xy.index
    land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
    xy=xy.merge(assign_points(xy,land,gpd.read_file(coast/"GSHHS_h_L2.shp"))[["unit_id","land_class","land_id"]],validate="1:1")
    east,north=Transformer.from_crs(4326,6933,always_xy=True).transform(xy.longitude,xy.latitude)
    xy["block"]=(np.floor(east/1e5)*1000+np.floor(north/1e5)).astype(int)
    d=d.merge(xy,validate="m:1")
    samples=[]
    for year,g in d.groupby("year"):
        pts=g[["unit_id","longitude","latitude"]].drop_duplicates()
        pts["dn"]=sample(pts,next(lights.glob(f"Harmonized_DN_NTL_{year}_*.tif")))
        pts["sos"]=sample(pts,swiss/f"GIMMS3g_v1_SOS_{year}_Method_4_30_90.tif")
        samples.append(pts.assign(year=year)[["year","unit_id","sos","dn"]])
    d=d.merge(pd.concat(samples),validate="m:1").query("1<=sos<=240 and 0<=dn<=63")
    d=d[d.land_class.isin(["island","mainland"])].copy()
    d["signed_gap"]=d.flower-d.sos
    d["gap"]=d.signed_gap.abs()
    d.to_parquet(cache/"flowering.parquet",index=False)
    print(d.groupby(["source","land_class"]).agg(rows=("gap","size"),stations=("unit_id","nunique"),taxa=("species","nunique")),flush=True)


if __name__=="__main__":main(*map(Path,sys.argv[1:6]))
