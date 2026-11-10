#!/usr/bin/env python3
"""Original 2026-08 GBIF gecko & skink timestamped day/night photo source.

Select original Squamata Parquet observations: 2013-25, HUMAN_OBSERVATION,
PRESENT, real coordinates, one of six gecko families or Scincidae, an actual
explicit-offset ISO 8601 eventtime and source year/month/day. Convert the
time to UTC, compute approximate NOAA solar altitude for actual source
lat/lon/time, then exact GSHHG 2.3.7 island point-in-polygon BEFORE
aggregating to original island × 0.05 degree cell × year × family ×
civil-night vs solar-day photographic counts.

This is observed photographic/reporting time, NOT continuous lizard activity,
insect consumption, light-at-gecko body or night flower pollination.
Original source clocks may be wrong or biased and only a subset has
time-zone-qualified times. Does not publish individual timestamp/coordinate.
"""
from __future__ import annotations
import argparse,collections,datetime,hashlib,json,math,re
from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pyarrow.fs import S3FileSystem
from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land
from scripts.analysis.alan_global_GBIF_gecko_island_original_stream_20261110 import (
    GECKO_FAMILIES)

SOURCE_DOI="10.15468/dl.8yrbe7"
PARTITION="taxonpartition=Animalia_Chordata_Reptilia_Squamata"
SCHEMA="ALAN_original_202608_GSHHG_island_gecko_skink_true_UTC_sun_alitude_photo_cellyear_v1"
FAMILIES=tuple(GECKO_FAMILIES)+("Scincidae",)
SOLAR_NIGHT=-6.0
SOLAR_DAY=6.0
YEAR_MIN,YEAR_MAX=2013,2025
GRID=.05

def require(c,why):
    if not bool(c):raise ValueError("ALAN_TRUE_GBIF_GECKO_SUN: "+str(why))

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def source_data(manifest):
    obj=json.loads(Path(manifest).read_text())
    require(obj["schema"]=="ALAN_original_GBIF_202608_Gekkota_reptile_Parquet_partition_S3_audit_v1"
            and obj["original_GBIF_snapshot_DOI"]==SOURCE_DOI,
            "Wrong genuine original 2026 GBIF source snapshot")
    return obj["matching_reptile_gecko_partition_folders"][PARTITION]["original_source_object_metadata"]

def bounds(land):
    coarse=np.zeros((180,360),dtype=bool)
    for a,b,c,d in land.geometry.bounds.itertuples(index=False,name=None):
        a,b,c,d=[float(v) for v in (a,b,c,d)]
        if not np.isfinite([a,b,c,d]).all():continue
        xmin=int(np.clip(np.floor(a+180),0,359))
        xmax=int(np.clip(np.floor(c+180),0,359))
        ymin=int(np.clip(np.floor(b+90),0,179))
        ymax=int(np.clip(np.floor(d+90),0,179))
        if xmax>=xmin and ymax>=ymin:
            coarse[ymin:ymax+1,xmin:xmax+1]=True
    return coarse

def noaa_solar_altitude(lat,lon,UTC):
    """Approximate solar geometric elevation from validated aware UTC times.

    Mean solar position (NOAA series), not cloudiness, irradiation, artificial
    luminance, photoreceptor exposure or a photo of gecko feeding.
    """
    la=np.deg2rad(np.asarray(lat,dtype=float))
    lo=np.asarray(lon,dtype=float)
    require(np.isfinite(la).all() and np.isfinite(lo).all(),
            "Original latitude/longitude nonfinite")
    utc=pd.DatetimeIndex(UTC)
    require(utc.tz is not None,
            "Must reject clock times without explicit timezone or UTC normalization")
    doy=utc.dayofyear.to_numpy(dtype=float)
    hh=(utc.hour.to_numpy(dtype=float)+
        utc.minute.to_numpy(dtype=float)/60+
        utc.second.to_numpy(dtype=float)/3600)
    g=2*np.pi/365*(doy-1+(hh-12)/24)
    eq=229.18*(.000075+.001868*np.cos(g)-.032077*np.sin(g)-
               .014615*np.cos(2*g)-.040849*np.sin(2*g))
    dec=(.006918-.399912*np.cos(g)+.070257*np.sin(g)-
         .006758*np.cos(2*g)+.000907*np.sin(2*g)-
         .002697*np.cos(3*g)+.00148*np.sin(3*g))
    solar=(60*hh+eq+4*lo)%1440.
    h=np.deg2rad(solar/4-180)
    sin_elevation=(np.sin(la)*np.sin(dec)+
                   np.cos(la)*np.cos(dec)*np.cos(h))
    return np.rad2deg(np.arcsin(np.clip(sin_elevation,-1.,1.)))

def original_clock(year,month,day,evtime):
    """Only parse a genuine ISO hour with explicit numeric offset/Z."""
    source=strings=evtime.fillna("").astype(str).str.strip()
    valid=source.str.match(
        r"^(?:[01]\d|2[0-3]):[0-5]\d(?::[0-5]\d(?:\.\d+)?)?(?:Z|[+-](?:0\d|1\d|2[0-3]):?[0-5]\d)$",
        na=False)
    yr=pd.to_numeric(year,errors="coerce")
    mo=pd.to_numeric(month,errors="coerce")
    da=pd.to_numeric(day,errors="coerce")
    valid &= yr.between(YEAR_MIN,YEAR_MAX)&mo.between(1,12)&da.between(1,31)
    if not valid.any():
        return pd.Series(pd.NaT,index=source.index,dtype="datetime64[ns, UTC]")
    y=yr.loc[valid].astype(int).astype(str)
    m=mo.loc[valid].astype(int).astype(str).str.zfill(2)
    d=da.loc[valid].astype(int).astype(str).str.zfill(2)
    times=pd.to_datetime(y+"-"+m+"-"+d+"T"+source.loc[valid],
                         format="mixed",utc=True,errors="coerce")
    out=pd.Series(pd.NaT,index=source.index,dtype="datetime64[ns, UTC]")
    out.loc[valid]=times
    return out

def run(manifest,coast,out):
    out=Path(out);require(not out.exists(),"Cannot overwrite authentic global sun-time reptile audit")
    sources=source_data(manifest)
    land=load_signed_land(coast)
    land=land.loc[land.land_class.eq("island")&
                  land.area_km2.ge(20),["land_id","geometry"]].copy()
    land["land_id"]=land.land_id.astype(str)
    require(len(land)==3696,
            "Original source GSHHG island geometry differs from existing moth study")
    bbox=bounds(land)
    fs=S3FileSystem(anonymous=True,region="us-east-1")
    tally=collections.Counter()
    family=collections.defaultdict(collections.Counter)
    cellcount=collections.Counter()
    files=[]
    for source in sources:
        file=fs.open_input_file("gbif-public-data/"+source["key"])
        pf=pq.ParquetFile(file)
        fields=list(pf.schema_arrow.names)
        need=("family","year","month","day","eventtime",
              "basisofrecord","occurrencestatus",
              "decimallatitude","decimallongitude")
        aliases={}
        for name in need:
            match=[c for c in fields if c.casefold()==name.casefold()]
            require(len(match)==1,
                    "Actual original GBIF Squamata clock/coordinate field missing "+name)
            aliases[name]=match[0]
        files.append({"source_key":source["key"],
                      "original_parquet_rows":pf.metadata.num_rows,
                      "source_row_groups":pf.metadata.num_row_groups})
        for rg in range(pf.metadata.num_row_groups):
            a=pf.read_row_group(rg,columns=list(aliases.values())).to_pandas()
            tally["original_all_family_all_year_reptile_rows_read"]+=len(a)
            year=pd.to_numeric(a[aliases["year"]],errors="coerce")
            fam=a[aliases["family"]].fillna("").astype(str)
            good=(year.between(YEAR_MIN,YEAR_MAX).to_numpy()&
                  fam.isin(FAMILIES).to_numpy()&
                  a[aliases["basisofrecord"]].fillna("").astype(str)
                    .str.upper().eq("HUMAN_OBSERVATION").to_numpy()&
                  a[aliases["occurrencestatus"]].fillna("").astype(str)
                    .str.upper().eq("PRESENT").to_numpy())
            d=a.loc[good].copy()
            tally["source_reptile_gecko_skink_HUMAN_PRESENT_2013to25"]+=len(d)
            if d.empty:continue
            lat=pd.to_numeric(d[aliases["decimallatitude"]],errors="coerce")
            lon=pd.to_numeric(d[aliases["decimallongitude"]],errors="coerce")
            ok=lat.between(-90,90)&lon.between(-180,180)
            d=d.loc[ok].copy()
            lat=lat.loc[ok].to_numpy(float);lon=lon.loc[ok].to_numpy(float)
            tally["genuine_georeferenced_reptile_source"]+=len(d)
            if d.empty:continue
            dt=original_clock(d[aliases["year"]],d[aliases["month"]],
                              d[aliases["day"]],d[aliases["eventtime"]])
            goodTime=dt.notna()
            tally["genuine_georeferenced_explicit_UTC_offset_eventtime"]+=int(goodTime.sum())
            if not goodTime.any():continue
            d=d.loc[goodTime].copy()
            lat=lat[goodTime.to_numpy()];lon=lon[goodTime.to_numpy()]
            timestamp=dt.loc[goodTime]
            near=bbox[np.minimum(179,np.maximum(0,np.floor(lat+90).astype(int))),
                      np.minimum(359,np.maximum(0,np.floor(lon+180).astype(int)))]
            d=d.loc[near].copy()
            lat=lat[near];lon=lon[near]
            timestamp=timestamp.loc[near]
            tally["valid_time_coarse_near_original_20km2_islands"]+=len(d)
            if d.empty:continue
            elev=noaa_solar_altitude(lat,lon,timestamp)
            d["true_photo_source_family"]=d[aliases["family"]].astype(str)
            d["recorded_source_year"]=year.loc[d.index].astype(int).to_numpy()
            d["latbin"]=np.minimum(3599,np.floor((lat+90)/GRID).astype(int))
            d["lonbin"]=np.minimum(7199,np.floor((lon+180)/GRID).astype(int))
            d["source_solar_class"]=np.where(elev<SOLAR_NIGHT,"solar_night",
                np.where(elev>SOLAR_DAY,"solar_day","twilight_excluded"))
            points=gpd.GeoDataFrame(d[[
               "true_photo_source_family","recorded_source_year",
               "latbin","lonbin","source_solar_class"]],
               geometry=gpd.points_from_xy(lon,lat),crs=4326)
            matched=gpd.sjoin(points,land,how="inner",predicate="within")
            dup=matched.index.duplicated(keep=False)
            tally["original_reptile_records_ambiguous_island_polygon"]+=int(dup.sum())
            matched=matched.loc[~dup].copy()
            matched["island"]=matched.land_id.astype(str)
            tally["genuine_island_source_photo_clock_records"]+=len(matched)
            for fam_name,g in matched.groupby("true_photo_source_family"):
                v=g.source_solar_class.value_counts()
                for category,count in v.items():
                    family[fam_name][category]+=int(count)
            group=["island","recorded_source_year","latbin",
                   "lonbin","true_photo_source_family","source_solar_class"]
            for key,n in matched.groupby(group,sort=False).size().items():
                cellcount[key]+=int(n)
            if tally["original_all_family_all_year_reptile_rows_read"]%1000000<150000:
                print("ALAN_TRUE_GECKO_ISLAND_CLOCK_SOURCE_PROGRESS",json.dumps({
                  "original_reptile_rows":tally["original_all_family_all_year_reptile_rows_read"],
                  "true_island_clock_source":tally["genuine_island_source_photo_clock_records"]},
                  sort_keys=True),flush=True)
        file.close()
    require(tally["original_all_family_all_year_reptile_rows_read"]==11634053 and
            sum(x["original_parquet_rows"] for x in files)==11634053,
            "Full original 2026-08 GBIF Squamata scan not complete")
    out.mkdir(parents=True)
    rows=[{"island":island,"year":int(year),"latbin":int(y),
           "lonbin":int(x),"source_taxon_family":fam,
           "verified_solar_photo_class":class_name,
           "source_GBIF_photo_record_count":int(n)}
          for (island,year,y,x,fam,class_name),n in sorted(cellcount.items())]
    tbl=pd.DataFrame(rows,columns=["island","year","latbin","lonbin",
         "source_taxon_family","verified_solar_photo_class",
         "source_GBIF_photo_record_count"])
    require(not tbl.duplicated(["island","year","latbin","lonbin",
       "source_taxon_family","verified_solar_photo_class"]).any(),
       "A true GSHHG reptile photo source-year/time/category was duplicated")
    path=out/"gecko_skink_real_202608_original_island_solar_daynight_photo_cells.csv.gz"
    tbl.to_csv(path,index=False,compression="gzip")
    stats={fam:dict(counts) for fam,counts in family.items()}
    gecko_total=sum(sum(r.values()) for k,r in stats.items() if k in GECKO_FAMILIES)
    skink_total=sum(stats.get("Scincidae",{}).values())
    result={
       "schema":SCHEMA,
       "actual_original_GBIF_DOI":SOURCE_DOI,
       "original_full_Squamata_source_rows":11634053,
       "original_all_three_parquet_source_files":files,
       "real_source_record_quality_gates":dict(tally),
       "gecko_and_skink_original_family_clock_counts":stats,
       "all_gecko_family_island_observations_with_valid_solar_clock":gecko_total,
       "all_scink_island_observations_with_valid_solar_clock":skink_total,
       "observed_islands_with_actual_solar_night_gecko_photo_source":
           int(tbl.loc[tbl.source_taxon_family.isin(GECKO_FAMILIES)&
              tbl.verified_solar_photo_class.eq("solar_night"),"island"].nunique()),
       "observed_islands_with_actual_solar_day_scink_photo_source":
           int(tbl.loc[tbl.source_taxon_family.eq("Scincidae")&
              tbl.verified_solar_photo_class.eq("solar_day"),"island"].nunique()),
       "source_clock_is_verified_UTC_offset_but_original_provider_clock_may_be_wrong":True,
       "solar_day_altitude_gt_plus6_and_civil_night_lt_minus6":True,
       "not_light_at_gecko_body_or_streetlight_association":True,
       "no_direct_gecko_moth_predation_or_nocturnal_pollination":True,
       "zero_time_photographs_are_not_animal_absence":True,
       "source_solar_photo_CSV_SHA256":sha(path)}
    (out/"ALAN_true_original_gecko_Squamata_island_solar_record_receipt.json").write_text(
         json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_TRUE_GECKO_ISLAND_DIEL_OBSERVATION_SOURCE",
       json.dumps({k:result[k] for k in (
        "all_gecko_family_island_observations_with_valid_solar_clock",
        "all_scink_island_observations_with_valid_solar_clock",
        "observed_islands_with_actual_solar_night_gecko_photo_source",
        "observed_islands_with_actual_solar_day_scink_photo_source",
        "gecko_and_skink_original_family_clock_counts")},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--coast",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.manifest,a.coast,a.out)
