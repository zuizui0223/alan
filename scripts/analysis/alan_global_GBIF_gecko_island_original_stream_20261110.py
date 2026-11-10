#!/usr/bin/env python3
"""Complete original 2026-08 GBIF Squamata -> actual island gecko/skink records.

Streaming exact version-fixed original GBIF Parquet source (source DOI
10.15468/dl.8yrbe7) all 11,634,053 original Squamata rows, reusing
SHA-verified GSHHG 2.3.7 islands >=20km2. Restrict HUMAN_OBSERVATION,
PRESENT, 2013–25, valid coordinates. Original GBIF source Gekkonidae,
other explicitly named gecko families and Scincidae (non-gecko lizard
PHOTO-record negative control). No genus is assumed nocturnal, invasive,
moth-eating, pollinating or locally absent without direct evidence.

Outputs ORIGINAL source photographed/georeferenced encounter counts by
island×0.05 degree source grid×true observation year×observed taxon group,
not insect/reptile population abundance. No individual contributor
coordinates, GBIF IDs or observer identities redistributed.
"""
from __future__ import annotations
import argparse,hashlib,json,math,datetime,collections
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from pyarrow.fs import S3FileSystem
import pyarrow.parquet as pq

from scripts.analysis.alan_isolation_area_light_20261108 import load_signed_land

DOI="10.15468/dl.8yrbe7"
PARTITION="taxonpartition=Animalia_Chordata_Reptilia_Squamata"
SCHEMA="ALAN_true_GBIF_202608_Squamata_gecko_skink_GSHHG_island_source_v1"
REGION="gbif-public-data"
GECKO_FAMILIES=("Gekkonidae","Phyllodactylidae","Sphaerodactylidae",
                "Eublepharidae","Diplodactylidae","Carphodactylidae")
CONTROL_FAMILIES=("Scincidae",)
ALL_FAMILIES=GECKO_FAMILIES+CONTROL_FAMILIES
KEY_GENERA=("Hemidactylus","Gehyra","Gekko","Phelsuma","Tarentola",
            "Cyrtodactylus","Lepidodactylus")
ORIGINAL_ROWS=11634053
YEAR_MIN,YEAR_MAX=2013,2025
BIN=.05

def require(ok,reason):
    if not bool(ok):raise ValueError("ALAN_REAL_GECKO_ISLAND_SOURCE: "+str(reason))

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as stream:
        for block in iter(lambda:stream.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def source_files(manifest):
    x=json.loads(Path(manifest).read_text())
    require(x.get("schema")=="ALAN_original_GBIF_2026_08_Gekkota_reptile_Parquet_partition_S3_audit_v1" and
            x.get("original_GBIF_snapshot_DOI")==DOI,
            "Not original frozen 2026-08 public reptile Parquet manifest")
    d=x["matching_reptile_gecko_partition_folders"][PARTITION]
    files=d["original_source_object_metadata"]
    require(len(files)==3 and
            sum(f["compressed_original_bytes"] for f in files)>2*1024**3,
            "Original 3-file Squamata source changed")
    return files

def geo_bbox_landmask(land):
    coarse=np.zeros((180,360),dtype=bool)
    for z in land.geometry.bounds.itertuples(index=False,name=None):
        a,b,c,d=[float(t) for t in z]
        if not np.isfinite([a,b,c,d]).all():continue
        west=max(0,min(359,int(math.floor(a+180))))
        east=max(0,min(359,int(math.floor(c+180))))
        south=max(0,min(179,int(math.floor(b+90))))
        north=max(0,min(179,int(math.floor(d+90))))
        if east>=west and north>=south:
            coarse[south:north+1,west:east+1]=True
    return coarse

def get_col(fields,name):
    items=[x for x in fields if x.casefold()==name.casefold()]
    require(len(items)==1,"Original GBIF Squamata source field missing "+name)
    return items[0]

def run(manifest,coast,out):
    out=Path(out);require(not out.exists(),"Do not overwrite original real gecko source")
    files=source_files(manifest)
    land=load_signed_land(coast)
    land=land.loc[land.land_class.eq("island")&
                  land.area_km2.ge(20.),["land_id","area_km2","geometry"]].copy()
    land["land_id"]=land.land_id.astype(str)
    require(len(land)==3696 and land.land_id.is_unique,
            "Original GSHHG coastal Gekkonidae source differs from real moth island source")
    bbox=geo_bbox_landmask(land)
    s3=S3FileSystem(anonymous=True,region="us-east-1")
    counters=collections.Counter()
    genus_counters=collections.Counter()
    species_cache=collections.defaultdict(set)
    stat=collections.Counter()
    file_receipts=[]
    for i,item in enumerate(files):
        f=s3.open_input_file(REGION+"/"+item["key"])
        pf=pq.ParquetFile(f)
        cols=list(pf.schema_arrow.names)
        wanted=["gbifid","year","family","genus","species",
                "specieskey","basisofrecord","occurrencestatus",
                "decimallatitude","decimallongitude","datasetkey"]
        rename={name:get_col(cols,name) for name in wanted}
        expected_rows=pf.metadata.num_rows
        require(expected_rows>0 and pf.metadata.num_row_groups>0,
                "Actual reptile Parquet source empty")
        file_receipts.append({"source_original_key":item["key"],
                              "original_source_rows":expected_rows,
                              "original_row_groups":pf.metadata.num_row_groups,
                              "original_bytes":item["compressed_original_bytes"]})
        for rg in range(pf.metadata.num_row_groups):
            source=pf.read_row_group(rg,columns=list(rename.values())).to_pandas()
            stat["all_original_reptile_rows_read"]+=len(source)
            stat["real_source_row_groups_read"]+=1
            family=source[rename["family"]].fillna("").astype(str)
            year=pd.to_numeric(source[rename["year"]],errors="coerce")
            ok=(family.isin(ALL_FAMILIES).to_numpy()&
                year.between(YEAR_MIN,YEAR_MAX).to_numpy()&
                source[rename["basisofrecord"]].fillna("").astype(str).str.upper()
                     .eq("HUMAN_OBSERVATION").to_numpy()&
                source[rename["occurrencestatus"]].fillna("").astype(str).str.upper()
                     .eq("PRESENT").to_numpy())
            d=source.loc[ok].copy()
            stat["valid_speciesfamily_human_observation_year_rows"]+=len(d)
            if d.empty:continue
            latitude=pd.to_numeric(d[rename["decimallatitude"]],errors="coerce").to_numpy()
            longitude=pd.to_numeric(d[rename["decimallongitude"]],errors="coerce").to_numpy()
            good=np.isfinite(latitude)&np.isfinite(longitude)&(
                latitude>=-90)&(latitude<=90)&(longitude>=-180)&(longitude<=180)
            d=d.loc[good].copy()
            latitude=latitude[good];longitude=longitude[good]
            stat["actual_GBIF_gecko_skink_original_georeferenced_photo_rows"]+=len(d)
            if d.empty:continue
            ix=np.minimum(359,np.maximum(0,np.floor(longitude+180).astype(int)))
            iy=np.minimum(179,np.maximum(0,np.floor(latitude+90).astype(int)))
            nearby=bbox[iy,ix]
            d=d.loc[nearby].copy()
            latitude=latitude[nearby];longitude=longitude[nearby]
            stat["coarse_near_GSHHG_source_islands"]+=len(d)
            if d.empty:continue
            d["year"]=year.loc[d.index].astype(int).to_numpy()
            d["source_taxon_family"]=d[rename["family"]].astype(str)
            d["source_taxon_genus"]=d[rename["genus"]].fillna("").astype(str)
            d["source_taxon_species"]=d[rename["species"]].fillna("").astype(str)
            d["source_specieskey"]=d[rename["specieskey"]].fillna("").astype(str)
            d["latbin"]=np.minimum(3599,np.floor((latitude+90)/BIN).astype(int))
            d["lonbin"]=np.minimum(7199,np.floor((longitude+180)/BIN).astype(int))
            geography=gpd.GeoDataFrame(d[[
                  "year","latbin","lonbin","source_taxon_family",
                  "source_taxon_genus","source_taxon_species","source_specieskey"]],
                  geometry=gpd.points_from_xy(longitude,latitude),crs="EPSG:4326")
            matched=gpd.sjoin(geography,land[["land_id","geometry"]],
                              how="inner",predicate="within")
            ambiguous=matched.index.duplicated(keep=False)
            stat["ambiguous_island_overlap_original_source_rows"]+=int(ambiguous.sum())
            matched=matched.loc[~ambiguous].copy()
            matched["island"]=matched.land_id.astype(str)
            stat["actual_gecko_or_skink_original_island_observation_rows"]+=len(matched)
            if matched.empty:continue
            group=["island","year","latbin","lonbin","source_taxon_family"]
            n=matched.groupby(group,sort=False).size()
            for key,value in n.items():counters[key]+=int(value)
            gg=matched.loc[matched.source_taxon_genus.isin(KEY_GENERA)].groupby(
                  ["island","year","latbin","lonbin","source_taxon_genus"],sort=False).size()
            for key,value in gg.items():genus_counters[key]+=int(value)
            p=matched.loc[matched.source_specieskey.ne("")].drop_duplicates(
                group+["source_specieskey"])
            for key,g in p.groupby(group,sort=False):
                species_cache[key].update(g.source_specieskey.tolist())
            if stat["real_source_row_groups_read"]%25==0:
                print("ALAN_GLOBAL_GECKO_ACTUAL_GBIF_SOURCE_STREAM",
                   json.dumps({"original_reptile_rows_read":
                       stat["all_original_reptile_rows_read"],
                       "actual_island_gecko_skink_source_records":
                       stat["actual_gecko_or_skink_original_island_observation_rows"],
                       "island_gridyearfamily_rows":len(counters)},
                       sort_keys=True),flush=True)
        f.close()
    require(stat["all_original_reptile_rows_read"]==ORIGINAL_ROWS and
            sum(x["original_source_rows"] for x in file_receipts)==ORIGINAL_ROWS,
            "Original reptile GBIF snapshot was not completely read")
    rows=[]
    for (island,year,latbin,lonbin,family),count in sorted(counters.items()):
        rows.append({"island":island,"year":int(year),
                    "latbin":int(latbin),"lonbin":int(lonbin),
                    "actual_source_taxon_family":family,
                    "actual_GBIF_human_observation_records":int(count),
                    "actual_distinct_species_keys":len(
                        species_cache[(island,year,latbin,lonbin,family)])})
    frame=pd.DataFrame(rows,columns=[
        "island","year","latbin","lonbin","actual_source_taxon_family",
        "actual_GBIF_human_observation_records","actual_distinct_species_keys"])
    require(not frame.duplicated(["island","year","latbin","lonbin",
                                  "actual_source_taxon_family"]).any(),
            "Source reptile family island year key repeated")
    genus=pd.DataFrame([{"island":island,"year":int(year),
           "latbin":int(latbin),"lonbin":int(lonbin),
           "original_genus":genus,"original_source_records":int(value)}
           for (island,year,latbin,lonbin,genus),value in sorted(genus_counters.items())])
    out.mkdir(parents=True)
    csv=out/"GBIF_202608_actual_gecko_skink_island_geo_cells2013_2025.csv.gz"
    genus_csv=out/"GBIF_202608_actual_gecko_genera_island_geo_cells2013_2025.csv.gz"
    frame.to_csv(csv,index=False,compression="gzip")
    genus.to_csv(genus_csv,index=False,compression="gzip")
    photo_gecko=frame.loc[frame.actual_source_taxon_family.isin(GECKO_FAMILIES)]
    photo_skinks=frame.loc[frame.actual_source_taxon_family.eq("Scincidae")]
    report={
        "schema":SCHEMA,
        "source_GBIF_frozen_snapshot_DOI":DOI,
        "full_original_Squamata_all_year_and_basis_rows":ORIGINAL_ROWS,
        "original_source_reptile_Parquet_files":file_receipts,
        "original_Squamata_row_group_source_audit":dict(stat),
        "source_true_GSHHG_20plus_km2_islands":len(land),
        "source_Gekkonidae_family_islands":int(
             frame.loc[frame.actual_source_taxon_family.eq("Gekkonidae"),"island"].nunique()),
        "source_any_named_gecko_families_islands":int(photo_gecko.island.nunique()),
        "source_any_named_gecko_families_island_cell_years":int(
             photo_gecko[["island","year","latbin","lonbin"]].drop_duplicates().shape[0]),
        "source_gecko_photo_observation_records":int(
             photo_gecko.actual_GBIF_human_observation_records.sum()),
        "source_Scincidae_negative_control_islands":int(photo_skinks.island.nunique()),
        "source_Scincidae_photo_observation_records":int(
             photo_skinks.actual_GBIF_human_observation_records.sum()),
        "original_gecko_families":list(GECKO_FAMILIES),
        "excluded_from_main_gecko_family_proxy":["Pygopodidae"],
        "original_gecko_named_genera":list(KEY_GENERA),
        "genus_source_record_summary":{str(k):int(v) for k,v in
            genus.groupby("original_genus").original_source_records.sum().items()}
            if len(genus) else {},
        "source_gecko_family_cell_year_CSV_SHA256":sha256(csv),
        "source_gecko_genus_cell_year_CSV_SHA256":sha256(genus_csv),
        "no_observer_ID_or_original_precise_coordinates_redistributed":True,
        "source_documented_gecko_photos_not_occupancy_or_diet":True,
        "source_other_reptile_or_gecko_zero_is_not_confirmed_real_absence":True,
        "not_direct_ALAN_gecko_moth_interaction_traction":True,
        "no_effective_nocturnal_pollination_or_insect_population_measured":True,
        "source_frozen_snapshot_not_same_as_LIVE_2026_GBIF_API":True}
    (out/"ALAN_2026_11_true_global_gecko_family_source_receipt.json").write_text(
         json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_TRUE_GLOBAL_GECKO_ISLAND_ORIGINAL_SOURCE",
          json.dumps({k:report[k] for k in (
              "source_Gekkonidae_family_islands",
              "source_any_named_gecko_families_islands",
              "source_any_named_gecko_families_island_cell_years",
              "source_gecko_photo_observation_records",
              "source_Scincidae_negative_control_islands",
              "source_Scincidae_photo_observation_records",
              "genus_source_record_summary")},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-manifest",type=Path,required=True)
    p.add_argument("--coast",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.source_manifest,a.coast,a.out)
