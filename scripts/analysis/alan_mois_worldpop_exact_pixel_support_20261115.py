#!/usr/bin/env python3
"""Audit exact WorldPop Global2 R2024B UC same-year pixel reuse at ALAN sites.

MOIS original 0.05-degree photo-cell source and ALAN original 0.05-degree
photo-cell source can overlap, but one 0.05-degree cell has ~36 1km
WorldPop pixels. Never join 1km population solely by a 0.05-degree bin.

Require original year, WGS84 1km pixel index and valid UC population status.
Report aggregate joins and attrition only, not species/coordinates.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import pandas as pd
import numpy as np
from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    LIGHT,load_sites,require
)
YEARS=(2019,2021,2023,2024,2025)
WP_2024_SHA="a07af611836d215fb7997c7f574729df10787f5a7cc9bfa1b67bf5a831d613dc"

def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(8<<20),b""):h.update(b)
    return h.hexdigest()

def bins(lon,lat):
    lon=np.asarray(pd.to_numeric(lon,errors="coerce"),dtype=float)
    lat=np.asarray(pd.to_numeric(lat,errors="coerce"),dtype=float)
    # Exact same WGS84 grid as MOIS original photographed-source cell creator.
    return (np.floor((np.minimum(lon,np.nextafter(180.,-np.inf))+180.)/.05+1e-8).astype(int),
            np.floor((np.minimum(lat,np.nextafter(90.,-np.inf))+90.)/.05+1e-8).astype(int))

def pixels(lon,lat):
    lon=np.asarray(pd.to_numeric(lon,errors="coerce"),dtype=float)
    lat=np.asarray(pd.to_numeric(lat,errors="coerce"),dtype=float)
    # Official 30-arcsecond WGS84 Global2 mosaics 43200×17280, bounds lat [-54,90].
    # Actual latitude index must start at 90° north; source ds height is 17280.
    return (np.floor((np.minimum(lon,np.nextafter(180.,-np.inf))+180)*120+1e-8).astype(int),
            np.floor((90-np.maximum(lat,np.nextafter(-54.,np.inf)))*120+1e-8).astype(int))

def inputs(root,year):
    if year in (2019,2021,2023):
        base=Path(root)/"mois_sites"/f"mois-original-{year}-global-mainland-island-monthly-VIIRS-source-SHA"
        coords=base/"source_year"/str(year)/f"mois_{year}_knownland_sites.csv"
        pop=Path(root)/f"year{year}"/f"2024B_UC_{year}_worldpop_sameyear_photo_landmass_popdensity.csv"
        source_sha={2019:"01f4e2d2e1c8bc9a1392f4143f693563bfa8ea1a0f65ba800e1caa450803831c",
        2021:"f0d7cb134c09c74cccfb69199396ac30bc5b663d39dc5369d6e296d8c3a1cdca",
        2023:"e2825386f45a77da6e3dacdf714fc940cf9ab9cc5174671672b4d57466d04fba"}[year]
        require(digest(coords)==source_sha,"MOIS exact annual source digest changed")
        z=pd.read_csv(coords,usecols=["target_id","query_longitude","query_latitude"],dtype={"target_id":str})
        q=pd.read_csv(pop,usecols=["target_id","population_status","population_density_per_km2"],
                      dtype={"target_id":str})
        q=q.rename(columns={"population_status":"pop_status",
                            "population_density_per_km2":"pop_density"})
    elif year==2024:
        coords=Path(root)/"mois2024/result_rest2024/all2024_known_land_combined_values.csv"
        pop=Path(root)/"mois2024pop/2024_66k_UC_R2024B_fullcohort_values.csv"
        require(digest(coords)=="6e0a077d496e90848d05d043bffae57e4772519f4572a457221c5ef4da0fb813",
                "2024 original MOIS photo coordinates changed")
        require(digest(pop)==WP_2024_SHA,"2024 original MOIS UC population output changed")
        z=pd.read_csv(coords,usecols=["target_id","query_longitude","query_latitude"],dtype={"target_id":str})
        q=pd.read_csv(pop,usecols=["target_id","population_UC_R2024B_status",
                         "population_UC_R2024B_density_per_km2"],dtype={"target_id":str})
        q=q.rename(columns={"population_UC_R2024B_status":"pop_status",
                            "population_UC_R2024B_density_per_km2":"pop_density"})
    elif year==2025:
        coords=Path(root)/"mois2025/numeric_2025/values.csv"
        z=pd.read_csv(coords,usecols=["target_id","query_longitude","query_latitude",
                                      "population_source_status","population_density"],
                      dtype={"target_id":str})
        q=z[["target_id","population_source_status","population_density"]].copy()
        z=z[["target_id","query_longitude","query_latitude"]]
        q=q.rename(columns={"population_source_status":"pop_status",
                            "population_density":"pop_density"})
    else: raise ValueError(year)
    require(z.target_id.is_unique and q.target_id.is_unique,"Duplicate genuine source target ID")
    require(z.target_id.str.startswith(f"{year}|").all() and
            q.target_id.str.startswith(f"{year}|").all(),"Mixed population/photo years")
    m=z.merge(q,on="target_id",how="inner",validate="one_to_one")
    require(len(m)==len(z)==len(q),"Original photo versus population source join lost targets")
    require(m.loc[m.pop_status.ne("OK"),"pop_density"].isna().all(),
            "Missing population silently imputed to zero")
    require(m.loc[m.pop_status.eq("OK"),"pop_density"].notna().all() and
            m.loc[m.pop_status.eq("OK"),"pop_density"].ge(0).all(),
            "Bad valid WorldPop input")
    m["lonbin"],m["latbin"]=bins(m.query_longitude,m.query_latitude)
    m["pop_col"],m["pop_row"]=pixels(m.query_longitude,m.query_latitude)
    m["year"]=year
    return m[["year","lonbin","latbin","pop_col","pop_row","pop_status","pop_density","target_id"]]

def join(sites,root,*,return_data=False):
    require(set(sites.year.unique())==set(YEARS),"Wrong source years")
    source=[]
    for y in YEARS:
        m=inputs(root,y)
        source.append(m)
    vals=pd.concat(source,ignore_index=True)
    require(vals.target_id.is_unique,"Repeated MOIS target across years")
    a=sites.copy()
    a["site_id"]=np.arange(len(a))
    require(a.site_id.is_unique,"Source site IDs must stay distinct even across overlapping 0.05-degree multi-island cells")
    a["pop_col"],a["pop_row"]=pixels(a.lon,a.lat)
    common=["year","lonbin","latbin"]
    matched=a[["site_id","island","year","lonbin","latbin","pop_col","pop_row",LIGHT]].merge(
        vals,on=common,how="left",validate="many_to_many",suffixes=("_alan","_mois"))
    matched["photo_bin_match"]=matched.target_id.notna()
    matched["same_pop_pixel"]=(matched.photo_bin_match &
        matched.pop_col_alan.eq(matched.pop_col_mois) &
        matched.pop_row_alan.eq(matched.pop_row_mois))
    matched["usable"]=(matched.same_pop_pixel & matched.pop_status.eq("OK") &
                       matched.pop_density.notna())
    # Exact matching 1km pixel should have one unique population estimate.
    q=matched.loc[matched.usable].copy()
    conflict=q.groupby("site_id").pop_density.agg(["count","min","max"])
    require((conflict["max"]-conflict["min"]).abs().lt(1e-8).all(),
            "Different population densities assigned to same exact year-pixel site")
    accepted=q.drop_duplicates("site_id").copy()
    require(accepted.site_id.is_unique,"Duplicate eligible ALAN photo site")
    original_year=a.groupby("year").size().to_dict()
    per_year=[]
    for y in YEARS:
        z=matched.loc[matched.year.eq(y)]
        original=int(original_year[y])
        groups=z.groupby("site_id",observed=True)
        photo=int(groups.photo_bin_match.any().sum())
        samepixel=int(groups.same_pop_pixel.any().sum())
        eligible=int(groups.usable.any().sum())
        per_year.append({"year":y,"original_ALAN_site_years":original,
                         "same_005_degree_MOIS_photo_bin":photo,
                         "same_original_WorldPop_1km_pixel":samepixel,
                         "valid_sameyear_R2024B_UC_population":eligible,
                         "valid_fraction":round(eligible/original,6),
                         "unique_MOIS_photo_source_target_bins":int(vals.loc[vals.year.eq(y),
                            ["lonbin","latbin"]].drop_duplicates().shape[0])})
    relevant=a.merge(accepted[["site_id","pop_density"]],on="site_id",how="left",
                     validate="one_to_one")
    require(relevant.pop_density.notna().sum()==len(accepted),"Joined source population lost")
    # Smallish annual support may be a photograph-source selection issue.
    groups=relevant.loc[relevant.pop_density.notna()].groupby("island",observed=True)
    summ={
      "schema":"ALAN_MOIS_WorldPop_exact_pixel_overlap_audit_20261115_v1",
      "status":"ORIGINAL_YEAR_AND_ORIGINAL_1KM_PIXEL_MATCHED_COVERAGE_ONLY",
      "source_population":"MOIS original annual 2019/2021/2023/2024/2025 WorldPop Global2 R2024B UC 1km",
      "original_ALAN_island_photo_sites":len(sites),
      "ALAN_original_years":[int(y) for y in YEARS],
      "matched_usable_site_years":len(accepted),
      "matched_distinct_islands":int(accepted.island.nunique()),
      "source_years":per_year,
      "unmatched_source_site_years":len(sites)-len(accepted),
      "source_zero_population_site_years":int(relevant.pop_density.eq(0).sum()),
      "original_population_no_imputed_zero":True,
      "same_year_required":True,
      "same_worldpop_pixel_required":True,
      "bin_only_nonpixel_matching_rejected":True,
      "no_ALAN_light_response_or_pop_control_fitted":True,
      "no_source_species_or_sensitive_coordinates_exported":True,
      "matched_population_quantiles":relevant.pop_density.dropna().quantile([0,.25,.5,.75,1]).to_dict(),
    }
    return (summ,relevant) if return_data else summ

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites",type=Path,required=True)
    p.add_argument("--sources",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    require(not a.out.exists(),"Refuse original audit overwrite")
    s=load_sites(a.sites)
    require(len(s)==77647 and s.island.nunique()==562,"Original 562-island source changed")
    q=join(s,a.sources)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(q,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_AUTHENTIC_SAMEYEAR_WORLDPOP_PIXEL_MATCH",json.dumps(q,allow_nan=False),flush=True)
if __name__=="__main__":main()
