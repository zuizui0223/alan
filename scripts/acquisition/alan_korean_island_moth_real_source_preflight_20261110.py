#!/usr/bin/env python3
"""Inspect genuine Korean island moth light-trap abundance + site metadata.

Independent of ALAN iNaturalist flora source. Retain only author-file
MD5-verified aggregate schema/coverage counts; do NOT copy original raw
species-by-site matrix or coordinates into git repository or artifacts.
"""
from __future__ import annotations
import argparse, hashlib, io, json, pathlib, urllib.request
import numpy as np
import pandas as pd

RECORD="https://zenodo.org/records/21438943"
SOURCES={
 "occurrences":("data1-islands_species_occurrences.xlsx",
               "1bb467aadecda9ecca0c45460a7f1646"),
 "taxa":("data2-species_metadata.csv","c787ba7e8ca474be9bdf6a901c3e6de2"),
 "sites":("data3-islands_site_metadata.csv","552b24db19a42c5e08a260a7da8d696f"),
 "islands":("data4-island_statistics.xlsx","aadc4eefc6d479d4f13ce80f598553ed"),
 "readme":("Readme.md","f5494f8387a896ea7ed0a6d0246faed8"),
 "metadata":("metadata.csv","b871e9aeb860d6db09b85abdd4fd5291")
}
SCHEMA="ALAN_Korean_island_moth_genuine_2007_2025_original_source_preflight_v1"

def fail(v,msg):
    if not v:raise ValueError("ALAN_REAL_KOREAN_MOTHS: "+msg)

def get_original(name,md5):
    req=urllib.request.Request(f"{RECORD}/files/{name}?download=1",headers={
        "User-Agent":"ALAN research original Korean insular moth source audit/1.0"})
    with urllib.request.urlopen(req,timeout=120) as r:
        fail(r.status==200,"Source unauthorized or missing for "+name)
        body=r.read()
    fail(len(body)>50, "Empty original source "+name)
    got=hashlib.md5(body).hexdigest()
    fail(got==md5,f"Original published MD5 wrong for {name}: {got}")
    return body

def inspect_xlsx(raw):
    result=[]
    with pd.ExcelFile(io.BytesIO(raw)) as excel:
        for sheet in excel.sheet_names:
            df=pd.read_excel(excel,sheet_name=sheet,header=0)
            cols=[str(x) for x in df.columns]
            result.append({
                "sheet":sheet,"n_rows":len(df),"n_columns":len(df.columns),
                "first_columns":cols[:35],
                "first_rows_safe":df.iloc[:4,:min(len(cols),8)].fillna("NA").astype(str).values.tolist(),
                "numeric_columns_count":int(df.select_dtypes(include="number").shape[1]),
                "explicit_year_column_candidates":[c for c in cols if "year" in c.lower()],
                "explicit_coordinate_candidates":[c for c in cols
                    if any(z in c.lower() for z in ["lat","lon","coord","geog","utm"])],
            })
    return result

def inspect_csv(raw):
    try:
        df=pd.read_csv(io.BytesIO(raw),encoding="utf-8-sig",low_memory=False)
    except UnicodeDecodeError:
        # Original authors explicitly document CP949 for Korean taxonomy.
        df=pd.read_csv(io.BytesIO(raw),encoding="cp949",low_memory=False)
    cols=[str(x) for x in df.columns]
    return {
        "n_rows":len(df),"n_columns":len(df.columns),
        "all_columns":cols[:75],
        "first_rows_safe":df.iloc[:5,:min(len(cols),10)].fillna("NA").astype(str).values.tolist(),
        "numeric_columns_count":int(df.select_dtypes(include="number").shape[1]),
        "missing_fraction_by_column":{str(k):round(float(v),3)
              for k,v in df.isna().mean().items()},
        "candidate_site_keys":[c for c in cols if
                any(v in c.lower() for v in ["site","station","point","sample","id"])],
        "candidate_dates":[c for c in cols if
                any(v in c.lower() for v in ["year","date","month","day","time"])],
        "candidate_coordinates":[c for c in cols if
                any(v in c.lower() for v in ["lat","lon","coord","utm"])],
    }

def run(out):
    outcome={"schema":SCHEMA,"zenodo_record":RECORD,
        "source_DOI":"10.5281/zenodo.21438943",
        "not_iNaturalist_plant_photo_data":True,
        "not_ground_lux_nor_flower_pollen_deposition":True,
        "genuine_original_moth_capture_dataset":True,
        "retains_original_source_file_binaries":False,
        "sources":{}}
    for key,(name,expected_md5) in SOURCES.items():
        try:
            raw=get_original(name,expected_md5)
            rec={"filename":name,"published_MD5_verified":expected_md5,
                 "downloaded_original_bytes":len(raw),
                 "original_sha256":hashlib.sha256(raw).hexdigest(),
                 "source_access":"original_bytes_verified"}
            if name.endswith(".xlsx"):
                rec["sheet_schemas"]=inspect_xlsx(raw)
            elif name.endswith(".csv"):
                rec["table_schema"]=inspect_csv(raw)
            else:
                rec["readme_excerpt"]=raw.decode("utf8","replace")[:6500]
            outcome["sources"][key]=rec
            print("ALAN_TRUE_ISLAND_MOTH_SOURCE",json.dumps({key:rec},sort_keys=True),flush=True)
        except Exception as exc:
            outcome["sources"][key]={"filename":name,
                "source_access":"unverified_original_source",
                "error":str(exc)[:450]}
            print("ALAN_ISLAND_MOTH_SOURCE_ACCESS_GAP",json.dumps(
                {key:outcome["sources"][key]},sort_keys=True),flush=True)
    # Independently reconstruct actual 2007-2025 source opportunity (not a
    # model). Each original row is a site × recorded year, not a moth
    # specimen-level independent sample and not a repeated annual census.
    if all(outcome["sources"].get(k,{}).get("source_access")==
           "original_bytes_verified" for k in ("occurrences","taxa","sites","islands")):
        m=pd.read_excel(io.BytesIO(get_original(*SOURCES["occurrences"])))
        loc=pd.read_csv(io.BytesIO(get_original(*SOURCES["sites"])))
        taxa=pd.read_csv(io.BytesIO(get_original(*SOURCES["taxa"])),encoding="cp949")
        island_meta=pd.read_excel(io.BytesIO(get_original(*SOURCES["islands"])))
        cols=[x for x in m.columns if str(x).startswith("sp")]
        fail(len(m)==377 and len(cols)==1179 and len(loc)==332 and
             len(island_meta)==30 and m.year.between(2007,2025).all(),
             "Original author-file row/column/year contract failed")
        fail(not m.duplicated(["island","year","sites"]).any(),
             "True source site-year is duplicated")
        fail(not loc.duplicated(["island","site_id"]).any(),
             "Original location was duplicated")
        fail(set(cols)==set(taxa["Species id"].astype(str)),
             "Original taxonomic species ID source does not reconcile")
        counts=m[cols].apply(pd.to_numeric,errors="coerce")
        fail(counts.notna().all().all() and (counts>=0).all().all() and
             np.all(counts.to_numpy()==np.floor(counts.to_numpy())),
             "Original moth data are not complete integer capture matrix")
        m["moth_individuals"]=counts.sum(axis=1)
        m["moth_recorded_species"]=(counts>0).sum(axis=1)
        v=m.merge(loc,left_on=["island","sites"],
                  right_on=["island","site_id"],how="left",
                  validate="many_to_one",indicator=True)
        missing=v._merge.ne("both")
        nsite_year=m.groupby("year",sort=True).agg(
            site_years=("sites","size"),
            islands=("island","nunique"),
            moth_captures=("moth_individuals","sum"),
            recorded_species_site_sum=("moth_recorded_species","sum"))
        il_years=m.groupby("island").year.nunique()
        x2025=m.loc[m.year.eq(2025)]
        eligible=m.loc[m.year.ge(2012)]
        outcome["actual_original_capture_support"]={
            "original_source_rows":len(m),
            "original_source_site_coordinates":len(loc),
            "original_source_islands":int(m.island.nunique()),
            "original_independent_site_IDs":int(m[["island","sites"]].drop_duplicates().shape[0]),
            "original_taxonomy_species_count":len(cols),
            "total_individual_moth_counts":int(m.moth_individuals.sum()),
            "species_seen_at_least_once":int((counts.sum(axis=0)>0).sum()),
            "zero_individual_source_site_years":int(m.moth_individuals.eq(0).sum()),
            "site_rows_missing_author_location":int(missing.sum()),
            "source_sites_with_repeated_observation_years":int(
                m.groupby(["island","sites"]).year.nunique().ge(2).sum()),
            "source_islands_with_multiple_years":int(il_years.ge(2).sum()),
            "source_islands_with_3_or_more_years":int(il_years.ge(3).sum()),
            "source_year_counts":[{"year":int(y),**{str(k):int(v)
              for k,v in rec.items()}} for y,rec in nsite_years.to_dict("index").items()],
            "years_2012_or_later":{"site_years":int(len(eligible)),
                "islands":int(eligible.island.nunique())},
            "genuinely_observed_2025":{"site_years":int(len(x2025)),
                "islands":int(x2025.island.nunique()),
                "individuals":int(x2025.moth_individuals.sum()),
                "nonzero_sites":int(x2025.moth_individuals.gt(0).sum())},
            "location_coordinate_schema":{"Latitude": "WGS84 degrees",
                "Longitude":"WGS84 degrees"},
            "has_sampling_month_or_actual_visit_date":False,
            "has_complete_trap_nights_or_days_effort_column":False,
            "light_traps_catchability_is_light_attraction_sensitive":True,
            "no_ecological_ALAN_effect_estimated":True}
        print("ALAN_REAL_KOREA_MOTH_STUDY_YEAR_SITE_SUPPORT",
              json.dumps(outcome["actual_original_capture_support"],
                         sort_keys=True),flush=True)

    dest=pathlib.Path(out)
    fail(not dest.exists(),"Cannot overwrite an existing run's provenance")
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(outcome,indent=2,sort_keys=True,default=str)+"\n")
    print("ALAN_KOREAN_ISLAND_MOTH_PREFLIGHT",json.dumps({
        "source_files_verified":[k for k,v in outcome["sources"].items()
                                 if v["source_access"]=="original_bytes_verified"],
        "source_files_failed":[k for k,v in outcome["sources"].items()
                               if v["source_access"]!="original_bytes_verified"]},
        sort_keys=True),flush=True)
    return outcome

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=pathlib.Path,required=True)
    a=p.parse_args()
    run(a.out)
