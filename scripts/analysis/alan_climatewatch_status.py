"""Fetch public project categories and dated flowering visits from the ClimateWatch map."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from alan_public_phenology_world import OUT, YEARS, wcvp_species

LANG = "https://www.spotteron.com/assets/map/config/climatewatch/lang/en-AU.json?v=124"
SPOTS = "https://www.spotteron.com/api/v2.6/spots"


def source_rows(folder, category):
    cache = folder / f"climatewatch_{category}"
    cache.mkdir(exist_ok=True)
    params = {"filter[topic_id]":40,"filter[category]":category,"limit":100}
    records = []
    page = 1
    while True:
        path = cache/f"p{page:03d}.json"
        if not path.exists():
            resp = requests.get(SPOTS,params={**params,"page":page},timeout=30)
            resp.raise_for_status()
            path.write_text(json.dumps(resp.json()))
        batch = json.loads(path.read_text())
        records.extend(batch["data"])
        if page >= batch["meta"]["page_count"]:
            break
        page += 1
    assert len({item["id"] for item in records}) == len(records) == batch["meta"]["total"]
    assert all(i["attributes"]["topic_id"]==40 and i["attributes"]["category"]==category for i in records)
    return records


def normalized(name):
    return " ".join(str(name or "").strip().lower().split()).capitalize()


def fetch(folder, wcvp):
    resp = requests.get(LANG, timeout=30)
    resp.raise_for_status()
    labels = resp.json()
    categories = []
    for key, value in labels.items():
        match = re.fullmatch(r"CONFIG_FIELD_764_OPTION_(\d+)_META", key)
        if match:
            categories.append((int(match[1]), value))
    categories = pd.DataFrame(categories, columns=["category","scientific_name"])
    other = source_rows(folder,3948)
    fields = {i["id"]:normalized(i["attributes"].get("fld_01_00000769")) for i in other}
    categories["scientific_name"] = categories.scientific_name.map(normalized)
    names = wcvp_species(wcvp, set(categories.scientific_name)|set(fields.values()))
    cohort = pd.read_parquet(folder/"target.parquet")[["plant","wcvp_accepted_plant_name_id"]].drop_duplicates()
    cohort["accepted"] = cohort.wcvp_accepted_plant_name_id.astype(str)
    resolved = names.merge(cohort,on="accepted")
    focal = categories.merge(resolved,left_on="scientific_name",right_on="taxon_name")
    lookup = dict(zip(resolved.taxon_name,resolved.plant))
    candidates = [(i,lookup[fields[i["id"]]]) for i in other if fields[i["id"]] in lookup]
    other_shared = len(candidates)
    totals = {}
    for r in focal.itertuples():
        records = source_rows(folder,r.category)
        totals[r.plant] = len(records)
        candidates.extend((i,r.plant) for i in records)
    visits = []
    for item,plant in candidates:
        obs = item["attributes"]
        flags = set(obs["fld_04_00000772"] or [])
        pos = bool(flags & {3978,3979})
        neg = bool(flags & {3981,4821})
        if pos == neg or 3982 in flags or (3980 in flags and neg):
            continue
        date = obs["spotted_at"][:10]
        if int(date[:4]) in YEARS and obs.get("latitude") is not None and obs.get("longitude") is not None:
            visits.append((plant,item["id"],obs.get("root_id") or item["id"],date,int(pos),
                           obs["latitude"],obs["longitude"]))
    x = pd.DataFrame(visits,columns=["plant","obs_id","station","date","status","latitude","longitude"])
    x["year"] = x.date.str[:4].astype(int)
    x["region5"] = (np.floor((x.longitude+180)/5).astype(int).astype(str)+"_"+
                    np.floor((x.latitude+90)/5).astype(int).astype(str))
    x["key"] = (x.year.astype(str)+"|"+np.floor((x.longitude+180)*20).astype(int).astype(str)+"_"+
                np.floor((x.latitude+90)*20).astype(int).astype(str))
    cells = pd.read_parquet(folder/"universe.parquet")[["region5","year","key","land_class","light"]].drop_duplicates()
    reg = x.merge(cells[["region5","year"]].drop_duplicates(),on=["region5","year"],validate="m:1")
    reg["individual"] = None
    reg["sampling_unit"] = "spot_update_thread"
    reg.to_parquet(folder/"australia_standard_visits.parquet",index=False)
    exact = reg.merge(cells,on=["key","region5","year"],validate="m:1")
    summary = dict(named_source_categories=len(categories),record_totals=totals, other_plant_rows=len(other),other_plant_matched=other_shared, focal_year_status_rows=len(x),
                   region_year_status_rows=len(reg),region_year_positive=int(reg.status.sum()),
                   region_year_negative=int(reg.status.eq(0).sum()),distinct_target_regions=reg.region5.nunique(),
                   distinct_shared_plants=reg.plant.nunique(),distinct_spot_threads=reg.station.nunique(),
                   same_frozen_cell_year_rows=len(exact))
    (OUT/"climatewatch_live_screen.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(summary)


if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder",type=Path)
    ap.add_argument("wcvp",type=Path)
    args = ap.parse_args()
    fetch(args.folder,args.wcvp)
