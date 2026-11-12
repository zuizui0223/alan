"""Test positive predator-record overlap at seasonal monitoring cells. No absence inference."""
import argparse
import gzip
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

from alan_repeat_flower_contrasts import GROUP, UNIT, extract_units, four, read_source
from mois_global_2020_GSHHG_scope_20261107 import assign_points, prepare_land

NAMES = {9537238:"Hemidactylus frenatus",5959942:"Hemidactylus mabouia",
         1311477:"Vespa velutina",5170994:"Argiope bruennichi",6258028:"Mantis religiosa"}


def run(v, predator, coast, cache, out):
    cache.mkdir(parents=True,exist_ok=True)
    supported = four(v).reset_index()[GROUP]
    supported.columns = GROUP
    subset = v.merge(supported, on=GROUP)
    loc = subset[["key","land_id"]].drop_duplicates()
    records = []
    for year in sorted(v.year.unique()):
        for taxon,name in NAMES.items():
            directory = "private" if taxon in (9537238,5959942) else "multiguild_private"
            with gzip.open(predator/directory/f"{year}_{taxon}.json.gz","rt") as f:
                page=json.load(f)
            assert len(page["rows"]) == page["count"]
            for r in page["rows"]:
                if (r["speciesKey"] != taxon or r["species"] != name or r["year"] != year or
                    r["basisOfRecord"] != "HUMAN_OBSERVATION" or r["occurrenceStatus"] not in ["PRESENT","DETECTED"] or
                    not 0 < (r.get("coordinateUncertaintyInMeters") or 0) <= 1000):
                    continue
                gx = int(np.floor((r["decimalLongitude"]+180)*20))
                gy = int(np.floor((r["decimalLatitude"]+90)*20))
                records.append((r["key"],f"{year}|{gx}_{gy}",r["decimalLongitude"],r["decimalLatitude"],name))
    det = pd.DataFrame(records,columns=["unit_id","key","longitude","latitude","predator"]).drop_duplicates("unit_id")
    det = det[det.key.isin(loc.key)].copy()
    if len(det):
        land=assign_points(det,prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp")),gpd.read_file(coast/"GSHHS_h_L2.shp"))
        det=det.merge(land[["unit_id","land_id"]],on="unit_id").merge(loc,on=["key","land_id"])
    else:
        det["land_id"] = pd.Series(dtype=str)

    def day(uid):
        file=cache/f"{uid}.json"
        for attempt in range(3):
            try:
                if file.exists():
                    rec=json.loads(file.read_bytes())
                else:
                    res=requests.get(f"https://api.gbif.org/v1/occurrence/{uid}",timeout=15)
                    res.raise_for_status()
                    rec=res.json()
                    if rec.get("key") != uid:
                        return None
                    file.write_text(json.dumps({k:rec.get(k) for k in ["key","year","month","day","eventDate"]}))
                if rec.get("day") and rec.get("month") and "/" not in str(rec.get("eventDate")):
                    return pd.Timestamp(int(rec["year"]),int(rec["month"]),int(rec["day"]))
                return None
            except (ValueError,requests.RequestException):
                if attempt < 2:
                    time.sleep(2 ** attempt)
        return None
    with ThreadPoolExecutor(max_workers=6) as pool:
        det["predator_date"] = pd.to_datetime(list(pool.map(day, det.unit_id.tolist())))
    x=subset[subset.status.eq(1)].drop(columns=["unit_id"],errors="ignore").merge(det[["key","land_id","unit_id","predator","predator_date"]],on=["key","land_id"])
    x["date"]=pd.to_datetime(x.date)
    complete = extract_units(subset).dropna(subset=["lo","hi"])
    bouts = x.merge(complete[UNIT+["first","last"]],on=UNIT)
    in_bout = bouts[bouts.predator_date.dt.dayofyear.between(bouts["first"],bouts["last"])]
    results=[]
    for scale,group in [("year",x),("within_bracketed_open_bout",in_bout),("same_open_observation_date",x[x.date==x.predator_date])]:
        for ix,g in group.groupby(["land_class","light","predator"]):
            results.append(dict(zip(["land_class","light","predator"],ix)) | dict(scale=scale,cell_years=g.key.nunique(),
                       taxa=g.plant.nunique(),species_region_years=g[GROUP].drop_duplicates().shape[0],animal_records=g.unit_id.nunique()))
    pd.DataFrame(results,columns=["land_class","light","predator","scale","cell_years","taxa","species_region_years","animal_records"]).to_csv(out/"predator_overlap.csv",index=False)
    (out/"predator_support.json").write_text(json.dumps(dict(eligible_candidate_cells=len(loc),qualifying_animal_records=len(det),
              animal_records_with_verified_day=int(det.predator_date.notna().sum()),
              positive_cell_years=x.key.nunique(),positive_bounded_open_bout_cell_years=in_bout.key.nunique(),positive_same_open_day_cell_years=x.loc[x.date==x.predator_date,"key"].nunique(),
              counts_are_positive_records_not_rates=True),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",action="append",nargs=3,required=True)
    for k in ["predator","coast","cache","out"]:
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    run(pd.concat([read_source(*r) for r in a.source],ignore_index=True),a.predator,a.coast,a.cache,a.out)
