"""Acquire public Budburst adult pollinator result pages and the UGA public results layer."""
import argparse
import json
from pathlib import Path
from urllib.parse import unquote

import requests

YEARS = (2019, 2021, 2023, 2024, 2025)
BASE = "https://budburst.org"
UGA = "https://services2.arcgis.com/PYn6bWCjT6bhw1z3/arcgis/rest/services/survey123_956c41fc2fa34638bb3342ac06905588_results/FeatureServer/0/query"


def acquire(folder):
    folder.mkdir(parents=True,exist_ok=True)
    session = requests.Session()
    session.get(BASE+"/data",timeout=45).raise_for_status()
    for year in YEARS:
        rows = []
        page = 1
        while True:
            filters = dict(page=page,report_types=["pollinators"],plant_groups=[],state_id=[],
                           from_date=f"{year}-01-01",to_date=f"{year}-12-31",
                           include_youth_observations=False,budburst_species_only=False)
            headers = {"X-Requested-With":"XMLHttpRequest","Accept":"application/json","Referer":BASE+"/data",
                       "X-XSRF-TOKEN":unquote(session.cookies.get("XSRF-TOKEN",""))}
            response = session.post(BASE+"/data/filter",json=filters,headers=headers,timeout=45)
            response.raise_for_status()
            data = response.json()
            for item in data["data"]:
                site = item["site_species"]
                rows.append({key:item[key] for key in ["id","observation_date","species_scientific_name","report_type_id"]}
                            | {"site_species":{key:site.get(key) for key in ["longitude","latitude"]}})
            if page==data["meta"]["last_page"]:
                break
            page += 1
        assert len(rows)==data["meta"]["total"]==len({r["id"] for r in rows})
        (folder/f"bud_{year}.json").write_text(json.dumps(rows),encoding="utf8")
    fields = ["objectid","date_month_day_year","start_time","end_time","scientific_plant_name"]
    fields += ["how_many_"+name for name in ("carpenter_bees","bumble_bees","honey_bees","small_bees","wasps","flies")]
    query = dict(where="1=1",outFields=",".join(fields),f="geojson",outSR=4326,resultRecordCount=2000)
    r = requests.get(UGA,params=query,timeout=45)
    r.raise_for_status()
    expected = requests.get(UGA,params=dict(where="1=1",f="json",returnCountOnly="true"),timeout=45).json()["count"]
    assert len(r.json()["features"])==expected
    (folder/"uga_times.json").write_bytes(r.content)


if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder",type=Path)
    acquire(ap.parse_args().folder)
