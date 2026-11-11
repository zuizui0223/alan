"""Screen US timed visitor reports against the frozen cohort; cache public biological fields only."""
import argparse
import concurrent.futures
import json
import re
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from statsmodels.formula.api import wls
import requests
from bs4 import BeautifulSoup

from alan_shared_visitors_regional import YEARS, grid, taxonomy, assign_points, prepare_land

ROOT = Path(__file__).resolve().parents[2]
GROUPS = {"Small bees & flies", "Honey bees", "Bumble bees", "Large bees & wasps", "Beetles"}


def detail(report_id, cache):
    file = cache / f"{report_id}.json"
    if file.exists():
        return json.loads(file.read_text())
    url = f"https://budburst.org/data/{report_id}"
    for attempt in range(4):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 3:
                raise
            time.sleep(2**attempt)
    page = BeautifulSoup(response.content, "html.parser")
    text = " ".join(page.stripped_strings)
    groups = {}
    for tag in page.select("div.ml-3.text-sm.font-sans-secondary"):
        label = tag.get_text(" ", strip=True)
        value = tag.parent.get_text(" ", strip=True).split(" ", 1)[0]
        if value.isdigit():
            groups[label] = int(value)
    span = re.search(r"\((\d+)\s+minutes?\)", text)
    event = dict(id=report_id, minutes=int(span[1]) if span else None, groups=groups)
    file.write_text(json.dumps(event))
    return event


def screen(p, inputs, out):
    out.mkdir(parents=True, exist_ok=True)
    cache = p / "bud_details"
    cache.mkdir(exist_ok=True)
    frozen = pd.read_parquet(p / "s.parquet")
    keys = ["year", "region5", "plant"]
    local = pd.read_parquet(p / "shared.parquet").rename(columns={"wcvp_accepted_name":"plant"})
    pool = set(pd.read_csv(ROOT/"results/shared_visitors_regional/shared_544.csv").plant)
    assert len(pool)==544 and frozen.region5.nunique()==78 and frozen.region_year.nunique()==195
    rows = []
    for year in sorted(YEARS):
        for r in json.loads((p/f"bud_{year}.json").read_text()):
            site = r["site_species"]
            assert int(r["observation_date"][:4]) == year
            rows.append(dict(source="Budburst", id=r["id"], year=year, plant_raw=r["species_scientific_name"], date=r["observation_date"],
                             x=site.get("longitude"), y=site.get("latitude")))
    usa = json.loads((p/"uga_times.json").read_text(encoding="utf8"))
    for f in usa["features"]:
        a = f["properties"]
        year = pd.to_datetime(a["date_month_day_year"], unit="ms").year
        if year not in YEARS or not f["geometry"]:
            continue
        start, end = [pd.to_datetime(a[k], format="%H:%M", errors="coerce") for k in ("start_time","end_time")]
        x,y = f["geometry"]["coordinates"]
        rows.append(dict(source="UGA", id=a["objectid"], year=year, plant_raw=a["scientific_plant_name"],
                         x=x, y=y, minutes=(end-start).total_seconds()/60))
    d = pd.DataFrame(rows)
    for col in ["x","y"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d = d.dropna(subset=["x","y"]).copy()
    d["region5"] = grid(d.x, d.y)
    d = d.merge(frozen[["year","region5"]].drop_duplicates(),on=["year","region5"])
    # Allow named cultivars and botanical authors after an unambiguous binomial.
    d["binomial"] = d.plant_raw.str.strip().str.extract(r"^([A-Z][a-z]+ [a-z-]+)(?:$| ['‘(]| [A-Z][a-z]*\.)")[0]
    w = taxonomy(inputs,set(d.binomial.dropna()),set())
    labels = w.loc[w.resolution.eq("resolved"),"accepted_name"].to_dict()
    labels.update({n:n for n in pool})
    d["plant"] = d.binomial.map(labels)
    d = d[d.plant.isin(pool)].copy()
    pts = d[["x","y"]].drop_duplicates().rename(columns={"x":"longitude","y":"latitude"}).reset_index(drop=True)
    pts["unit_id"] = pts.index
    coast = inputs/"coast"
    a = assign_points(pts,prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp")),gpd.read_file(coast/"GSHHS_h_L2.shp"))
    d = d.merge(pts.merge(a,on="unit_id").rename(columns={"longitude":"x","latitude":"y"})[["x","y","land_class"]],on=["x","y"])
    d = d.merge(local[keys].drop_duplicates().assign(local_photo=True),on=keys,how="left")
    d["local_photo"] = d.local_photo.eq(True)
    d.to_csv(p/"bud_uga_frame.csv",index=False)
    d = d[d.local_photo & d.land_class.isin(["island","mainland"])].copy()
    print("Candidates",d.groupby(["source","land_class"]).size().to_dict(),flush=True)
    ids = d.loc[d.source.eq("Budburst"),"id"].tolist()
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
        details = list(ex.map(lambda i:detail(i,cache), ids))
    for record in details:
        selected = d.source.eq("Budburst") & d.id.eq(record["id"])
        d.loc[selected,"minutes"] = record["minutes"]
        if GROUPS.issubset(record["groups"]):
            d.loc[selected,"confirmed_nonmoth"] = sum(record["groups"][k] for k in GROUPS)
        if "Bees & other" in record["groups"]:
            d.loc[selected,"broad_nonlep"] = record["groups"]["Bees & other"]
        elif GROUPS.issubset(record["groups"]):
            d.loc[selected,"broad_nonlep"] = sum(record["groups"][k] for k in GROUPS)
        if record["groups"]:
            d.loc[selected,"total_visitors"] = sum(record["groups"].values())
    d.drop(columns=["x","y"]).to_csv(out/"us_candidate_events.csv",index=False)
    rate = d[d.minutes.gt(0) & d.broad_nonlep.notna()].copy()
    keys = ["region5","year"]
    sides = rate.groupby(keys).land_class.nunique()
    paired = rate.merge(sides.loc[sides.eq(2)].reset_index()[keys],on=keys)
    paired.drop(columns=["x","y"]).to_csv(out/"bilateral_events.csv",index=False)
    ledger = frozen.groupby(keys).key.size().rename("primary_site_years").reset_index()
    for side in ["island","mainland"]:
        eligible = rate[rate.land_class.eq(side)].groupby(keys).id.size().rename(side+"_budburst_reports")
        ledger = ledger.merge(eligible,on=keys,how="left")
        ledger[eligible.name] = ledger[eligible.name].fillna(0).astype(int)
    ledger.to_csv(out/"regional_frame.csv",index=False)
    frozen["original"] = 100*frozen.all_moth/frozen.n_species
    frozen["shared"] = 100*frozen.key.map(local.groupby("key").wcvp_accepted_plant_name_id.nunique()).fillna(0)/frozen.n_species
    results = []
    for (region,year), reports in paired.groupby(keys):
        g = frozen[frozen.region5.eq(region) & frozen.year.eq(year)]
        row = dict(region5=region,year=year)
        counts = g.groupby("land_class").key.transform("size")
        for side, r in reports.groupby("land_class"):
            row[side+"_visits"] = int(r.broad_nonlep.sum())
            row[side+"_observer_minutes"] = int(r.minutes.sum())
            row[side+"_per_10min"] = 10*r.broad_nonlep.sum()/r.minutes.sum()
            row[side+"_site_years"] = int(g.land_class.eq(side).sum())
            row[side+"_plants"] = "|".join(sorted(r.plant.unique()))
        ranges = g.groupby("land_class").loglight.agg(["min","max","nunique"])
        if ranges["nunique"].min() >= 2 and ranges["min"].max() < ranges["max"].min():
            for y in ("original","shared"):
                model = wls(y+" ~ island * loglight + logpop + logphotos + logspecies + bio5 + lograin + elevation_km",
                            data=g,weights=counts.min()/counts).fit()
                if model.model.rank == len(model.params):
                    row[y+"_gap"] = model.params["island:loglight"]
        results.append(row)
    pd.DataFrame(results).to_csv(out/"bilateral_results.csv",index=False)
    return d


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    for name in ["survey_cache","inputs","out"]:
        ap.add_argument(name,type=Path)
    args = ap.parse_args()
    screen(args.survey_cache,args.inputs,args.out)
