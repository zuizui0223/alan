"""Summarize source-verified regional flower-visit observations for two focal plants."""
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

from alan_gecko_primary_fiveyear import assign_points, prepare_land

PLANTS = ["Asclepias syriaca", "Cirsium vulgare"]
YEARS = {2019, 2021, 2023, 2024, 2025}
FIELDS = ["source_taxon_name", "target_taxon_name", "latitude", "longitude", "event_date",
          "study_url", "source_specimen_occurrence_id", "source_taxon_path"]


def download(work):
    folder = work / "local21_vis"
    folder.mkdir(exist_ok=True)
    pages = []
    for plant in PLANTS:
        path = folder / f"{plant}.parquet"
        if not path.exists():
            blocks = []
            while True:
                params = dict(interactionType="visitsFlowersOf", bbox="-75,40,-70,45", includeObservations="true",
                              limit=1024, field=FIELDS, targetTaxon=plant, offset=1024 * len(blocks))
                response = requests.get("https://api.globalbioticinteractions.org/interaction.csv", params=params, timeout=90)
                response.raise_for_status()
                block = pd.read_csv(BytesIO(response.content))
                blocks.append(block)
                if len(block) < 1024:
                    break
            pd.concat(blocks).to_parquet(path)
        pages.append(pd.read_parquet(path))
    records = pd.concat(pages)
    dated = records.loc[records.event_date.str[:4].isin(map(str, YEARS)) & records.target_taxon_name.isin(PLANTS)]
    ids = sorted(dated.study_url.str.extract(r"inaturalist.org/observations/(\d+)")[0].dropna().astype(int).unique())
    payloads = []
    for k in range(0, len(ids), 100):
        expected = set(ids[k:k+100])
        path = folder / f"inat_{k:05}.json"
        if not path.exists():
            response = requests.get("https://api.inaturalist.org/v1/observations",
                                    params=dict(id=",".join(map(str, sorted(expected))), per_page=100), timeout=90)
            response.raise_for_status()
            fields = ["id","observed_on","positional_accuracy","obscured","geoprivacy","taxon_geoprivacy",
                      "location","quality_grade","ofvs","taxon"]
            extracts = [{c: row.get(c) for c in fields} for row in response.json()["results"]]
            path.write_text(json.dumps(extracts), encoding="utf8")
            time.sleep(1)
        part = json.loads(path.read_text(encoding="utf8"))
        assert all(row["id"] in expected for row in part), "stale iNaturalist cache"
        payloads.extend(part)
    return payloads, len(ids)


def summarize(work, out):
    payloads, nids = download(work)
    rows = []
    for row in payloads:
        year = int(row["observed_on"][:4])
        if year not in YEARS or row["obscured"] or row["geoprivacy"] in {"private","obscured"}:
            continue
        if row["quality_grade"] != "research" or not 0 < (row["positional_accuracy"] or 0) <= 1000:
            continue
        y,x = map(float, row["location"].split(","))
        if not (-75 <= x < -70 and 40 <= y < 45):
            continue
        taxon = row["taxon"]
        if taxon["rank_level"] > 10 or 47158 not in taxon["ancestor_ids"]:
            continue
        ancestors = set(taxon["ancestor_ids"])
        if 47157 in ancestors:
            group = "butterfly" if 47224 in ancestors else "moth"
        else:
            group = next((name for key,name in [(47201,"Hymenoptera"),(47822,"Diptera"),(47208,"Coleoptera")]
                          if key in ancestors),"other_insect")
        for field in row["ofvs"]:
            plant = (field.get("taxon") or {}).get("name")
            if field["field_id"] == 3126 and plant in PLANTS:
                rows.append(dict(unit_id=row["id"], plant=plant, year=year, date=row["observed_on"],
                                 taxon=" ".join(taxon["name"].split()[:2]), group=group, longitude=x, latitude=y))
    visits = pd.DataFrame(rows).drop_duplicates(["unit_id","plant"])
    coast = work / "v011/coast"
    points = visits[["unit_id","longitude","latitude"]].drop_duplicates()
    classified = assign_points(points, prepare_land(gpd.read_file(coast / "GSHHS_h_L1.shp")),
                               gpd.read_file(coast / "GSHHS_h_L2.shp"))
    v = visits.merge(classified[["unit_id","land_id","land_class"]], on="unit_id", validate="m:1")
    v = v[v.land_class.isin(["island","mainland"])].copy()
    v["key"] = (v.year.astype(str) + "|" + np.floor((v.longitude+180)*20).astype(int).astype(str) + "_" +
                np.floor((v.latitude+90)*20).astype(int).astype(str))
    photo = pd.read_parquet(work / "phenology_run/target.parquet").query("region5 == '21_26' and plant in @PLANTS")
    land = pd.read_parquet(work / "phenology_run/strict_site_light.parquet")[["key","land_id"]].drop_duplicates()
    photo = photo.merge(land, on="key", validate="m:1")
    for frame in (v, photo):
        frame["land_id"] = frame.land_id.astype(str)
    keys = ["plant","land_class","land_id","key"]
    exact = v.merge(photo[keys].drop_duplicates(), on=keys, validate="m:1")
    same_day = exact.merge(photo[keys+["observed_on"]].drop_duplicates().rename(columns={"observed_on":"date"}),
                           on=keys+["date"],validate="m:1")
    output = []
    for level, frame in [("region",v),("same_focal_cell_year",exact),("same_focal_cell_date",same_day)]:
        outcounts=frame.groupby(["plant","land_class","year","group"]).agg(
            records=("unit_id","nunique"), visitor_taxa=("taxon","nunique"), occupied_cell_years=("key","nunique"))
        output.append(outcounts.reset_index().assign(scale=level))
        print(level,frame.groupby(["plant","land_class"]).size().to_dict())
        print(frame.groupby(["plant","land_class","group"]).size().to_dict())
        print("pairs",frame.groupby(["plant","year"]).land_class.nunique().loc[lambda x:x==2].index.tolist())
    out.mkdir(exist_ok=True,parents=True)
    pd.concat(output,ignore_index=True).to_csv(out / "verified_visitors.csv",index=False)
    taxon_counts = v.groupby(["plant","group","taxon","land_class"]).unit_id.nunique().unstack(fill_value=0)
    taxon_counts.to_csv(out / "visitor_taxon_counts.csv")
    for plant in PLANTS:
        both = taxon_counts.loc[plant].query("island > 0 and mainland > 0").sort_values("island",ascending=False)
        print(plant,"shared taxa",len(both),both.head(12).to_string())
        print(both.loc["moth"].head(5) if "moth" in both.index.get_level_values(0) else "")
    print("links", v[v.taxon.isin(["Bombus impatiens","Hemaris thysbe","Hemaris diffinis"])].groupby(["plant","taxon","land_class"]).unit_id.first().to_dict())
    digest=sha256(("\n".join(map(str,sorted(v.unit_id.unique())))+"\n").encode()).hexdigest()
    meta=dict(candidate_inat_ids=nids,api_records_returned=len(payloads),included_visits=len(v),
              occurrence_id_set_sha256=digest, exact_same_date_records=len(same_day))
    (out / "visitor_receipt.json").write_text(json.dumps(meta,indent=2)+"\n",encoding="utf8")
    print(meta)
    print("exact examples",exact.groupby(["plant","land_class","group"]).first()[["taxon","unit_id"]])
    return v,exact


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path)
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    summarize(args.work,args.out)
