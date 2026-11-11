"""Download gecko human observations for the five original years. Keep coordinates private."""
import argparse
import concurrent.futures
import gzip
import json
import time
from pathlib import Path

import requests

YEARS = (2019, 2021, 2023, 2024, 2025)
SPECIES = {9537238: "Hemidactylus frenatus", 5959942: "Hemidactylus mabouia"}
MULTIGUILD = {2433194:"Glossophaga soricina", 2432339:"Antrozous pallidus", 2433346:"Mystacina tuberculata",
              2231197:"Thereuopoda clunifera", 5179562:"Scolopendra subspinipes", 1311477:"Vespa velutina",
              5170994:"Argiope bruennichi", 6258028:"Mantis religiosa"}
LINKED_MOTHS = {1886251:"Spoladea recurvalis", 7400814:"Feltia subterranea"}
FIELDS = ("key", "speciesKey", "species", "year", "basisOfRecord", "occurrenceStatus",
          "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "datasetKey", "lifeStage")


URL = "https://api.gbif.org/v1/occurrence/search"


def page(session, params, offset=0):
    for retry in range(6):
        try:
            response = session.get(URL, params={**params, "offset":offset}, timeout=45)
            response.raise_for_status()
            return response.json()
        except requests.RequestException:
            if retry == 5:
                raise
            time.sleep(2**retry)


def fetch_query(session, params):
    data = page(session, params)
    count = data["count"]
    if count > 100000:
        raise ValueError("GBIF paging ceiling exceeded for a single geographic partition")
    rows = []
    while True:
        assert data["count"] == count and data["offset"] == len(rows)
        rows.extend({k:r.get(k) for k in FIELDS} for r in data["results"])
        if data["endOfRecords"]:
            break
        if not data["results"]:
            raise ValueError("Incomplete GBIF page")
        data = page(session, params, len(rows))
    assert len(rows) == count == len({r["key"] for r in rows})
    return rows


def collect(job):
    root, year, taxon, scope = job
    dest = root / f"{year}_{taxon}.json.gz"
    if dest.exists():
        return
    params = dict(taxonKey=taxon, year=year, hasCoordinate="true",
                  basisOfRecord="HUMAN_OBSERVATION", limit=300)
    queries = [params]
    with requests.Session() as session:
        total = page(session, {**params, "limit":0})["count"]
        if scope and total > 12000:
            queries = [{**params, "decimalLongitude":f"{x},{x+5}", "decimalLatitude":f"{y},{y+5}"}
                       for x,y in scope["blocks"][str(year)]]
            partitions = []
            for query in queries:
                size = page(session, {**query, "limit":0})["count"]
                if size <= 100000:
                    partitions.append(query)
                    continue
                xmin = int(query["decimalLongitude"].split(",")[0])
                ymin = int(query["decimalLatitude"].split(",")[0])
                for x,y in scope["one_degree"][str(year)]:
                    if xmin <= x < xmin+5 and ymin <= y < ymin+5:
                        partitions.append({**params,"decimalLongitude":f"{x},{x+1}","decimalLatitude":f"{y},{y+1}"})
            queries = partitions
        records = {}
        counts = []
        for query in queries:
            rows = fetch_query(session,query)
            counts.append(len(rows))
            for row in rows:
                # Bounds are inclusive: neighbouring boxes may share exact-boundary events.
                if row["key"] in records:
                    assert records[row["key"]] == row
                records[row["key"]] = row
    with gzip.open(dest.with_suffix(".tmp"), "wt", encoding="utf-8") as f:
        json.dump({"url":URL, "query":params, "queries":queries, "query_counts":counts,
                   "world_query_count":total, "count":len(records), "rows":list(records.values())}, f)
    dest.with_suffix(".tmp").replace(dest)
    print(year, taxon, len(records), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("private_directory", type=Path)
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--multiguild", action="store_true", help="Use the eight taxa from the upstream predator-role pilot")
    target.add_argument("--linked-moths", action="store_true", help="Use the two moth species with a documented prey and flower link")
    parser.add_argument("--blocks", type=Path, help="Complete set of 5-degree boxes enclosing frozen sites by year")
    args = parser.parse_args()
    scope = None
    if args.blocks:
        scope = json.loads(args.blocks.read_text())
        assert scope["target_hash"] == "37f3c2c519871b6905ffa81c5a14536b248611861ad143b1d4d77d95b7c6b00d"
    root = args.private_directory
    root.mkdir(parents=True, exist_ok=True)
    taxa = LINKED_MOTHS if args.linked_moths else MULTIGUILD if args.multiguild else SPECIES
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(collect, [(root, year, taxon, scope) for year in YEARS for taxon in taxa]))
