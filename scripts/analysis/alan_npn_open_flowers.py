"""Acquire USA-NPN individual open-flower status for the frozen years and WCVP cohort."""
import argparse
import calendar
import concurrent.futures
import json
import time
from pathlib import Path

import pandas as pd
import requests

from alan_globi_v011_same_fiveyear import YEARS, taxonomy

BASE = "https://services.usanpn.org/npn_portal/"
ROOT = Path(__file__).resolve().parents[2]
COLUMNS = ["observation_id","site_id","latitude","longitude","species_id","individual_id",
           "phenophase_id","observation_date","phenophase_status","intensity_value","intensity_category_id"]


def get_catalog(folder, endpoint, filename):
    path = folder / filename
    if not path.exists():
        result = requests.get(BASE + endpoint + ".json", timeout=90)
        result.raise_for_status()
        path.write_bytes(result.content)
    return pd.read_json(path)


def run(folder, wcvp):
    folder.mkdir(parents=True, exist_ok=True)
    catalog = get_catalog(folder, "species/getSpecies", "species.json")
    get_catalog(folder, "stations/getAllStations", "getAllStations.json")
    phases = get_catalog(folder, "phenophases/getPhenophases", "getPhenophases.json")
    phase_ids = sorted(phases.loc[phases.phenophase_name.str.strip().str.lower().eq("open flowers"), "phenophase_id"])
    plants = catalog[catalog.kingdom.eq("Plantae")].copy()
    plants["plant_name"] = plants.genus.str.strip() + " " + plants.species.str.strip()
    tax = taxonomy(wcvp, set(plants.plant_name), set())
    target = set(pd.read_csv(ROOT / "results/moth_pollen_evidence/shared_544_evidence.csv").plant)
    allowed = plants.merge(tax[["accepted_name","resolution"]], on="plant_name")
    allowed = allowed[allowed.accepted_name.isin(target) & allowed.resolution.eq("resolved")]
    crosswalk = allowed[["species_id","plant_name","accepted_name"]].drop_duplicates().sort_values("species_id")
    crosswalk.to_csv(folder / "target_names.csv", index=False)
    species_ids = list(crosswalk.species_id)
    print("species",len(species_ids),"target",crosswalk.accepted_name.nunique(),"phase ids",phase_ids,flush=True)

    def fetch(period):
        year, month = period
        path = folder / f"open_{year}_{month:02}.parquet"
        if path.exists():
            return
        start = f"{year}-{month:02}-01"
        end = f"{year}-{month:02}-{calendar.monthrange(year,month)[1]}"
        params = {"request_src":"ALAN island mainland flowering phenology research", "start_date": start, "end_date": end}
        for key, ids in [("phenophase_id",phase_ids),("species_id",species_ids)]:
            params.update({f"{key}[{i}]": int(value) for i,value in enumerate(ids, 1)})
        for attempt in range(4):
            try:
                result = requests.post(BASE + "observations/getObservations.json", data=params, timeout=(15,150))
                result.raise_for_status()
                data = pd.DataFrame(result.json(), columns=COLUMNS)
                assert set(data.species_id) <= set(species_ids) and set(data.phenophase_id) <= set(phase_ids)
                assert data.observation_date.between(start, end).all() and data.observation_id.is_unique
                data["intensity_value"] = data.intensity_value.astype("string")
                tmp = path.with_suffix(".tmp")
                data.to_parquet(tmp,index=False)
                tmp.replace(path)
                print(start,len(data),flush=True)
                return
            except (ValueError, requests.RequestException):
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        list(pool.map(fetch, [(year, month) for year in sorted(YEARS) for month in range(1,13)]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("private_cache", type=Path)
    parser.add_argument("wcvp_directory", type=Path)
    args = parser.parse_args()
    run(args.private_cache, args.wcvp_directory)
