"""Cache public observation-level dates and phenology for the frozen photo IDs."""
import argparse
import json
import time
from pathlib import Path

import pandas as pd
import requests

FIELDS = ("id,uuid,observed_on,photos.id,photos.url,annotations.controlled_attribute_id,"
          "annotations.controlled_value_id,annotations.vote_score")


def run(folder):
    photos = pd.read_parquet(folder/"target.parquet")
    ids_all = set(photos.observation_id)
    cache = folder/"responses_by_id_type"
    cache.mkdir(exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "alan-research-phenology/1.0"
    batches = [(param,group[start:start+100])
               for param,group in [("id",sorted(x for x in ids_all if x.isdecimal())),
                                   ("uuid",sorted(x for x in ids_all if not x.isdecimal()))]
               for start in range(0,len(group),100)]
    for param,ids in batches:
        path = cache/(ids[0]+".json")
        if path.exists() and json.loads(path.read_text())["requested"] == ids:
            continue
        for attempt in range(5):
            try:
                r = session.get("https://api.inaturalist.org/v2/observations",
                                params={param:",".join(ids),"per_page":100,"fields":FIELDS},timeout=35)
                r.raise_for_status()
                response = r.json()
                objects = response["results"]
                assert len(objects)==response["total_results"]
                assert set(str(x[param]) for x in objects)<=set(ids)
                break
            except (requests.RequestException,ValueError):
                if attempt==4:
                    raise
                time.sleep(5*2**attempt)
        tmp=path.with_suffix(".tmp")
        tmp.write_text(json.dumps({"param":param,"requested":ids,"results":objects}))
        tmp.replace(path)
        time.sleep(1.1)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("private_cache",type=Path)
    run(parser.parse_args().private_cache)
