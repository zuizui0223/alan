"""Join dated focal observations to the original regional contributions and annual land cover."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from pyproj import Geod

from alan_shared_species_contributions import primary, influence

NAMES = ("Asclepias syriaca", "Cirsium vulgare")
REGION = "21_26"


def annotations(folder, focal, cache):
    observations = {}
    for path in (folder / "responses_by_id_type").glob("*.json"):
        page = json.loads(path.read_bytes())
        observations.update((str(obs[page["param"]]), obs) for obs in page["results"])
    uuids = sorted({observations[i]["uuid"] for i in focal.observation_id})
    if not cache.exists():
        response = requests.get("https://api.inaturalist.org/v2/observations",
                                params={"uuid": ",".join(uuids), "per_page": 100,
                                        "fields": "uuid,user.id,positional_accuracy,geoprivacy,obscured,observed_on,time_observed_at,photos.id"},
                                timeout=60)
        response.raise_for_status()
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(response.content)
    context = {o["uuid"]: o for o in json.loads(cache.read_bytes())["results"]}
    assert set(context) == set(uuids)
    result = []
    for row in focal.itertuples():
        obs = observations.get(row.observation_id, {})
        ok = obs.get("observed_on") == row.observed_on and row.photo_id in {p["id"] for p in obs.get("photos", [])}
        values = {a["controlled_value_id"] for a in obs.get("annotations", [])
                  if a["controlled_attribute_id"] == 12 and (a.get("vote_score") or 0) >= 0} if ok else set()
        meta = context[obs["uuid"]]
        assert meta["observed_on"] == row.observed_on and row.photo_id in {p["id"] for p in meta["photos"]}
        result.append((ok, 13 in values, 21 in values, meta["positional_accuracy"], meta["user"]["id"], meta["obscured"]))
    return pd.DataFrame(result, columns=["exact_date_photo", "flowering", "no_evidence_of_flowering", "accuracy_m", "observer", "obscured"], index=focal.index)


def run(work, out, cover_path):
    folder = work / "phenology_run"
    all_photos = pd.read_parquet(folder / "universe.parquet").query("region5 == @REGION")
    photos = pd.read_parquet(folder / "target.parquet").query("region5 == @REGION and plant in @NAMES").copy()
    assert photos.photo_id.is_unique and photos.observation_id.is_unique
    assert len(photos) == 73 and photos.key.nunique() == 37
    photos = photos.join(annotations(folder, photos, work / "local_21_26_private/photo_context.json"))
    sites, _, _, _ = primary(work / "v011", work / "mois/signed")
    sites["contribution_pp"] = 100 * influence(sites)[1] / sites.n_species
    sites = sites[sites.region5.eq(REGION)]
    assert len(sites) == all_photos.key.nunique() == 835
    exposure = sites[["key", "n_species", "contribution_pp", "photo_attempts", "pop_density"]]
    background = all_photos.groupby("key").observed_on.nunique().rename("plant_photography_days")
    counts = photos.groupby(["plant", "year", "key", "land_class", "viirs_median"]).agg(
        observations=("photo_id", "size"), observers=("observer", "nunique"),
        accuracy_reported=("accuracy_m", "count"), max_reported_accuracy_m=("accuracy_m", "max"),
        obscured_observations=("obscured", "sum"), observation_dates=("observed_on", lambda x: "|".join(sorted(set(x)))),
        exact_photo_dates=("exact_date_photo", "sum"), flowering_annotations=("flowering", "sum"),
        no_evidence_of_flowering=("no_evidence_of_flowering", "sum")).reset_index()
    counts = counts.merge(exposure, on="key", validate="m:1").join(background, on="key")
    np.testing.assert_allclose(counts.groupby("plant").contribution_pp.sum().loc[list(NAMES)],
                               [-0.19643862937088152, 0.13959816725703159], atol=1e-10)
    obsgeo = pd.read_csv(work / "mois/signed/canonical/canonical_observations_climate.csv.gz",
                         usecols=["photo_id", "latitude", "longitude"])
    shared = photos[photos.key.isin(counts[counts.key.duplicated(False)].key)].merge(
        obsgeo[obsgeo.photo_id.isin(photos.photo_id)].drop_duplicates(), on="photo_id", validate="1:1")
    assert len(shared) == 2
    distance = Geod(ellps="WGS84").inv(*shared.iloc[0][["longitude", "latitude"]],
                                       *shared.iloc[1][["longitude", "latitude"]])[2]
    summary = dict(all_regional_site_years=835, occupied_site_years=37, focal_species_site_years=38,
                   physical_grids=photos.key.str.split("|", regex=False).str[-1].nunique(),
                   shared_cell_year=shared.year.iloc[0].item(), shared_land=shared.land_class.iloc[0],
                   shared_dates=shared.set_index("plant").observed_on.to_dict(),
                   shared_distance_m=round(distance), date_match=int(photos.exact_date_photo.sum()),
                   observation_annotation_counts=photos.groupby("plant")[["flowering","no_evidence_of_flowering"]]
                       .sum().to_dict("index"))
    strata = ["year", "land_class", "month"]
    month = photos.groupby(strata + ["plant"]).agg(site_years=("key", "nunique"),
                                                   photos=("photo_id", "size"),
                                                   flowering=("flowering", "sum")).reset_index()
    month = month.merge(all_photos.groupby(strata).observed_on.nunique().rename("background_days").reset_index(),
                        on=strata, validate="m:1")
    out.mkdir(parents=True, exist_ok=True)
    precision = photos.groupby(["plant", "land_class"]).agg(
        observations=("photo_id", "size"), distinct_observers=("observer", "nunique"),
        obscured=("obscured", "sum"), reported_accuracy=("accuracy_m", "count"),
        median_reported_accuracy_m=("accuracy_m", "median"),
        within_100m=("accuracy_m", lambda x: x.le(100).sum()))
    precision.to_csv(out / "photo_precision.csv")
    month.to_csv(out / "photographic_months.csv", index=False)
    if cover_path is not None:
        cover = pd.read_csv(cover_path)
        cols = [c for c in cover if c.startswith("cell_") or c.startswith("photo_")]
        counts = counts.merge(cover[["key", "plant", *cols]], on=["key","plant"], validate="1:1", how="left",
                              indicator=True)
        assert counts.pop("_merge").eq("both").all()
    counts = counts.sort_values(["plant", "year", "key"])
    cell = counts.key.str.split("|", regex=False).str[-1]
    counts["site"] = cell
    assert counts.site.nunique() == cell.nunique()
    counts.to_csv(out / "focal_dated_sites.csv", index=False)
    (out / "photo_support.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    print(month[month.duplicated(strata, keep=False)].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--cover", type=Path)
    args = parser.parse_args()
    run(args.work, args.out, args.cover)
