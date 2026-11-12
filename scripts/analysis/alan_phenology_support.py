"""Summarize date, flowering and effort support without imputing unobserved visits."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "results/phenology_support"
ARMS = ["land_class", "light"]


def support(data, groups):
    days = data.groupby(groups + ARMS).observed_on.nunique().unstack(ARMS)
    assert days.shape[1] == 4
    return days.dropna().min(axis=1).rename("min_days").astype(int).reset_index()


def run(folder):
    x = pd.read_parquet(folder / "target.parquet")
    u = pd.read_parquet(folder / "universe.parquet")
    objects = {}
    requested = []
    for file in (folder / "responses_by_id_type").glob("*.json"):
        batch = json.loads(file.read_text())
        requested.extend(batch["requested"])
        objects.update({str(r[batch["param"]]): r for r in batch["results"]})
    assert sorted(requested) == sorted(x.observation_id.unique())
    matched, open_at_obs, explicit_none, uuid = [], [], [], []
    for row in x.itertuples():
        obs = objects.get(row.observation_id, {})
        ok = obs.get("observed_on") == row.observed_on and row.photo_id in {p["id"] for p in obs.get("photos", [])}
        values = {a["controlled_value_id"] for a in obs.get("annotations", [])
                  if a["controlled_attribute_id"] == 12 and (a.get("vote_score") or 0) >= 0} if ok else set()
        matched.append(ok)
        open_at_obs.append(13 in values)
        explicit_none.append(21 in values)
        uuid.append(obs.get("uuid"))
    x = x.assign(matched=matched, flowering=open_at_obs, explicit_none=explicit_none)
    x["photo_hash"] = x.photo_id.map(lambda i: hashlib.sha256(str(i).encode()).hexdigest())
    keys = ["plant", "region5", "year"]
    annual = support(x, keys)
    candidates = x.merge(annual, on=keys)
    review = pd.read_csv(OUT / "photo_review.csv")
    assert set(review.photo_hash) == set(candidates.photo_hash)
    candidates = candidates.merge(review, on=["photo_hash","plant"], validate="1:1")
    monthly = support(x, keys + ["month"])
    background = support(u, ["region5", "year", "month"])
    effort = background[background.min_days >= 2].rename(columns={"min_days":"background_min_days"})
    monthly = monthly.merge(effort, on=["region5","year","month"], how="left")
    monthly.to_csv(OUT / "monthly_overlap.csv", index=False)
    annual.to_csv(OUT / "annual_overlap.csv", index=False)
    species = x.groupby("plant").agg(photos=("photo_hash","size"),regions=("region5","nunique"),
            date_photo_matched=("matched","sum"),flowering_annotations=("flowering","sum"),no_flowering_signs=("explicit_none","sum"))
    species.to_csv(OUT / "species_544.csv")
    visible = support(candidates[candidates.visible.eq("open")], keys)
    bg_months = effort.groupby(["region5","year"]).size()
    checked = monthly.dropna()
    results = dict(photos=len(x), regions=x.region5.nunique(), species=len(species),
        universe_photos=len(u), universe_region_years=u.region_year.nunique(), universe_regions=u.region5.nunique(),
        api_observations_found=len(objects), distinct_public_uuids=pd.Series(uuid).nunique(),
        exact_photo_and_date=sum(matched), flowering_annotations=sum(open_at_obs),
        plants_with_flowering_annotation=int((species.flowering_annotations>0).sum()), no_signs=sum(explicit_none),
        conflicting=sum(a and b for a,b in zip(open_at_obs,explicit_none)),
        four_arm_species=annual.plant.nunique(),four_arm_regions=annual.region5.nunique(),
        four_arm_photos=len(candidates),four_arm_year_min_days=annual.min_days.value_counts().to_dict(),
        common_month_species=monthly.plant.nunique(),common_month_regions=monthly.region5.nunique(),
        with_background_effort_species=checked.plant.nunique(),with_background_effort_regions=checked.region5.nunique(),
        max_common_months=int(monthly.groupby(keys).size().max()),
        background_two_day_region_years=len(bg_months),background_two_day_regions=effort.region5.nunique(),
        three_month_background_region_years=int((bg_months>=3).sum()),visual=candidates.visible.value_counts().to_dict(),
        all_four_visibly_open=visible.to_dict("records"),exif=candidates.original_date_check.value_counts().to_dict())
    (OUT / "summary.json").write_text(json.dumps(results, indent=2, default=int)+"\n")
    print(json.dumps(results,default=int))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("private_cache", type=Path)
    run(parser.parse_args().private_cache)
