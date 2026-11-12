"""Summarize archived assays and observation support for the six largest signed contributions."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

from alan_gecko_primary_fiveyear import assign_points, prepare_land, LINKED_MOTHS, SPECIES
from alan_repeat_flower_contrasts import read_source
from alan_shared_visitors_regional import grid, uk

RESULTS = Path(__file__).resolve().parents[2] / "results"
ARMS = ["land_class", "light"]


def read_occurrences(paths, named_keys, cells, *, adult, counts=None):
    records = {}
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            page = json.load(stream)
        expected = counts[path.name] if counts is not None else page["count"]
        assert expected == len(page["rows"])
        for row in page["rows"]:
            if named_keys.get(row["speciesKey"]) != row["species"]:
                continue
            if row["basisOfRecord"] != "HUMAN_OBSERVATION" or row["occurrenceStatus"] not in ("PRESENT", "DETECTED"):
                continue
            if not 0 < (row["coordinateUncertaintyInMeters"] or 0) <= 1000:
                continue
            if adult and str(row["lifeStage"]).lower() not in ("adult", "imago"):
                continue
            x, y = row["decimalLongitude"], row["decimalLatitude"]
            if (row["year"], int(np.floor(x)), int(np.floor(y))) in cells:
                records[row["key"]] = (row["key"], row["year"], row["species"], x, y)
    return pd.DataFrame(records.values(), columns=["unit_id", "year", "species", "longitude", "latitude"])


def main(work, out):
    ledger = pd.read_csv(RESULTS / "species_interpretation/shared_544.csv").sort_values("gap_pp")
    focal = pd.concat([ledger.head(3), ledger.tail(3).iloc[::-1]]).plant
    pred = work / "predator"
    plants = pd.read_parquet(pred / "arthropod/plants.parquet").rename(columns={"wcvp_accepted_name": "plant"})
    sites = pd.read_parquet(pred / "dual80/sites.parquet")
    assert len(sites) == sites.key.nunique() == 10698
    obs = sites.merge(plants[plants.plant.isin(focal)], on="key", validate="1:m")
    cells = set(obs[["year", "x1", "y1"]].itertuples(index=False, name=None))
    v = pd.read_parquet(pred / "arthropod/visits.parquet").rename(columns={"wcvp_accepted_name": "plant"})
    edges = pd.read_csv(pred / "prey/species_prey.csv").rename(columns={"gecko": "predator"})
    arth = pd.read_csv(pred / "arthropod/edges/species_prey.csv")
    pairs = pd.concat([edges, arth]).query("source_stage == 'adult'").merge(v[v.plant.isin(focal)], on="moth")
    pairs = pairs[["predator", "moth", "plant"]].drop_duplicates()
    query_dir = pred / "arthropod"
    manifest = json.loads((query_dir / "moth_queries.json").read_text())
    keys = {q["taxon"]: q["moth"] for q in manifest if q["moth"] in set(pairs.moth)}
    moths = pd.concat([
        read_occurrences((query_dir / "moth_queries" / q["file"] for q in manifest if q["taxon"] in keys), keys, cells, adult=True,
                         counts={q["file"]: q["count"] for q in manifest}),
        read_occurrences((pred / "moths_private").glob("*.json.gz"), LINKED_MOTHS, cells, adult=True)])
    geckos = read_occurrences((pred / "private").glob("*.json.gz"), SPECIES, cells, adult=False)
    coast = work / "v011/coast"
    polygons = prepare_land(gpd.read_file(coast / "GSHHS_h_L1.shp"))
    lakes = gpd.read_file(coast / "GSHHS_h_L2.shp")
    animal = pd.concat([geckos, moths]).drop_duplicates("unit_id")
    animal = animal.merge(assign_points(animal, polygons, lakes), on="unit_id", validate="1:1")
    animal["x1"] = np.floor(animal.longitude).astype(int)
    animal["y1"] = np.floor(animal.latitude).astype(int)
    animal["gx"] = np.floor((animal.longitude + 180) * 20).astype(int)
    animal["gy"] = np.floor((animal.latitude + 90) * 20).astype(int)
    arth_sites = pd.read_parquet(query_dir / "predator_sites.parquet")
    hits = []
    for scale, axes in [("1deg", ["x1", "y1"]), ("005deg", ["gx", "gy"])]:
        cols = ["year", "land_id", "land_class", *axes]
        geo = animal[cols + ["species"]].drop_duplicates()
        row = obs.merge(pairs, on="plant").merge(geo.rename(columns={"species": "moth"}), on=cols + ["moth"])
        gp = row.merge(geo.rename(columns={"species": "predator"}), on=cols + ["predator"])
        ap = row.merge(arth_sites.loc[arth_sites.scale.eq("1deg" if scale == "1deg" else "0.05deg"),
                                     ["key", "predator"]].drop_duplicates(), on=["key", "predator"])
        hits.append(pd.concat([gp, ap]).assign(scale=scale))
    matched = pd.concat(hits)
    folder = work / "phenology_run"
    phen = pd.concat([read_source(n, folder / visits, folder / exp) for n, visits, exp in [
        ("NPN", "npn/individual_visits.parquet", "npn_exposure.parquet"),
        ("Sweden", "sweden_all_visits.parquet", "sweden_all_exposure.parquet"),
        ("ClimateWatch", "australia_standard_visits.parquet", "australia_exposure.parquet")]])
    assert len(phen) == 82369
    phen = phen[phen.plant.isin(focal)]
    survey = pd.read_csv(RESULTS / "bilateral_visitors/us_candidate_events.csv")
    survey = survey[survey.local_photo & survey.minutes.gt(0) & survey.broad_nonlep.notna()].copy()
    pom = uk(work / "visitor_effort").dropna(subset=["x", "y", "original"])
    pom["plant"] = pom.original.str.strip().str.casefold().map({x.casefold(): x for x in focal})
    pom["region5"] = grid(pom.x, pom.y)
    local = pd.read_parquet(work / "visitor_effort/shared.parquet").rename(columns={"wcvp_accepted_name": "plant"})
    pom = pom.merge(local[["plant", "year", "region5"]].drop_duplicates(), on=["plant", "year", "region5"])
    pom = pom.rename(columns={"x": "longitude", "y": "latitude"}).reset_index(drop=True)
    pom["unit_id"] = pom.index
    pom = pom.merge(assign_points(pom, polygons, lakes), on="unit_id", validate="1:1")
    pom["minutes"] = 10
    survey = pd.concat([survey, pom], ignore_index=True)
    survey = survey[survey.land_class.isin(["island", "mainland"])]
    print(survey[survey.plant.isin(focal)].groupby(["plant", "year", "region5", "land_class", "source"]).minutes.agg(["sum", "size"]).to_string())
    reviewed = pd.read_csv(RESULTS / "dual_role_80/island_cohort_reviewed_links.csv")
    reviewed = reviewed[reviewed.assessment.eq("predation_recorded")]
    rows = []
    for plant in focal:
        f = phen[phen.plant.eq(plant)]
        e = survey[survey.plant.eq(plant)]
        pair = pairs[pairs.plant.eq(plant)]
        h = matched[matched.plant.eq(plant)]
        r = dict(plant=plant, named_adult_prey_links=" | ".join(pair.predator + " > " + pair.moth),
                 source_checked_same_plant_links_in_island_dualrole_cohort=int(reviewed.plant.eq(plant).sum()))
        for scale in ["1deg", "005deg"]:
            for side in ["island", "mainland"]:
                r[f"triad_{scale}_{side}_siteyears"] = h[h.scale.eq(scale) & h.land_class.eq(side)].key.nunique()
        for side in ["island", "mainland"]:
            z = f[f.land_class.eq(side)]
            r[f"flower_status_records_{side}"] = len(z)
            r[f"flower_open_records_{side}"] = int(z.status.sum())
            r[f"timed_visitor_minutes_{side}"] = int(e.loc[e.land_class.eq(side), "minutes"].sum())
        r["phenology_region_years"] = len(f[["region5", "year"]].drop_duplicates())
        r["phenology_max_arms_in_one_region_year"] = int(f.drop_duplicates(["year", "region5", *ARMS]).groupby(["year", "region5"]).size().max()) if len(f) else 0
        r["bilateral_timed_same_plant_region_years"] = int(e.groupby(["year", "region5"]).land_class.nunique().eq(2).sum())
        rows.append(r)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.work, args.output)
