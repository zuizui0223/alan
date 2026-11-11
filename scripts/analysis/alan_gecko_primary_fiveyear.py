"""Condition the original plant composition model on documented same-year house gecko overlap."""
import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box
from shapely.ops import unary_union

from alan_globi_v011_same_fiveyear import ROOT, ID, TERM, environment, site_species
from alan_globi_crossland_light_robust_20261111 import prepare, fit

sys.path.insert(0, str(ROOT))
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import assign_points, prepare_land
from scripts.acquisition.alan_gecko_gbif_fiveyear import SPECIES, MULTIGUILD, LINKED_MOTHS, YEARS


def primary(inputs, signed):
    reference = json.loads((ROOT / "results/globi_v011_fiveyear/comparison.json").read_text(encoding="utf-8"))
    sites, species = site_species(signed, environment(inputs))
    legacy = pd.read_csv(inputs / "old/GloBI_stage2_signed_extension.csv.gz")
    names = legacy.loc[legacy.guild.eq("moth_family_visit"), "plant_name"]
    species["all_moth"] = species.wcvp_accepted_name.isin(names)
    counts = species.groupby("key").agg(n_species=(ID, "size"), all_moth=("all_moth", "sum"))
    s = sites.merge(counts, on="key", validate="one_to_one")
    s["island"] = s.land_class.eq("island").astype(float)
    s["region_year"] = s.year.astype(str) + "|" + s.region5
    for out, col in [("loglight", "viirs_median"), ("logpop", "pop_density"), ("logphotos", "photo_attempts"),
                     ("logspecies", "n_species"), ("lograin", "rain")]:
        s[out] = np.log1p(s[col].clip(lower=0))
    s["elevation_km"] = s.elevation / 1000
    s["light_island"] = s.loglight * s.island
    s = prepare(s, 3, 2).dropna(subset=["bio5", "lograin", "elevation_km", "logpop", "loglight"])
    s = s[s.groupby("region_year").land_class.transform("nunique").eq(2)].copy()
    digest = hashlib.sha256(("\n".join(sorted(s.target_id)) + "\n").encode()).hexdigest()
    assert digest == reference["target_key_set_sha256"]
    result = fit(s, "all", "moth", balance=True)
    np.testing.assert_allclose(result[TERM], reference["fits"]["published"][TERM], rtol=0, atol=1e-10)
    return s, digest, species[species.key.isin(s.key)], legacy


def positive_records(folder, coast, taxa=SPECIES, expected_blocks=None):
    data, hashes, queries = [], {}, []
    for year in YEARS:
        for taxon in taxa:
            path = folder / f"{year}_{taxon}.json.gz"
            with gzip.open(path, "rt", encoding="utf-8") as f:
                obj = json.load(f)
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            assert obj["count"] == len(obj["rows"]) and obj["query"]["taxonKey"] == taxon and obj["query"]["year"] == year
            partitions = obj.get("queries", [obj["query"]])
            spatial = any("decimalLatitude" in q for q in partitions)
            if spatial:
                rectangles = []
                for q in partitions:
                    x, xmax = map(float, q["decimalLongitude"].split(","))
                    y, ymax = map(float, q["decimalLatitude"].split(","))
                    assert 0 < xmax-x <= 5 and abs((xmax-x)-(ymax-y)) < 1e-8
                    rectangles.append(box(x,y,xmax,ymax))
                coverage = unary_union(rectangles)
                assert all(coverage.covers(box(x,y,x+1,y+1)) for x,y in expected_blocks[str(year)])
            queries.append(dict(year=year,taxon=taxa[taxon],rows=obj["count"],
                                coverage="boxes covering all frozen-site 1-degree grids" if spatial else "world"))
            if not obj["rows"]:
                continue
            d = pd.DataFrame(obj["rows"])
            assert d.year.eq(year).all()
            data.append(d)
    records = pd.concat(data, ignore_index=True)
    assert records.key.is_unique
    exact = records.speciesKey.map(taxa).eq(records.species)
    good = records[exact & records.basisOfRecord.eq("HUMAN_OBSERVATION") &
                   records.occurrenceStatus.isin(["PRESENT", "DETECTED"]) &
                   records.coordinateUncertaintyInMeters.gt(0) &
                   records.coordinateUncertaintyInMeters.le(1000)].rename(
                       columns={"key":"unit_id","decimalLongitude":"longitude","decimalLatitude":"latitude"})
    receipt = json.loads((coast / "source_receipt.json").read_text(encoding="utf-8"))
    for file, digest in receipt["extracted_sha256"].items():
        assert hashlib.sha256((coast / file).read_bytes()).hexdigest() == digest
    assert len(good) > 0, "No precise positive records; inspect GBIF fields"
    assignment = assign_points(good, prepare_land(gpd.read_file(coast / "GSHHS_h_L1.shp")),
                               gpd.read_file(coast / "GSHHS_h_L2.shp"))
    land = good.merge(assignment, on="unit_id", validate="one_to_one")
    land = land[land.land_class.isin(["island","mainland"])].copy()
    land["x1"] = np.floor(land.longitude).astype(int)
    land["y1"] = np.floor(land.latitude).astype(int)
    land["gx"] = np.floor((land.longitude+180)*20).astype(int)
    land["gy"] = np.floor((land.latitude+90)*20).astype(int)
    return land, dict(downloaded=len(records), precise_present=len(good), assigned=len(land),
                      input_sha256=hashes, query_scope=queries, dataset_keys=sorted(land.datasetKey.dropna().unique()))


def linked_prey(s, members, legacy, land, prey, local_moths=None):
    """Exact named prey -> flower visitor -> plant, with the same gecko recorded nearby."""
    receipt = json.loads((prey / "receipt.json").read_text(encoding="utf-8"))
    path = prey / "species_prey.csv"
    assert receipt["scanned"] == 24453954
    assert receipt["source_sha256"] == "670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3"
    assert receipt["prey_table_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    edges = pd.read_csv(path, keep_default_na=False)
    visits = legacy.loc[legacy.guild.eq("moth_family_visit"), ["insect_name","plant_name"]].drop_duplicates()
    keys = ["year","land_class","land_id","x1","y1"]
    observations = s[["key",*keys]].merge(land[[*keys,"species"]].drop_duplicates(), on=keys)
    counts = []
    for gecko in SPECIES.values():
        animal = edges[edges.gecko.eq(gecko)]
        for stage in ["any","adult"]:
            eligible = animal if stage == "any" else animal[animal.source_stage.eq("adult")]
            prey_names = set(eligible.moth) - {""}
            partners = visits[visits.insect_name.isin(prey_names)]
            observed = members[members.wcvp_accepted_name.isin(partners.plant_name)]
            documented = observed.merge(observations.loc[observations.species.eq(gecko)], on="key", validate="many_to_one")
            counts.append(dict(gecko=gecko, prey_stage=stage, relationship_rows=int(eligible.relationship_rows.sum()),
                               named_prey=len(prey_names), named_prey_with_flower_link=partners.insect_name.nunique(),
                               study_plant_names=partners.plant_name.nunique(), plants_on_primary_sites=observed[ID].nunique(),
                               plants_on_gecko_positive_sites=documented[ID].nunique(),
                               site_years_by_land=documented.groupby("land_class").key.nunique().to_dict()))
    chains = []
    if local_moths is not None:
        adult_links = edges[edges.source_stage.eq("adult") & edges.moth.isin(visits.insect_name)]
        assert set(adult_links.moth) == set(LINKED_MOTHS.values())
        for pair in adult_links[["gecko","moth"]].drop_duplicates().itertuples(index=False):
            plants = visits.loc[visits.insect_name.eq(pair.moth),"plant_name"]
            observed = members[members.wcvp_accepted_name.isin(plants)]
            for scale, dims in [("1deg",["x1","y1"]),("0.05deg",["gx","gy"])]:
                cols = ["year","land_class","land_id",*dims]
                sites = s.merge(observed[["key",ID]],on="key")
                same_gecko = land.loc[land.species.eq(pair.gecko),cols].drop_duplicates()
                sites = sites.merge(same_gecko,on=cols,validate="many_to_one")
                m = local_moths[local_moths.species.eq(pair.moth)]
                for tier in ("any_stage","adult_recorded"):
                    subset = m if tier == "any_stage" else m[m.lifeStage.fillna("").str.lower().isin(["adult","imago"])]
                    matched = sites.merge(subset[cols].drop_duplicates(),on=cols,validate="many_to_one")
                    chains.append(dict(gecko=pair.gecko,moth=pair.moth,gbif_stage=tier,scale=scale,
                                       plants=matched[ID].nunique(), years=sorted(map(int,matched.year.unique())),
                                       site_years_by_land=matched.groupby("land_class").key.nunique().to_dict()))
    return counts, chains


def multiguild_support(s, land, taxa):
    rows = []
    for name in taxa.values():
        detections = land[land.species.eq(name)]
        for scale, dims in [("1deg",["x1","y1"]),("0.05deg",["gx","gy"])]:
            keys = ["year","land_class","land_id",*dims]
            overlap = s.merge(detections[keys].drop_duplicates(), on=keys, validate="many_to_one")
            for land_class in ("island","mainland"):
                group = overlap[overlap.land_class.eq(land_class)]
                rows.append(dict(species=name, scale=scale, land_class=land_class,
                                 documented_site_years=len(group), years=group.year.nunique(),
                                 physical_landmasses=group.land_id.nunique(),
                                 sites_with_legacy_moth_plant=int(group.all_moth.gt(0).sum())))
    return rows


def run(inputs, signed, private, out, prey=None, multiguild=None, moths=None, write_blocks=None):
    s, digest, members, legacy = primary(inputs, signed)
    land, provenance = positive_records(private, inputs / "coast")
    parts = s.target_id.str.split("|", regex=False, expand=True)
    grid = parts[1].str.split("_", expand=True).astype(int)
    s["land_id"] = parts[2]
    s["gx"], s["gy"] = grid[0], grid[1]
    s["x1"], s["y1"] = s.gx // 20 - 180, s.gy // 20 - 90
    coarse = s.assign(x=s.x1//5*5, y=s.y1//5*5)
    blocks = {str(year):g[["x","y"]].drop_duplicates().values.tolist() for year,g in coarse.groupby("year")}
    one_degree = {str(year):g[["x1","y1"]].drop_duplicates().values.tolist() for year,g in s.groupby("year")}
    if write_blocks is not None:
        write_blocks.write_text(json.dumps(dict(target_hash=digest,blocks=blocks,one_degree=one_degree)),encoding="utf-8")
    for label, dims in [("same_1deg",["x1","y1"]), ("same_005deg",["gx","gy"])]:
        cols = ["year","land_id","land_class",*dims]
        s[label] = pd.MultiIndex.from_frame(s[cols]).isin(pd.MultiIndex.from_frame(land[cols]))
    assert (~s.same_005deg | s.same_1deg).all()
    positive = s[s.same_1deg].copy()
    shared = positive[positive.groupby("region_year").land_class.transform("nunique").eq(2)].copy()
    model = fit(shared,"all","moth",balance=True)
    tables = s.groupby(["year","land_class"]).agg(primary=("target_id","size"),
                 documented_1deg=("same_1deg","sum"),documented_005deg=("same_005deg","sum")).reset_index()
    result = dict(frozen_primary_hash=digest, gbif=provenance, support=tables.to_dict("records"),
                  landmasses_by_class=positive.groupby("land_class").land_id.nunique().to_dict(),
                  positive_same_1deg=len(positive), matched_year_region_counts=shared.groupby("land_class").size().to_dict(),
                  model_conditional_on_positive_documentation=model)
    if prey is not None:
        local_moths = None
        if moths is not None:
            local_moths, result["moth_gbif"] = positive_records(moths, inputs / "coast", LINKED_MOTHS)
        result["gecko_prey_legacy_plant_link"], result["three_way_local_support"] = linked_prey(
            s, members, legacy, land, prey, local_moths)
    if multiguild is not None:
        extra, source = positive_records(multiguild, inputs / "coast", MULTIGUILD, one_degree)
        result["additional_gbif"] = source
        result["additional_species_support"] = multiguild_support(s, extra, MULTIGUILD)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs","signed","private","out"):
        p.add_argument(name,type=Path)
    p.add_argument("--prey", type=Path)
    p.add_argument("--multiguild", type=Path)
    p.add_argument("--moths", type=Path)
    p.add_argument("--write-blocks", type=Path)
    a = p.parse_args()
    run(a.inputs,a.signed,a.private,a.out,a.prey,a.multiguild,a.moths,a.write_blocks)
