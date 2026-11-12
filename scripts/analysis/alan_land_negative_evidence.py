"""Link the screened species to archived adult-prey chains and genetic mating estimates."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from alan_globi_v011_same_fiveyear import taxonomy
from alan_gecko_primary_fiveyear import LINKED_MOTHS, positive_records, primary

ROOT = Path(__file__).resolve().parents[2] / "results"
OUT = ROOT / "species_interpretation"


def gecko_matches(inputs, signed, private, chains):
    sites, _, members, _ = primary(inputs, signed)
    target = sites.target_id.str.split("|", regex=False, expand=True)
    grid = target[1].str.split("_", expand=True).astype(int)
    sites = sites.assign(land_id=target[2], gx=grid[0], gy=grid[1])
    predators, _ = positive_records(private / "private", inputs / "coast")
    moths, _ = positive_records(private / "moths_private", inputs / "coast", LINKED_MOTHS)
    moths = moths[moths.lifeStage.fillna("").str.lower().isin(["adult", "imago"])]
    keys = ["year", "land_id", "land_class", "gx", "gy"]
    expected = json.loads((ROOT / "alan_gecko_primary_fiveyear_20261111.json").read_text())
    legacy = pd.read_csv(inputs / "old/GloBI_stage2_signed_extension.csv.gz")
    for row in expected["three_way_local_support"]:
        if row["scale"] != "0.05deg" or row["gbif_stage"] != "adult_recorded":
            continue
        predator, moth = row["gecko"], row["moth"]
        visited = legacy.loc[legacy.guild.eq("moth_family_visit") & legacy.insect_name.eq(moth), "plant_name"]
        linked = members[members.wcvp_accepted_name.isin(visited)].merge(sites[["key", *keys]], on="key")
        for frame, name in [(predators, predator), (moths, moth)]:
            linked = linked.merge(frame.loc[frame.species.eq(name),keys].drop_duplicates(), on=keys)
        assert linked.groupby("land_class").key.nunique().to_dict() == row["site_years_by_land"]
        mask = chains.predator.eq(predator) & chains.moth.eq(moth)
        for side in ["island", "mainland"]:
            observed = linked[linked.land_class.eq(side)].groupby("wcvp_accepted_name").key.nunique()
            chains.loc[mask,"gecko_005deg_site_years_"+side] = chains.loc[mask,"plant"].map(observed).fillna(0)
    return chains


def run(inputs, mating, signed, private):
    plants = set(pd.read_csv(OUT / "land_negative_species.csv").plant)
    assert hashlib.md5(mating.read_bytes()).hexdigest() == "c6b816501a98d336b6fa43fa60ae7ca4"
    rates = pd.read_csv(mating, encoding="latin1")
    rates["source_name"] = rates.genus.str.strip() + " " + rates.species.str.strip()
    names = set(rates.source_name.dropna())
    crosswalk = taxonomy(inputs, names, set())
    accepted = crosswalk.loc[crosswalk.resolution.eq("resolved"), "accepted_name"]
    rates["plant"] = rates.source_name.map(accepted)
    exact = rates.source_name.isin(plants)
    rates.loc[exact, "plant"] = rates.loc[exact, "source_name"]
    rates = rates[rates.plant.isin(plants) & rates["mean.tm"].notna()].copy()
    rates["selfing_from_tm"] = (1 - rates["mean.tm"]).where(rates["mean.tm"].between(0, 1))
    rates[["plant", "source_name", "year", "reference", "latitude", "hemisphere", "mean.tm",
           "selfing_from_tm"]].sort_values(["plant","year"]).to_csv(OUT / "land_negative_mating_records.csv", index=False)
    geckos = pd.read_csv(ROOT / "alan_two_house_gecko_species_prey_20261111.csv").rename(columns={"gecko":"predator"})
    bugs = pd.read_csv(ROOT / "arthropod_moth_source_20261111/species_prey.csv")
    prey = pd.concat([geckos, bugs])
    prey = prey.loc[prey.source_stage.eq("adult"), ["predator", "moth"]].dropna().drop_duplicates()
    visits = pd.read_csv(inputs / "old/GloBI_stage2_signed_extension.csv.gz")
    visits = visits.loc[visits.guild.eq("moth_family_visit"), ["insect_name", "plant_name"]].drop_duplicates()
    chains = prey.merge(visits.rename(columns={"insect_name":"moth", "plant_name":"plant"}), on="moth")
    dual = pd.read_csv(ROOT / "dual_role_80/island_cohort_reviewed_links.csv")
    confirmed = dual.loc[dual.assessment.eq("predation_recorded"), ["predator","moth","plant"]].assign(reviewed_same_plant=True)
    chains = chains.merge(confirmed, on=["predator","moth","plant"],how="outer")
    chains = chains[chains.plant.isin(plants)].copy()
    archived = json.loads((ROOT / "alan_arthropod_adult_moth_flora_20261111.json").read_text())
    for side in ["island", "mainland"]:
        by_pair = {(r["predator"],r["moth"]):set(r["plants_at_matched_sites"][side]) for r in archived["pairs"]}
        chains["in_archived_arthropod_005deg_"+side] = [
            (r.plant in by_pair[r.predator,r.moth]) if (r.predator,r.moth) in by_pair else None
            for r in chains.itertuples()]
    chains = gecko_matches(inputs, signed, private, chains)
    chains.sort_values(["plant","predator","moth"]).to_csv(OUT / "land_negative_predator_chains.csv",index=False)
    profiles = pd.read_csv(OUT / "land_negative_species.csv")
    summary = rates.groupby("plant").agg(genetic_tm_records=("mean.tm", "count"),
                                         selfing_min=("selfing_from_tm", "min"),
                                         selfing_max=("selfing_from_tm", "max"))
    profiles = profiles.join(summary, on="plant")
    profiles["genetic_tm_records"] = profiles.genetic_tm_records.fillna(0).astype(int)
    profiles["adult_prey_chains"] = profiles.plant.map(chains.groupby("plant").size()).fillna(0).astype(int)
    profiles["predator_species"] = profiles.plant.map(chains.groupby("plant").predator.agg(
        lambda names: " | ".join(sorted(set(names)))))
    for side in ["island", "mainland"]:
        supported = chains["gecko_005deg_site_years_" + side].gt(0) | chains[
            "in_archived_arthropod_005deg_" + side].eq(True)
        profiles["supported_005deg_chains_" + side] = profiles.plant.map(
            chains[supported].groupby("plant").size()).fillna(0).astype(int)
    profiles.to_csv(OUT / "land_negative_species_evidence.csv", index=False)
    print("Mating records/species",len(rates),rates.plant.nunique(),"adult-prey chains/species",len(chains),chains.plant.nunique())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs",type=Path)
    parser.add_argument("mating",type=Path)
    parser.add_argument("signed",type=Path)
    parser.add_argument("private",type=Path)
    args = parser.parse_args()
    run(args.inputs,args.mating,args.signed,args.private)
