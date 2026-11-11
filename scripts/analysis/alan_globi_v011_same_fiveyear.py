"""Same-design relabelling of the frozen five-year GloBI moth interaction."""
import argparse
import hashlib
import json
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import alan_globi_crossland_light_robust_20261111 as original

ROOT = Path(__file__).resolve().parents[2]
YEARS = (2019, 2021, 2023, 2024, 2025)
ID = "wcvp_accepted_plant_name_id"
TERM = "beta_island_minus_mainland_slope_pp_per_log1p"


def read_one(folder, pattern, **kwargs):
    path, = folder.rglob(pattern)
    return pd.read_csv(path, low_memory=False, **kwargs)


def environment(folder):
    frames = []
    for year in YEARS:
        if year < 2024:
            e = read_one(folder / f"g{year}", "values.csv")
            p = read_one(folder / f"p{year}", "*popdensity.csv")
            pop = "population_density_per_km2"
        elif year == 2024:
            e = read_one(folder / "c2024", "*combined_values.csv")
            p = read_one(folder / "p2024", "*fullcohort_values.csv")
            pop = "population_UC_R2024B_density_per_km2"
        else:
            e = read_one(folder / "v2025", "values.csv")
            e["pop_density"] = e.population_density
        if year < 2025:
            e = e.merge(p[["target_id", pop]], validate="one_to_one", on="target_id")
            e["pop_density"] = e[pop]
            # Historical inputs include mixed coast/water cells in the ambiguity check.
            geo = read_one(folder / f"g{year}", "*full_photo_landmass_source.csv")
            e = e[e.target_id.isin(geo.loc[~geo.cell_id.duplicated(False), "target_id"])].copy()
        # The archived 2025 numeric population contains known-land targets only.
        e["key"] = e.target_id.str.rsplit("|", n=1).str[0]
        e = e[~e.key.duplicated(False) & e.viirs_valid_months.ge(3)].copy()
        e["year"] = year
        e["region5"] = (np.floor((e.query_longitude + 180) / 5).astype(int).astype(str)
                        + "_" + np.floor((e.query_latitude + 90) / 5).astype(int).astype(str))
        frames.append(e[["target_id", "key", "year", "land_class", "viirs_median", "pop_density", "region5"]])
    return pd.concat(frames, ignore_index=True)


def taxonomy(inputs, old, new):
    flags = pd.DataFrame(index=sorted(old | new))
    flags.index.name = "plant_name"
    flags["old"] = flags.index.isin(old)
    flags["new"] = flags.index.isin(new)
    with zipfile.ZipFile(inputs / "wcvp.zip") as archive:
        w = pd.read_csv(archive.open("wcvp_names.csv"), sep="|", dtype=str, keep_default_na=False,
                        usecols=["taxon_name", "taxon_rank", "taxon_status", "plant_name_id", "accepted_plant_name_id"])
    species = w[w.taxon_rank.eq("Species")]
    accepted = species[species.taxon_status.eq("Accepted")].set_index("plant_name_id").taxon_name
    assert accepted.index.is_unique
    set_ids = set(accepted.index)
    flags["resolution"] = "no_exact_species_name"
    flags["accepted_id"] = ""
    flags["accepted_name"] = ""
    flags["accepted_priority_id"] = ""
    for name, candidates in species[species.taxon_name.isin(flags.index)].groupby("taxon_name"):
        targets = set(candidates.accepted_plant_name_id)
        direct = set(candidates.loc[candidates.taxon_status.eq("Accepted"), "plant_name_id"])
        status_ok = candidates.taxon_status.isin(["Accepted", "Synonym", "Homotypic Synonym", "Heterotypic Synonym"]).all()
        if len(targets) == 1 and targets <= set_ids and status_ok:
            target, = targets
            flags.loc[name, ["resolution", "accepted_id", "accepted_name", "accepted_priority_id"]] = [
                "resolved", target, accepted[target], target]
        else:
            flags.loc[name, "resolution"] = "ambiguous_or_unresolved"
            if len(direct) == 1:
                flags.loc[name, "accepted_priority_id"] = next(iter(direct))
    return flags


def site_species(signed, env):
    cols = ["photo_id", "source_subarm", "observed_on", "longitude", "latitude", "inat_taxon_id",
            "wc2_bio05_10m", "wc2_bio12_10m", "wc2_elev_2p5m"]
    p = read_one(signed / "canonical", "canonical*.gz", usecols=cols)
    p = p[p.source_subarm.eq("historical_main") & p.observed_on.str[:4].isin(map(str, YEARS))
          & p.latitude.between(-89.999, 89.999) & p.longitude.between(-179.999, 179.999)].copy()
    p["key"] = (p.observed_on.str[:4] + "|" + np.floor((p.longitude + 180) * 20).astype(int).astype(str)
                + "_" + np.floor((p.latitude + 90) * 20).astype(int).astype(str))
    p = p[p.key.isin(env.key)]
    effort = p.groupby("key").agg(photo_attempts=("photo_id", "size"), bio5=("wc2_bio05_10m", "mean"),
                                  rain=("wc2_bio12_10m", "mean"), elevation=("wc2_elev_2p5m", "mean"))
    w = read_one(signed / "taxonomy", "*resolution.csv.gz", dtype={ID: str})
    names = read_one(signed / "taxonomy", "*distribution*.gz", usecols=[ID, "wcvp_accepted_name"],
                     dtype={ID: str}).drop_duplicates()
    w = w[w.wcvp_match_state.eq("unique_exact_name_to_accepted_id")]
    sp = p[["key", "inat_taxon_id"]].merge(w[["inat_taxon_id", ID]], on="inat_taxon_id", validate="many_to_one")
    sp = sp[["key", ID]].drop_duplicates().merge(names, on=ID, how="left", validate="many_to_one")
    return env.merge(effort, on="key", validate="one_to_one"), sp


def delta_sets(old, new):
    return dict(old=len(old), new=len(new), retained=len(old & new), added=len(new-old), lost=len(old-new))


def run(inputs, signed, output):
    manifest = json.loads((ROOT / "config/alan_v011_fiveyear_input_sha256.json").read_text(encoding="utf-8"))
    for rel, expected in manifest.items():
        prefix, tail = rel.split("/", 1)
        path = (inputs if prefix == "inputs" else signed) / tail
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"SHA mismatch: {path.name}")
    old = read_one(inputs / "old", "GloBI_stage2_signed_extension.csv.gz")
    new = read_one(inputs, "GloBI_v011_positive_species_flags.csv.gz")
    old_names = set(old.loc[old.guild.eq("moth_family_visit"), "plant_name"])
    new_names = set(new.loc[new.all_moth_any_family.eq(1), "plant_name"])
    flags = taxonomy(inputs, old_names, new_names)
    sites, sp = site_species(signed, environment(inputs))
    study = old.study_title.fillna("").str.lower()
    for scope, mask in {"all": np.ones(len(old), dtype=bool), "non_inat": ~study.str.contains("inaturalist"),
                        "non_inat_trophich": ~study.str.contains("inaturalist|trophich"),
                        "non_inat_visits": ~study.str.contains("inaturalist") & old.relation.eq("visitsFlowersOf")}.items():
        for guild in ["moth", "butterfly"]:
            names = old.loc[mask & old.guild.eq(guild + "_family_visit"), "plant_name"]
            sp[f"{scope}_{guild}"] = sp.wcvp_accepted_name.isin(names)
    old_ids = set(flags.loc[flags.old, "accepted_id"]) - {""}
    new_ids = set(flags.loc[flags.new, "accepted_id"]) - {""}
    priority = set(flags.loc[flags.new, "accepted_priority_id"]) - {""}
    sp["old_resolved"] = sp[ID].isin(old_ids)
    sp["updated"] = sp[ID].isin(new_ids)
    sp["priority"] = sp[ID].isin(priority)
    sp["new_exact"] = sp.wcvp_accepted_name.isin(new_names)
    counts = sp.drop(columns=[ID, "wcvp_accepted_name"]).groupby("key").sum()
    sites = sites.merge(counts.assign(n_species=sp.groupby("key").size()), on="key", validate="one_to_one")
    with tempfile.TemporaryDirectory() as tmp:
        source = Path(tmp) / "sites.csv.gz"
        sites.to_csv(source, index=False)
        frame = original.prepare(original.cleaned(source), 3, 2)
    complete = frame.dropna(subset=["bio5", "lograin", "elevation_km", "logpop", "loglight", "logspecies", "logphotos"])
    complete = complete[complete.groupby("region_year").land_class.transform("nunique").eq(2)].copy()
    legacy = original.fit(complete, "all", "moth", balance=True)
    reference = json.loads((ROOT / "results/alan_globi_crossland_robust_light_20261111.json").read_text(encoding="utf-8"))
    ref = reference["result_pp_per_log1p_light"]["broad_all_moth"]
    np.testing.assert_allclose([legacy[TERM], *legacy["95CI_cluster_5deg"]], [ref["beta"], *ref["CI95"]], atol=1e-10, rtol=0)
    assert (len(complete), legacy["n_island_site_years"], legacy["source_5degree_clusters"], legacy["year_region_groups"]) == (10698,933,78,195)
    yearly = []
    annual = json.loads((ROOT / "results/alan_globi_allrecord_year_dynamics_20261111_source_verified.json").read_text(encoding="utf-8"))
    for ref_year in annual["all_GloBI_moth_primary"]["per_year"]:
        subset = complete[complete.year.eq(ref_year["year"])]
        fit = original.fit(subset, "all", "moth", balance=True)
        np.testing.assert_allclose([fit[TERM], *fit["95CI_cluster_5deg"]], [ref_year["beta"], *ref_year["ci95"]], atol=1e-8, rtol=0)
        yearly.append({"year": ref_year["year"], "old": fit,
                       "updated": original.fit(subset.assign(all_moth=subset.updated), "all", "moth", balance=True)})
    fits = {"published": legacy}
    for label in ("old_resolved", "updated", "priority", "new_exact"):
        fits[label] = original.fit(complete.assign(all_moth=complete[label]), "all", "moth", balance=True)
    diff = complete.updated - complete.all_moth
    fits["updated_minus_published"] = original.fit(complete.assign(all_moth=complete.updated, all_butterfly=complete.all_moth),
                                                   "all", "moth_minus_butterfly", balance=True)
    fits["updated_minus_old_resolved"] = original.fit(
        complete.assign(all_moth=complete.updated, all_butterfly=complete.old_resolved),
        "all", "moth_minus_butterfly", balance=True)
    np.testing.assert_allclose(fits["updated_minus_published"][TERM], fits["updated"][TERM] - legacy[TERM], atol=1e-10)
    used = sp[sp.key.isin(complete.key)]
    photographed = set(used[ID])
    baseline_ids = set(used.loc[used.all_moth, ID])
    result = {"fits": fits, "yearly": yearly, "global_wcvp_sets": delta_sets(old_ids, new_ids),
              "site_species_harmonized": delta_sets(old_ids & photographed, new_ids & photographed),
              "site_species_as_published_vs_updated": delta_sets(baseline_ids, new_ids & photographed),
              "dated_site_species_units": delta_sets(set(zip(used.loc[used.all_moth,"key"], used.loc[used.all_moth,ID])),
                                                     set(zip(used.loc[used.updated,"key"], used.loc[used.updated,ID]))),
              "n_sites_increased": int((diff>0).sum()), "n_sites_decreased": int((diff<0).sum()),
              "n_sites_same_count": int((diff==0).sum()), "n_species_denominator": len(photographed),
              "target_key_set_sha256": hashlib.sha256(("\n".join(sorted(complete.target_id))+"\n").encode()).hexdigest(),
              "unresolved": {col: flags.loc[flags[col],"resolution"].value_counts().to_dict() for col in ("old","new")}}
    output.mkdir(exist_ok=True, parents=True)
    flags.to_csv(output / "moth_name_crosswalk.csv")
    species_table = flags[flags.accepted_id.ne("")][["accepted_id", "accepted_name"]].drop_duplicates()
    for col, ids in [("legacy_wcvp",old_ids), ("v011_wcvp",new_ids)]:
        species_table[col] = species_table.accepted_id.isin(ids)
    species_table.sort_values("accepted_name").to_csv(output / "source_moth_species.csv", index=False)
    (output / "comparison.json").write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["inputs","signed","output"]:
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    run(args.inputs, args.signed, args.output)
