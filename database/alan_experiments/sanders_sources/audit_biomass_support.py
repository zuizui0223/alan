"""Audit field biomass sampling support; do not estimate selective effects."""
from pathlib import Path
import zipfile
import hashlib
import json
import pandas as pd

OUT = Path(__file__).resolve().parent
ARCHIVE = OUT.parents[3] / "gbif_annual/alan_experiment_sources/sanders2018/data.zip"
with zipfile.ZipFile(ARCHIVE) as z:
    population = pd.read_csv(z.open("data/ALANIntensityPopulationdynamics.csv"))
    design = population[["Cage", "Block", "Treat.lux"]].drop_duplicates()
    assert len(design) == 48 and design.Cage.is_unique
    assert design.groupby(["Block", "Treat.lux"]).size().eq(1).all()
    units, arms, pairs, conflicts, hashes = [], [], [], [], {}
    for species in ["bean", "barley"]:
        member = f"data/ALANIntensityPlantbiomass_field_{species}.csv"
        raw = z.read(member)
        hashes[member] = hashlib.sha256(raw).hexdigest()
        biomass = pd.read_csv(z.open(member))
        assert biomass.Cage.is_unique
        join = design.merge(biomass, on="Cage", how="left", suffixes=("", "_biomass"), validate="one_to_one", indicator=True)
        found = join._merge.eq("both")
        assert set(biomass.Cage) <= set(design.Cage)
        assert join.loc[found, "Block"].eq(join.loc[found, "Block_biomass"]).all()
        mismatch = found & ~join["Treat.lux"].eq(join["Treat.lux_biomass"])
        for _, row in join.loc[mismatch].iterrows():
            conflicts.append(dict(species=species, Cage=int(row.Cage), Block=int(row.Block), population_light=row["Treat.lux"], biomass_light=row["Treat.lux_biomass"]))
        assert join.loc[found, "Plantbiomass.g"].notna().all()
        join["species"] = species
        join["recorded"] = found
        units.append(join[["species", "Cage", "Block", "Treat.lux", "recorded", "Treat.lux_biomass", "Plantbiomass.g"]])
        for dose, subset in join.groupby("Treat.lux"):
            arms.append(dict(species=species, treatment=dose, planned_cages=len(subset), recorded_cages=int(subset.recorded.sum()), missing_cages=int((~subset.recorded).sum()), basis="population_design_ID_link_not_corrected_biomass_allocation"))
        dark = set(join.loc[join.recorded & join["Treat.lux_biomass"].eq("c"), "Block"])
        for dose in ["0.1", "1", "5", "10", "20", "50", "100"]:
            lit = set(join.loc[join.recorded & join["Treat.lux_biomass"].eq(dose), "Block"])
            complete = sorted(dark & lit)
            pairs.append(dict(species=species, light_lux=dose, complete_block_pairs=len(complete), block_ids=";".join(map(str, complete)), extraction_status="biomass_labels_only_allocation_conflict_and_missingness_hold"))
units = pd.concat(units, ignore_index=True)
units.to_csv(OUT / "biomass_sampling_units.csv", index=False)
pd.DataFrame(arms).to_csv(OUT / "biomass_arm_support.csv", index=False)
pd.DataFrame(conflicts).to_csv(OUT / "biomass_allocation_conflicts.csv", index=False)
pd.DataFrame(pairs).to_csv(OUT / "biomass_paired_support.csv", index=False)
common = units.pivot(index="Cage", columns="species", values="recorded")
audit = dict(design_cages=48, recorded_by_species=units.groupby("species").recorded.sum().astype(int).to_dict(),
             both_species_recorded=int(common.all(axis=1).sum()), neither_species_recorded=int((~common.any(axis=1)).sum()),
             archive_sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(), member_sha256=hashes,
             allocation_matches_population_data=not conflicts, allocation_conflict_records=len(conflicts), imputed_values=0, effect_estimates_added=0,
             source_reason="Author data description: missing biomass due to loss of labels; randomness not established.")
(OUT / "biomass_support_audit.json").write_text(json.dumps(audit, indent=2)+"\n", encoding="utf-8")
print(pd.DataFrame(pairs).to_string(index=False))
print(json.dumps(audit, indent=2))
