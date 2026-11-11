"""Validate tables saved by run 42952223131 and retain aggregate source evidence."""
import hashlib
import json
import sys
from pathlib import Path
import pandas as pd


def audit(root, out):
    r = json.loads((root/"source_receipt.json").read_text(encoding="utf-8"))
    assert r["status"] == "FULL_GLOBI_V011_VERIFIED" and r["all_original_interaction_rows"] == 24453954
    assert r["source_sha256"] == "670c34b667f0a13d2bfd99edf7a3aa1518b92fa5493614f9dbbcaa0c6821e0c3"
    assert r["limitations"]["one_species_dual_role_not_same_place_or_same_time"]
    assert sum(r["moth_prey_stage_counts"].values()) == r["qualified_role_counts"]["moth_predation"]
    guilds = pd.read_csv(root/"guild_summary.csv")
    assert len(guilds) == 36 and guilds.to_dict("records") == r["guild_summary"]
    animals = pd.read_csv(root/"animal_species_role_flags.csv.gz")
    plants = pd.read_csv(root/"animal_flower_plant_positive_flags.csv.gz")
    assert not animals.duplicated(["scope","animal_group","animal_species"]).any()
    assert not plants.duplicated(["scope","animal_group","plant_species"]).any()
    for g in guilds.itertuples():
        a = animals[(animals.scope==g.scope) & (animals.animal_group==g.predator_group)]
        p = plants[(plants.scope==g.scope) & (plants.animal_group==g.predator_group)]
        flower = a.flower_visit.eq(1) | a.reported_pollinates.eq(1)
        assert len(a) == g.animal_species_with_any_record
        assert a.moth_predation.sum() == g.animal_species_moth_predation_recorded
        assert a.flower_visit.sum() == g.animal_species_flower_visit_recorded
        assert a.reported_pollinates.sum() == g.animal_species_pollinates_relation_recorded
        assert flower.sum() == g.animal_species_flower_positive_either
        assert (flower & a.moth_predation.eq(1)).sum() == g.animal_species_both_moth_predation_and_flower_positive
        assert (a.reported_pollinates.eq(1) & a.moth_predation.eq(1)).sum() == g.animal_species_both_moth_predation_and_pollinates_relation
        assert a.plant_dispersal.sum() == g.animal_species_plant_dispersal_recorded
        assert (p.flower_visit.eq(1) | p.reported_pollinates.eq(1)).sum() == g.plants_with_any_flower_visit_from_group
        pairs = [e["distinct_named_species_pairs"] or 0 for e in r["predator_predator_group_links"]
                 if e["scope"]==g.scope and e["predator_group"]==g.predator_group]
        assert sum(pairs)==g.focused_predator_of_predator_species_edges
    report = dict(run_id=42952223131, workflow_conclusion="failure in post-scan key assertion",
                 artifact_files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.glob("*") if p.is_file()},
                 full_input_sha256=r["source_sha256"],all_interaction_rows=r["all_original_interaction_rows"],
                 n_group_scope_rows_verified=len(guilds), stage=r["moth_prey_stage_counts"],
                 all_sources=guilds[guilds.scope.eq("all")].to_dict("records"))
    out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    audit(Path(sys.argv[1]),Path(sys.argv[2]))
