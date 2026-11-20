"""Attach source adjudication and moderator coverage without changing estimates."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
FILES = ["effect_catalogue.csv", "studies.csv", "site_context/catalogue_coordinate_coverage.csv", "site_context/plant_moisture_effect_support.csv", "experiments.csv"]

def build():
    effects, studies, sites, moisture, experiments = [pd.read_csv(ROOT / name) for name in FILES]
    assert effects.effect_id.is_unique and studies.study_id.is_unique
    assert sites.programme_id.is_unique and moisture.effect_id.is_unique
    assert set(effects.study_id) <= set(studies.study_id)
    evidence = effects.merge(studies.rename(columns={"programme_id": "study_programme_id", "status": "study_review_status"}), on="study_id", how="left", validate="many_to_one")
    assert evidence.programme_id.eq(evidence.study_programme_id).all()
    evidence = evidence.drop(columns="study_programme_id")
    evidence["primary_article_url"] = "https://doi.org/" + evidence.doi
    site_columns = ["programme_id", "coordinate_records", "records_with_numeric_coordinates", "has_site_coordinate", "site_ids", "coordinate_role"]
    evidence = evidence.merge(sites[site_columns], on="programme_id", how="left", validate="many_to_one")
    assert evidence.has_site_coordinate.notna().all()
    linked_moisture = moisture.loc[moisture.effect_id.isin(effects.effect_id)]
    comparison = linked_moisture.merge(effects[["effect_id", "lnRR", "sampling_variance", "programme_id"]], on="effect_id", suffixes=("_moisture", "_catalogue"), validate="one_to_one")
    for column in ["lnRR", "sampling_variance"]:
        assert np.allclose(comparison[column+"_moisture"], comparison[column+"_catalogue"], equal_nan=True)
    assert comparison.programme_id_moisture.eq(comparison.programme_id_catalogue).all()
    cols = ["effect_id", "source_species", "taxon_status", "EIV_MOIST", "EIV_MOIST_raw", "join_method", "species_affinity_not_site_condition"]
    evidence = evidence.merge(linked_moisture[cols], on="effect_id", how="left", validate="one_to_one")
    evidence["species_moisture_status"] = np.where(evidence.EIV_MOIST.notna(), "numeric_species_affinity_not_site_wetness", np.where(evidence.effect_id.isin(linked_moisture.effect_id), "endpoint_mapped_numeric_moisture_unavailable", "no_single_species_mapping_in_current_join"))
    evidence["experiment_row_present"] = evidence.experiment_id.isin(experiments.experiment_id)
    evidence["environmental_moderator_readiness"] = "not_established_by_coordinate_presence"
    evidence["final_eligibility"] = "not_adjudicated"
    assert len(evidence) == len(effects)
    assert evidence.effect_id.tolist() == effects.effect_id.tolist()
    for column in ["lnRR", "sampling_variance", "mean_lit", "mean_dark", "n_lit", "n_dark"]:
        assert np.allclose(evidence[column], effects[column], equal_nan=True)
    evidence.to_csv(ROOT / "catalogue_evidence_view.csv", index=False)
    summary = dict(effect_rows=len(evidence), studies=evidence.study_id.nunique(), programmes=evidence.programme_id.nunique(),
        study_doi_linked_rows=int(evidence.doi.notna().sum()),
        rows_with_experiment_id=int(evidence.experiment_id.notna().sum()),
        rows_with_experiment_table_record=int(evidence.experiment_row_present.sum()),
        programme_labels_with_coordinates=int(evidence.loc[evidence.has_site_coordinate, "programme_id"].nunique()),
        rows_with_numeric_species_moisture=int(evidence.EIV_MOIST.notna().sum()),
        species_with_numeric_moisture=int(evidence.loc[evidence.EIV_MOIST.notna(), "source_species"].nunique()),
        programmes_with_numeric_species_moisture=int(evidence.loc[evidence.EIV_MOIST.notna(), "programme_id"].nunique()),
        study_status_counts=evidence.groupby("study_review_status").size().to_dict(),
        effect_rows_added=0, effect_values_changed=0,
        input_sha256={name:hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in FILES},
        limitation="Coverage is not approval or statistical independence. Unmapped EIV is not dry habitat. Coordinates do not establish measured site wetness. Structured experiment records remain incomplete.")
    (ROOT / "catalogue_evidence_audit.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k not in ["study_status_counts", "input_sha256"]}, indent=2))

if __name__ == "__main__":
    build()
