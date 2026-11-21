#!/usr/bin/env python3
"""Sanders 2021: source-exact 23-paper preliminary A process triage.

Original 1304-effect CSV/126 Study IDs fetched via original author Zenodo
MD5/SHA; Stage A broad arthropod pool is 23 source Study IDs, 152 Keep effects.
Published original 'Study' DOES NOT mean 23 independently randomized lamps.
The original paper/DOI membership is independently checked against the
manually labelled preliminary primary-study mechanism classes.

No one may fit A 'water-side vs terrestrial predator abundance' as though
fish behaviour, marine amphipod attraction, moth life history, spider prey
interception and stream insect drift were a common abundance response.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from scripts.analysis.alan_sanders2020_keep_delete_stageAB_primary_screen_20261119 import (
    audit_keep_delete,
)
from scripts.analysis.alan_sanders2020_stageAB_source_census_20261119 import (
    original_source_rows,
)
from scripts.acquisition.alan_sanders2020_source_126paper_Dryad_preflight_20261119 import (
    zenodo_original_author_mirror,
)

SOURCE=Path("config/alan_sanders2020_23paper_stageA_original_endpoint_title_triage_20261119.json")
SCHEMA="alan_sanders126_author_source_23_stageA_title_primary_process_review_v1"
FAMILIES={
    "trophic_predation_parasitism":7,
    "aquatic_behaviour_not_emergence":5,
    "aerial_light_attraction_activity":5,
    "other_noncomparable_behaviour_reproduction":5,
    "aquatic_drift_source_NOT_equivalent_to_emergence":1,
}
SOURCE_DIGEST="d11d1e29c05826f2b4a5c1c2d63b9c195dc2effbedd1ff531afa7b6bdf497693"


def audit(register,author_bibliography,strict=True):
    if register.get("schema")!=SCHEMA or register.get("source_original_CSV_SHA256")!=SOURCE_DIGEST:
        raise ValueError("Original source identity changed")
    if register.get("distinct_source_study_records")!=23:
        raise ValueError("23 source Study ID candidates changed")
    if not register["includes_title_screen_only_not_systematic_search"]:
        raise ValueError("Preliminary labels cannot count as systematic search")
    if not register["no_automated_habitat_determination_from_plant_or_insect_taxonomy"]:
        raise ValueError("Wetland habitat classification cannot be inferred from taxon")
    rows=register["source_study_entries"]
    if len(rows)!=23 or len(author_bibliography)!=23:
        raise ValueError("Wrong original Stage A preliminary source Study count")
    canonical=lambda x:tuple(sorted(x))
    original=[canonical(x["source_published_DOIs"]) for x in author_bibliography]
    manual=[canonical(x["original_source_DOIs"]) for x in rows]
    if len(set(manual))!=23 or sorted(manual)!=sorted(original):
        raise ValueError("DOI list does not exactly match original author source Study cluster")
    families=Counter(x["preliminary_outcome_family_for_screen_only"] for x in rows)
    if strict and dict(families)!=FAMILIES:
        raise ValueError("Original study source stage classification changed; recheck")
    if register["title_screen_only_families_counts_not_meta_ready"]!=dict(families):
        raise ValueError("Initial A family counts inconsistent with source records")
    if any(x.get("final_PRISMA_fulltext_decision")!="NOT_SCREENED" for x in rows):
        raise ValueError("Title/abstract labels cannot be automatically counted as final systematic inclusion")
    if any(x.get("primary_illumination_programme_id_verified") is not False for x in rows):
        raise ValueError("Physical source experiment identity must be separately audited")
    if any(x.get("original_research_light_LRR_SE_verified") is not False for x in rows):
        raise ValueError("No original source illumination-cluster-level lnRR variance available in this title screen")
    if register.get("strict_aquatic_emergence_moderator_meta_eligible_k")!="UNKNOWN":
        raise ValueError("Strict wetland emerging insect meta denominator remains unknown")
    if register.get("strict_predator_abundance_moderator_meta_eligible_k")!="UNKNOWN":
        raise ValueError("True predator abundance independent study count remains unknown")
    conflict=sum(len(x["source_published_DOIs"])>1 for x in author_bibliography)
    if conflict!=1:
        raise ValueError("Original Sanders two conflicting DOI strings for one Study must be quarantined")
    result={
        "schema":"alan_sanders_original23_stageA_separate_biological_processes_title_screen_gate_v1",
        "author_original_Sanders_CSV_SHA":SOURCE_DIGEST,
        "original_broad_stageA_Keep_StudyID_count":len(rows),
        "original_study_with_source_DOI_identity_conflict":conflict,
        "preliminary_outcome_family_counts_not_meta_sample":dict(sorted(families.items())),
        "single_original_stream_drift_preliminary_paper_count":families["aquatic_drift_source_NOT_equivalent_to_emergence"],
        "aquatic_mobile_behaviour_is_not_emergent_insect_abundance":True,
        "total_arthropod_predator_abundance_meta_ready_k":None,
        "riparian_terrestrial_modification_meta_ready_k":None,
        "independent_physical_lighting_installations_verified_for_all_23":False,
        "study_habitats_from_primary_full_text_fully_verified":False,
        "WoS_or_Scopus_PRISMA_search_paper_count":None,
        "no_new_global_lnRR_or_predation_causal_meta_fitted":True,
        "published_author_Source_StudyID_is_not_independent_light_programme":True,
        "do_not_inflate_site_levels_from_multiple_Penryn_Lizard_sources":True,
    }
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--register",type=Path,default=SOURCE)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    data=json.loads(args.register.read_text())
    source,_=zenodo_original_author_mirror()
    original_rows,base=original_source_rows(source,strict=True)
    paper_source=audit_keep_delete(original_rows,base,strict=True)
    pubs=paper_source["candidate_publication_bibliography_for_fulltext_screening"][
         "A_arthropod_abundance_predation_BROAD_SOURCE_SCREEN"]
    result=audit(data,pubs,strict=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("SANDERS23_BROAD_STAGEA_ORIGINAL_SOURCE_PRIMARY_PROCESS_FAMILIES",
          result["preliminary_outcome_family_counts_not_meta_sample"])
    print("SANDERS23_SOURCE_DOI_IDENTITY_CONFLICT",result["original_study_with_source_DOI_identity_conflict"])
    print("SANDERS23_ORIGINAL_STRICT_WETLAND_META_ESSENTIAL_GATE",
          result["riparian_terrestrial_modification_meta_ready_k"])

if __name__=="__main__":
    main()
