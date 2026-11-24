#!/usr/bin/env python3
"""Adjudicate the 30 high-priority DOI *strings* from pinned 62-repo scan.

These 30 are unmatched to narrowly selected ALAN DOI registries, not
necessarily new papers. Distinguish archival data persistent identifiers
and articles, then preserve physical-programme overlap and extraction gates.
Titles/designs for ten article DOI strings checked against original ALAN
source audit notes and publisher/article pages on 2026-11-24.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

DATA_PREFIXES=(
    "10.16904/envidat.",
    "10.5061/dryad.",
    "10.5281/zenodo.",
    "10.5285/",
    "10.6084/m9.figshare.",
)
ARTICLES={
    "10.1002/eap.1821":(
        "ALAN_BIOTIC_MANIPULATED_REPORT",
        "Sullivan2013_wetland_three_independent_blocks",
        "Riparian wetland arthropod energy source, original-source reanalysis already in ALAN PR101; 3 blocks, not a novel source programme.",
        "SOURCE_PROGRAMME_ALREADY_AUDITED_DRAFT_NO_NEW_K",
    ),
    "10.1002/ecs2.3227":(
        "ALAN_BIOTIC_MANIPULATED_REPORT",
        "Blandy2015_2016_four_randomized_plot_pairs_Firebaugh",
        "Two-year grassland canopy arthropod abundance and trophic structure; source rounded paired t3 report, not independent 2015+2016 studies.",
        "SOURCE_PAIRED_RAW_DENSITY_HOLD_LOGRR_VARIANCE",
    ),
    "10.1007/s00027-025-01161-7":(
        "ALAN_BIOTIC_MANIPULATED_REPORT",
        "Westhavelland2012_lit_one_site_dark_one_site",
        "Manfrin2025 four emergence/12 flying traps per site; time samples are pseudo-replicates; same physical installation as 2017.",
        "NONREPLICATED_LIT_DARK_SITE_HOLD_OUT_OF_INDEPENDENT_META",
    ),
    "10.1098/rsbl.2019.0877":(
        "NO_EXPERIMENTAL_ALAN_POLLEN_NETWORK",
        "Walton2020_nocturnal_moth_pollen_network",
        "Wildflower moth pollen-carrying interaction network; no randomized light/dark manipulation.",
        "BACKGROUND_ONLY_NO_ALAN_LIGHT_EFFECT",
    ),
    "10.1111/een.13446":(
        "NO_EXPERIMENTAL_ALAN_POLLEN_NETWORK",
        "Atwater2025_moth_pollen_network",
        "Flower-settling moth interaction network; not an ALAN light experiment.",
        "BACKGROUND_ONLY_NO_ALAN_LIGHT_EFFECT",
    ),
    "10.1111/jvs.13168":(
        "PLANT_EIV_TRAIT_REFERENCE",
        "Tichy2023_Ellenberg_values",
        "Ellenberg-type indicator values used to describe plant wetness, not an intervention response study.",
        "TRAIT_COVARIATE_SOURCE_ONLY",
    ),
    "10.1186/s12874-015-0091-1":(
        "META_STATISTICAL_METHOD_REFERENCE",
        "RoverKnappFriede2015_meta_method",
        "Meta-analysis variance / small-k statistical method, no ALAN field trial.",
        "METHOD_CITATION_ONLY",
    ),
    "10.1890/11-0423.1":(
        "META_STATISTICAL_METHOD_REFERENCE",
        "Lajeunesse2011_meta_covariance",
        "Meta-analytic correlated effect-size methodology, no ALAN field trial.",
        "METHOD_CITATION_ONLY",
    ),
    "10.3389/fenvs.2017.00061":(
        "ALAN_BIOTIC_MANIPULATED_REPORT",
        "Westhavelland2012_lit_one_site_dark_one_site",
        "Manfrin2017 BACI across same lit and dark riparian sites reused by Manfrin2025, one actual site per lighting condition.",
        "NONREPLICATED_REUSED_PROGRAMME_DO_NOT_DOUBLE_COUNT",
    ),
    # The source occurrence is a publisher URL ending /full, not canonical DOI.
    "10.3389/fpls.2023.1263795/full":(
        "ALAN_ABIOTIC_PLANT_PHYSIOLOGY_FACTORIAL",
        "Wei2023_ALAN_drought_urban_plants",
        "Field block ALAN x drought effects on photosynthesis/plant physiology, not biological interactions; URL suffix /full is not canonical DOI.",
        "COVARIATE_OR_ALTERNATIVE_MECHANISM_NOT_BIOTIC_META",
    ),
}

def triage(rows):
    if len(rows)!=30 or len({r.get("doi_string") for r in rows})!=30:
        raise ValueError("STOP: original P1 DOI source list must have exactly 30 unique strings")
    output=[]
    articles=set()
    for row in rows:
        key=row["doi_string"]
        if row["triage_priority"]!="P1_EXACT_UNMATCHED_WITH_ALAN_IWE_SOURCE":
            raise ValueError("STOP: input has a non-P1 source")
        if key.startswith(DATA_PREFIXES):
            kind="ORIGINAL_DATA_REPOSITORY_IDENTIFIER"
            pg="NOT_AN_INDEPENDENT_PAPER_DOI"
            evidence="Data/archive deposit ID: the linked publication may be separately present in ALAN. Dataset itself cannot be a second trial."
            status="CROSSREFERENCE_TO_SOURCE_PROGRAMME_NOT_A_NEW_STUDY"
            canonical=key
        else:
            if key not in ARTICLES:raise ValueError(f"STOP: unreviewed non-dataset P1 DOI string: {key}")
            kind,pg,evidence,status=ARTICLES[key]
            articles.add(key)
            canonical=key.removesuffix("/full")
        output.append({
            "source_doi_literal":key,
            "canonical_candidate_doi":canonical,
            "source_type_gate":kind,
            "known_or_candidate_physical_programme":pg,
            "design_or_status_basis":evidence,
            "meta_admission_gate":status,
            "source_repositories":row["source_repositories"],
            "source_locations":row["source_file_locations"],
            "independently_adoptable_effect_added":False,
            "source_attribution_note":"Source type based on DOI repository namespace for datasets and named publisher records/repo audits for article DOI; not an automated systematic primary fulltext read",
        })
    if articles!=set(ARTICLES):
        raise ValueError("STOP: planned 10 source-reviewed non-data DOI identities not all present")
    counts=Counter(x["source_type_gate"] for x in output)
    if counts["ORIGINAL_DATA_REPOSITORY_IDENTIFIER"]!=20:
        raise ValueError("STOP: published dataset namespace count changed")
    return output, {
        "snapshot":"2026-11-24",
        "p1_source_doi_strings":30,
        "data_repository_identifiers":20,
        "non_dataset_publication_identifiers":10,
        "type_counts":dict(counts),
        "ALAN_biotic_manipulation_reports":counts["ALAN_BIOTIC_MANIPULATED_REPORT"],
        "Westhavelland_reused_reports":2,
        "ALAN_abiotic_plant_physiology_study":counts["ALAN_ABIOTIC_PLANT_PHYSIOLOGY_FACTORIAL"],
        "new_meta_eligible_programmes_added":0,
        "documented_read_status":"Only ten publisher/repo-source article identities were manually classified; not all 676 DOI strings and not all scientific full texts were read.",
        "study_count_warning":"Even four ALAN-biotic paper reports include one shared Westhavelland infrastructure and previously known Sullivan/Firebaugh experiments. Data repositories and background papers are not independent lighting trials.",
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--p1",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,default=Path("results/alan_crossrepo_p1_source_triage_20261124"))
    args=ap.parse_args()
    with args.p1.open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    out,summary=triage(rows)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/"p1_source_type_reconciled.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(out[0]))
        writer.writeheader();writer.writerows(out)
    (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("P1_SOURCE_TYPE_REVIEW",json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
