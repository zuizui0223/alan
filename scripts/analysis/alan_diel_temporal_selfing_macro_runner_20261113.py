#!/usr/bin/env python3
"""Source-verified ALAN global functional compensation opportunity and temporal check.

Result is photographic plant composition and co-documentation, never locally
verified pollinator replacement or light-caused flowering plant decline.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    load_sites, load_labels, original_accepted_species, require
)
from scripts.analysis.alan_real_daynight_reproductive_PCNL_562_20261113 import (
    experimental_daynight
)
from scripts.analysis.alan_diel_39island_leaveout_selfing_20261113 import (
    original_counts, BASELINE, fit_sample
)
from scripts.analysis.alan_direct_selfing_portfolio_mothstudy_ALAN_20261113 import (
    run_selfing_portfolio
)
from scripts.analysis.alan_diel_repeatcell_temporal_20261113 import (
    global_source_and_codocumented_temporal
)


def run(sites_path,globi_path,groups,study,study_traits,island_traits,output):
    if output.exists():raise FileExistsError(str(output))
    sites=load_sites(sites_path)
    require((len(sites),sites.island.nunique())==(77647,562),
            "Original island×source PCNL frame drift")
    labels=load_labels(globi_path)
    flora=original_accepted_species(groups,sites,labels)
    require(len(flora)==2012925 and flora.species.nunique()==14330,
            "Original photographic species source drift")
    studies,origin=experimental_daynight(study,study_traits)
    require(len(studies)==134 and
            int(studies.night_stronger.sum())==75 and
            int(studies.day_stronger.sum())==59,"Original experimental functional classification changed")
    base,species,matched=original_counts(flora,studies,sites)
    baseline=fit_sample(base)
    require(baseline["status"]=="source_photographic_composition_observational" and
            baseline["support"]["evaluated_islands"]==39 and
            baseline["support"]["evaluated_sites"]==2373 and
            abs(baseline["estimate"]["beta_pp"]-BASELINE)<1.e-8,
            "Frozen original 39-island effect not reverified")
    temporal=global_source_and_codocumented_temporal(base,species)
    selfing=run_selfing_portfolio(flora,sites,studies,island_traits)
    result={
        "schema":"ALAN_true562island_independent_selfing_x_night_experiment_rollingpriorPCNL_and_temporal_source_v1",
        "status":"AUTHENTIC_SOURCE_93M_AND_33YEAR_LIGHT_RECOMPUTED_EXPLORATORY",
        "source":{
            "source_photo_ID_index":93209351,
            "original_plant_species_site_years":len(flora),
            "original_accepted_taxa":int(flora.species.nunique()),
            "original_islands":int(sites.island.nunique()),
            "original_source_photo_cell_years":len(sites),
            "original_historical_exposure":"calibrated 1992 to observed-year-1 PCNL mean log1p: NOT actual population trend",
            "independent_night_experiments":"Kendall & Nicholson 2025 131 studies, 228 fruit/seed experiment pairs",
            "independent_autonomous_selfing":"island source-DIRECT high/medium and high-only frozen approved plant trait snapshot",
            "study_species_excluded_from_SELFING_portfolio_numerator_and_denominator":True
        },
        "baseline_previous_analysis":baseline,
        "same_photo_cell_temporal_effect":temporal,
        "independent_direct_selfing_portfolio_and_documentation":selfing,
        "limitations":{
            "genuine_moth_visitation_not_identified_by_diel_exclusion":True,
            "no_actual_pollination_or_seedset_at_photographed_sites":True,
            "external_source_selfing_trait_not_measured_at_photo_site":True,
            "historical_light_rolling_mean_not_year_on_year_lamp_experiment":True,
            "unknown_mating_trait_never_negative_reference":True,
            "night_STUDY_species_record_at_photo_site_not_proven_night_pollinator_service":True,
            "external_night_STUDY_plant_species_excluded_from_selfing_response_for_independence":True,
            "source_species_selection_and_cultivated_vs_native_status_remain":True,
            "exploratory_posthoc_multiple_specs_and_late_hypothesis":True,
            "no_private_species_trait_by_location_exported":True
        }
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False,default=str)+"\n")
    digest={}
    for quality,data in selfing.items():
        interaction=data["external_night_reproduction_species_presence_moderation"]
        display={}
        for key in ("island_equal","site_year_equal"):
            if key in interaction:
                z=interaction[key]
                display[key]={
                    "status":z["status"],"n_islands":z.get("n_islands"),
                    "night_x_ALAN":z.get("coefficients",{}).get(
                        "light_x_night_presence")
                }
        digest[quality]={
            "source_trait_observed":data["source_species_direct_selfing_coverage"],
            "source_known_autonomous":data["observed_autonomous_selfing_ratio_among_directly_known"],
            "night_interaction_support":interaction.get("source_support"),
            "night_light_interaction":display,
            "SI_SC_control_status":data["independent_sc_si_documentation_control"].get("status")
        }
    print("REAL_ALAN_562_SOURCE_TEMPORAL_SELFING",json.dumps({
        "temporal":temporal,"selfing":digest
    },ensure_ascii=False),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for a in ("sites","globi","groups","study","study_traits","island_traits","out"):
        p.add_argument("--"+a,required=True,type=Path)
    a=p.parse_args()
    run(a.sites,a.globi,a.groups,a.study,a.study_traits,a.island_traits,a.out)
