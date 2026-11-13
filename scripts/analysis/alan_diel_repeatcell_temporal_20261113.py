#!/usr/bin/env python3
"""Original same-photo-cell repeated-year negative control for day/night reproduction × PCNL.

Question: Does the previously observed cross-sectional documented-reproductive
night-share versus historical nightlight slope persist within an unchanged real
island geographic cell across years? This does NOT measure real flower
population trend, actual moth pollination, or historical lighting intervention.

Use original fixed 562-island 77,647 source photo-cell × observation-year frame,
source-verified experiment comparisons and preceding-year calibrated PCNL.
No source locations, individual species identities or trait-by-site rows emitted.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from scripts.analysis.alan_real_genus_prevalence_bias_20261110 import (
    KEY, LIGHT, fixed_effects_regression, require,
)
from scripts.analysis.alan_diel_39island_leaveout_selfing_20261113 import (
    supported, BASELINE,
)

def repeat_panel(sites: pd.DataFrame, min_repeats=2):
    """Repeated geographic PHOTO source cell, not a revisited plant individual."""
    x = sites.copy()
    x["source_cell_id"] = (x.island.astype(str)+"|"+
                           x.latbin.astype(str)+"|"+x.lonbin.astype(str))
    n_years=x.groupby("source_cell_id",observed=True).year.nunique()
    keep=set(n_years.loc[n_years.ge(min_repeats)].index)
    x=x[x.source_cell_id.isin(keep)].copy()
    islandyear=x.island.astype(str).str.cat(x.year.astype(str),sep="|")
    light_by_site=x.groupby("source_cell_id",observed=True)[LIGHT].agg(["min","max"])
    site_span=(light_by_site["max"]-light_by_site["min"])
    return x,{
        "original_source_photo_site_years":len(sites),
        "repeat_source_photo_cells":len(keep),
        "repeated_photo_cell_years":len(x),
        "repeated_islands":int(x.island.nunique()),
        "repeated_island_by_year_groups":int(islandyear.nunique()),
        "repeated_cells_with_nonzero_prior_mean_PCNL_change":int(site_span.gt(0).sum()),
        "independent_islands_with_nonzero_within_cell_PCNL_change":int(
            x.loc[x.source_cell_id.isin(set(site_span.loc[site_span.gt(1e-9)].index)),
                  "island"].nunique()),
        "median_within_cell_historic_PCNL_mean_span":float(site_span.median()) if len(site_span) else None,
        "max_within_cell_historic_PCNL_mean_span":float(site_span.max()) if len(site_span) else None,
        "repeated_cell_location_not_original_population_identity":True
    }


def longitudinal_regression(sites:pd.DataFrame, response:np.ndarray, name:str,
                            *, weight:str="island_equal",
                            min_islands:int=30, min_rows:int=250):
    """Source panel with true photo CELL FE plus island×year FE.

    Exposure is rolling mean of PREVIOUS annual radiance values, so apparent
    within-cell change is not a direct illumination intervention or lamp lux.
    """
    x=sites.copy()
    require(len(x)==len(response),"Response must align original source photo site rows")
    x["outcome"]=np.asarray(response,float)
    repeated,metadata=repeat_panel(x)
    metadata["response"]=name
    metadata["weight"]=weight
    if len(repeated)<min_rows or repeated.island.nunique()<min_islands:
        return {"status":"NOT_ESTIMABLE_REPEATED_SOURCE_CELL_SUPPORT","support":metadata}
    if (metadata["repeated_cells_with_nonzero_prior_mean_PCNL_change"]<100 or
        metadata["independent_islands_with_nonzero_within_cell_PCNL_change"]<min_islands):
        return {"status":"NOT_ESTIMABLE_INSUFFICIENT_WITHIN_CELL_PRELIGHT_VARIATION",
                "support":metadata}
    # A zero median exposure change does not imply zero identifying support:
    # use the actual repeated-site fixed effects rank check below.
    site=repeated.source_cell_id.to_numpy()
    iy=repeated.island.astype(str).str.cat(repeated.year.astype(str),sep="|").to_numpy()
    if weight=="island_equal":
        counts=repeated.groupby("island",observed=True).island.transform("size")
        w=1.0/counts.to_numpy(float)
    elif weight=="site_equal":
        w=np.ones(len(repeated))
    else:raise ValueError("unknown source weight")
    try:
        m=fixed_effects_regression(
            repeated.outcome.to_numpy(float),
            repeated[LIGHT].to_numpy(float),
            repeated.log_obs.to_numpy(float),w,
            repeated.island.to_numpy(),[site,iy])
        return {
            "status":"EXPLORATORY_WITHIN_SOURCE_PHOTO_CELL_TEMPORAL_ASSOCIATION",
            "support":metadata,
            "model":m,
            "source_photo_observer_species_turnover_unmeasured":True,
            "no_local_plant_population_followed":True,
            "past_PCNL_rolling_mean_not_actual_LAMP_TREATMENT":True
        }
    except (ValueError,np.linalg.LinAlgError) as e:
        return {"status":"NOT_ESTIMABLE_SPARSE_FIXED_EFFECT_OR_RANK",
                "support":metadata,"reason":str(e)[:180]}



def changing_cell_temporal_regression(sites:pd.DataFrame,
                                      response:np.ndarray,
                                      name:str,
                                      *, minimum_within_cell_span=.01,
                                      weight:str="island_equal"):
    """Posthoc SOURCE EXPOSURE variation gate, not a treatment choice.

    Only repeated original photographed cells whose rolling preceding-years'
    log1p PCNL mean varies by at least .01. This is not a global randomized
    sample and zero-change cells are retained in the separate full cohort
    source-support counts. Do not claim illumination intervention or decline.
    """
    require(len(sites)==len(response),"Misaligned repeated-cell outcome")
    frame=sites.copy()
    frame["outcome"]=np.asarray(response,float)
    repeated,base_support=repeat_panel(frame)
    span=repeated.groupby("source_cell_id",observed=True)[LIGHT].agg(
        lambda z:float(z.max()-z.min()))
    allowed=set(span[span.ge(minimum_within_cell_span)].index)
    changed=repeated.loc[repeated.source_cell_id.isin(allowed)].copy()
    support={
        **base_support,
        "posthoc_exposure_span_threshold":minimum_within_cell_span,
        "changed_cell_count":len(allowed),
        "changed_cell_site_years":len(changed),
        "changed_cell_independent_islands":int(changed.island.nunique()),
        "changed_cell_observed_years":int(changed.year.nunique()),
        "photo_cell_years_selection_is_exposure_only":True,
        "not_original_full_global_floristic_sample":True,
    }
    if len(changed)<300 or changed.island.nunique()<30 or changed.year.nunique()<2:
        return {"status":"NOT_ESTIMABLE_INDEPENDENT_CHANGED_PHOTO_CELL_SUPPORT",
                "source_support":support}
    y=changed.outcome.to_numpy(float)
    w=(1/changed.groupby("island",observed=True).island.transform("size").to_numpy(float)
       if weight=="island_equal" else np.ones(len(changed)))
    iy=changed.island.astype(str).str.cat(changed.year.astype(str),sep="|").to_numpy()
    try:
        ans=fixed_effects_regression(
           y,changed[LIGHT].to_numpy(float),
           changed.log_obs.to_numpy(float),w,
           changed.island.to_numpy(),
           [changed.source_cell_id.to_numpy(),iy])
        return {"status":"EXPLORATORY_ORIGINAL_CHANGED_PHOTO_CELL_WITHIN_SOURCE_TEMPORAL",
                "fit":ans,
                "response":name,
                "weight":weight,
                "source_support":support,
                "L_variation_is_rolling_prior_mean_not_lamp_onset":True,
                "no_measured_local_reproductive_or_pollinator_function":True}
    except (ValueError,np.linalg.LinAlgError) as e:
        return {"status":"NOT_ESTIMABLE_ORIGINAL_TEMPORAL_DESIGN_OR_SOLVER",
                "reason":str(e)[:180],"source_support":support}


def global_source_and_codocumented_temporal(s:pd.DataFrame,study_source:pd.DataFrame):
    """Compare full source nighttime-vs-daytime documented plant fraction with
    39-island common-source bin contrast under repeatedly sampled cells."""
    require(not study_source.duplicated(KEY+["species"]).any(),
            "Original species source duplicated")
    g=s.copy()
    require(g[["night","day"]].notna().all().all(),"Functional count missing")

    # Full-source night-minus-day real documented plant fractions:
    # many genuinely present species have no experimental source; unknown,
    # NEVER presumed biologically no-night/day pollination.
    y=(g.night-g.day).to_numpy(float)/g.all_accepted_recorded_species.to_numpy(float)
    full=longitudinal_regression(g,y,"externally_documented_night_minus_day_all_photographed_taxa")
    full_site_equal=longitudinal_regression(g,y,
        "externally_documented_night_minus_day_all_photographed_taxa",
        weight="site_equal")
    q,info=supported(g)
    require(len(q)==2373 and q.island.nunique()==39,
            "Original 39-island common-source sites drift")
    qy=q.night.to_numpy(float)/(q.night+q.day).to_numpy(float)
    narrower=longitudinal_regression(
        q,qy,"night_share_among_both_diel_experiment_source_documented")
    changed=changing_cell_temporal_regression(
        g,y,"externally_documented_night_minus_day")
    changed_eq=changing_cell_temporal_regression(
        g,y,"externally_documented_night_minus_day",weight="site_equal")
    documented_y=(g.night+g.day).to_numpy(float)/g.all_accepted_recorded_species.to_numpy(float)
    changed_study=changing_cell_temporal_regression(
        g,documented_y,"control_total_photo_reproductive_experiment_species_documented")
    return {
        "full_562_island_original_photo_site_year_response":full,
        "changed_cell_exposure_supported_within_cell":changed,
        "changed_cell_site_equal_sensitivity":changed_eq,
        "control_experiment_study_documentation_on_identical_exposure_qualified_cells":changed_study,
        "full_562_island_source_site_equal_sensitivity":full_site_equal,
        "original_39island_common_source_response":narrower,
        "frozen_parent_comparison":{
            "original_39island_source_same_year_spatial_beta":BASELINE,
            "important":"Different site and year support and FE terms; do not compare numeric beta sizes as formal mediation"
        }
    }
