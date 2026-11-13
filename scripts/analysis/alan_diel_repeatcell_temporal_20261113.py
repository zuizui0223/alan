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
    median_span=metadata["median_within_cell_historic_PCNL_mean_span"]
    if median_span is None or median_span<=0:
        return {"status":"NOT_ESTIMABLE_NO_PRIOR_LIGHT_CHANGE","support":metadata}
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
    return {
        "full_562_island_original_photo_site_year_response":full,
        "full_562_island_source_site_equal_sensitivity":full_site_equal,
        "original_39island_common_source_response":narrower,
        "frozen_parent_comparison":{
            "original_39island_source_same_year_spatial_beta":BASELINE,
            "important":"Different site and year support and FE terms; do not compare numeric beta sizes as formal mediation"
        }
    }
