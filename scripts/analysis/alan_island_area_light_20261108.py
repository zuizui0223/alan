#!/usr/bin/env python3
"""ALAN: source-authenticated island area as moderator of observed floral whiteness.

The estimand is the DIFFERENCE IN WITHIN-ISLAND SPATIAL NIGHTLIGHT SLOPE
per tenfold change in GSHHG island landmass area. It is NOT an island-vs-
mainland effect, an isolation-distance effect, or causal light selection.

All 2019/2021/2023/2024/2025 original-photo ROI9 inputs are reused with their
existing SHA gates. No new image downloads, no imputed ROI colours, no 2020
Global1 and no incomplete 2022 VIIRS source. Photos/taxa are the sampled
assemblage, not botanical-confirmed native vegetation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import statsmodels.formula.api as smf
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as panel
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as tier
from scripts.analysis.mois_global_2020_GSHHG_scope_20261107 import prepare_land
from scripts.analysis.mois_global2025_original_source_20261107 import COAST_SHA

EXPECTED = {"strict_5p3r3s":7982, "measured_3p2r2s":20543}
MIN_SUPPORT_ISLANDS = 10
YEARS = (2019,2021,2023,2024,2025)

def source_land_area(coast):
    coast=Path(coast)
    rec=json.loads((coast/"source_receipt.json").read_text())
    panel.require(rec.get("original_archive_sha256")==COAST_SHA,
                  "GSHHG original archive source drift")
    for rel,digest in rec["extracted_sha256"].items():
        panel.require(panel.sha(coast/rel)==digest, "GSHHG extract drift: "+rel)
    land=prepare_land(gpd.read_file(coast/"GSHHS_h_L1.shp"))
    area=land[["land_id","land_class","area_km2"]].copy()
    area.land_id=area.land_id.astype(str)
    panel.require(area.land_id.is_unique and area.area_km2.gt(0).all(),
                  "Source landmass identity/area invalid")
    return area

def join_area_and_select(d, source_area):
    """Use real distinct 0.05-degree sites, never multiple years as new sites."""
    x=d.loc[d.is_island.eq(1)].copy()
    ids=x.target_id.astype(str).str.split("|",regex=False,expand=True)
    panel.require(ids.shape[1]==3 and ids[0].eq(x.year.astype(str)).all(),
                  "Original cell-year target IDs changed")
    x["land_id"]=x.land_id.astype(str)
    panel.require(ids[2].eq(x.land_id).all(), "Original landmass ID mismatch")
    x["site_id"]=ids[1].to_numpy()
    lookup=source_area.set_index("land_id")
    x["area_km2"]=x.land_id.map(lookup.area_km2)
    x["land_area_class"]=x.land_id.map(lookup.land_class)
    panel.require(x.area_km2.notna().all() and x.area_km2.gt(0).all() and
                  x.land_area_class.eq("island").all(),
                  "Unsupported original GSHHG island/area join")
    groups=x.groupby("land_id",observed=True).agg(
       n_sites=("site_id","nunique"), n_years=("year","nunique"),
       n_cell_years=("target_id","size"), n_distinct_light=("log_light","nunique"))
    spatial=x.groupby(["land_id","site_id"],observed=True).log_light.median().groupby(level=0).agg(
       lambda a:float(a.max()-a.min()))
    groups["within_island_spatial_loglight_span"]=spatial
    supported=groups.index[(groups.n_sites>=2)&(groups.n_distinct_light>=2)]
    within=x.loc[x.land_id.isin(supported)].copy()
    panel.require(not within.duplicated("target_id").any(),
                  "Repeated island observation used twice")
    return x,within,groups

def fit_area_moderator(within, centre_log10_area, equal_islands=False):
    count=int(within.land_id.nunique())
    panel.require(count>=MIN_SUPPORT_ISLANDS,
                  "Insufficient distinct islands with 2+ true sites and 2+ light levels")
    a=within.copy()
    a["log10_area_c"]=np.log10(a.area_km2)-centre_log10_area
    panel.require(a.log10_area_c.groupby(a.land_id).nunique().eq(1).all(),
                  "Island area must be invariant for all source cells and years")
    panel.require(a.log10_area_c.std()>0 and a.log_light.std()>0,
                  "Area or nightlight contrast absent")
    formula=("white_pp ~ log_light + log_light:log10_area_c + log_pop + "
             "bio5 + log_rain + elevation + log_attempts + C(land_id) + C(year)")
    # Two distinct pre-specified observational targets:
    # cell-equal: every actually measured photo cell-year contributes equally;
    # island-equal: each GSHHG island has total weight one, regardless of its
    # volunteer/photo coverage, preventing the largest photographed islands
    # from accounting for most model weight.
    weights=1.0/a.groupby("land_id",observed=True).land_id.transform("size").astype(float)
    panel.require(np.allclose(weights.groupby(a.land_id).sum().to_numpy(),1.0),
                  "Equal-island weights do not sum to one per landmass")
    if equal_islands:
        fit=smf.wls(formula,data=a,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":a.land_id},use_t=True)
    else:
        fit=smf.ols(formula,data=a).fit(
            cov_type="cluster",cov_kwds={"groups":a.land_id},use_t=True)
    terms={}
    for term in ("log_light","log_light:log10_area_c"):
        panel.require(term in fit.params and np.isfinite(fit.bse[term]),
                      "Area moderation unidentified or unstable")
        lo,hi=fit.conf_int().loc[term]
        terms[term]={
            "beta":float(fit.params[term]),"SE":float(fit.bse[term]),
            "ci_low":float(lo),"ci_high":float(hi),"p":float(fit.pvalues[term])}
    return {"model":formula,
            "estimand":"nightlight-whiteness slope change per tenfold island AREA (km^2); within-island site/years",
            "centre_log10_area_km2":float(centre_log10_area),
            "cluster_unit":"GSHHG landmass (island)",
            "weighting":"equal_total_per_island" if equal_islands else "equal_per_observed_cell_year",
            "distinct_islands":count,
            "cell_years":len(a),"terms":terms,
            "nominal_exploratory_CIs":True}

def run(paths,coast,out):
    out=Path(out)
    panel.require(not out.exists(),"Do not overwrite signed ALAN ecological result")
    frames=tier.build_frozen_five_years(paths)
    cohorts=tier.by_year_tier(frames)
    land=source_land_area(coast)
    current={}
    for name,n in EXPECTED.items():
        d=cohorts[name]
        panel.require(len(d)==n and set(d.year)==set(YEARS) and d.target_id.is_unique,
                      "Original frozen source cohort changed")
        all_island,within,groups=join_area_and_select(d,land)
        current[name]=(all_island,within,groups)
    strict_supported=current["strict_5p3r3s"][1]
    strict_unit_areas=strict_supported.groupby("land_id").area_km2.first()
    centre=float(np.log10(strict_unit_areas).median())
    main={}
    for name,(all_island,within,groups) in current.items():
        unit=within.groupby("land_id").area_km2.first()
        q1,q2=np.quantile(np.log10(unit),[1/3,2/3])
        strata=pd.cut(np.log10(unit),bins=[-np.inf,q1,q2,np.inf],
                      labels=["smaller_third","middle_third","larger_third"],
                      include_lowest=True)
        island_groups=strata.rename("area_bin").to_frame()
        island_groups["n_cells_years"]=within.groupby("land_id").size()
        totals=island_groups.groupby("area_bin",observed=True).agg(
            n_islands=("n_cells_years","size"),n_cell_years=("n_cells_years","sum"))
        model=fit_area_moderator(within,centre,equal_islands=False)
        equal_island_model=fit_area_moderator(within,centre,equal_islands=True)
        main[name]={
           "all_island_photo_cell_years":len(all_island),
           "all_island_landmasses":int(all_island.land_id.nunique()),
           "eligible_spatially_repeated_islands":int(len(unit)),
           "selected_original_photo_cell_years":len(within),
           "islands_with_spatial_loglight_span_ge_0p05":int(
               groups.loc[unit.index,"within_island_spatial_loglight_span"].ge(0.05).sum()),
           "island_area_km2_min":float(unit.min()),
           "island_area_km2_median_equal_island":float(unit.median()),
           "island_area_km2_max":float(unit.max()),
           "observed_island_area_thirds":[
               {"area_bin":str(k),
                "independent_islands":int(v.n_islands),
                "photographic_cell_years":int(v.n_cell_years)}
               for k,v in totals.iterrows()],
           "fit":model,
           "equal_island_fit":equal_island_model
        }
    panel.require(main["strict_5p3r3s"]["all_island_photo_cell_years"]==710 and
                  main["strict_5p3r3s"]["all_island_landmasses"]==101 and
                  main["strict_5p3r3s"]["eligible_spatially_repeated_islands"]==44,
                  "Previously signed 710-cell/101-landmass/44-island support changed")
    report={
       "schema":"ALAN_source_real_island_area_nightlight_moderation_v1",
       "years":list(YEARS),
       "strict_photo_cell_years":7982,
       "measured_broad_photo_cell_years":20543,
       "area_source":"signed GSHHG 2.3.7 L1 EPSG:6933 polygon area_km2",
       "island_definition":"original source land_class island, landmass area <= 7 million km^2",
       "fixed_centre_equal_island_strict_median_log10_area":centre,
       "results":main,
       "not_isolation_or_pollinator_or_causal_effect":True,
       "source_photo_white_ROI_not_true_flower_species_percentage":True,
       "never_imputed_unknown_colour":True,
       "new_images_or_ROI_inferences":0,
       "source_SHA256":{**{k:panel.sha(v) for k,v in paths.items()},
                       "coast_original_archive_sha256":COAST_SHA},
    }
    out.mkdir(parents=True)
    (out/"island_area_light_source_results.json").write_text(
       json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_ACTUAL_GSHHG_ISLAND_AREA_LIGHT_RESULTS",
          json.dumps({name:{
            "n_islands":data["eligible_spatially_repeated_islands"],
            "n_cell_years":data["selected_original_photo_cell_years"],
            "area_moderator":data["fit"]["terms"]["log_light:log10_area_c"],
            "nightlight_at_strict_median_area":data["fit"]["terms"]["log_light"],
            "island_equal_area_moderator":data["equal_island_fit"]["terms"]["log_light:log10_area_c"],
            "island_equal_light_at_strict_median_area":data["equal_island_fit"]["terms"]["log_light"]}
            for name,data in main.items()},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    source_keys=("canonical","wcvp","geo2019","viirs2019","pop2019",
                 "geo2021","viirs2021","pop2021","geo2023","viirs2023","pop2023",
                 "geo2024","cn2024","uc2024","uc2024receipt",
                 "geo2025","geo2025receipt","viirs2025")
    for key in source_keys:
        parser.add_argument("--"+key.replace("_","-"),dest=key,type=Path,required=True)
    parser.add_argument("--coast",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    run({k:getattr(args,k) for k in source_keys},args.coast,args.out)
