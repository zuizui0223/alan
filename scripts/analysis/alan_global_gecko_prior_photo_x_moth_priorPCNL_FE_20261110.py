#!/usr/bin/env python3
"""Exploratory ALAN × historically documented gecko on global moth PHOTO records.

Actual Gekkonidae+5 other gecko-family human observations 2013..year-1
from complete 2026-08 frozen GBIF Squamata source (11,634,053 original rows)
moderate true 2015-25 global source moth/family photo proportion. The
reference group is "gecko not documented in available earlier GBIF photos",
NEVER gecko biologically absent, predator absent, not nocturnal or not
feeding on moths. Island×year fixed intercept, genuine island cluster SE,
island-equal weight. Scincidae and documented Hemidactylus/Phelsuma
comparison retain observer/method confounding disclosure. No causal effect.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import (
    fit_source)
from scripts.analysis.alan_global_moth_priorPCNL_darkpatch_conditionalFE_20261110 import (
    SCHEMA as ORIGINAL_DARK_MODEL_SCHEMA)
from scripts.analysis.alan_global_gecko_moth_historical_source_cosupport_20261110 import (
    SCHEMA as GECKO_COSUPPORT_SCHEMA)
from scripts.analysis.alan_global_GBIF_gecko_island_original_stream_20261110 import sha256

SCHEMA="ALAN_2026_11_true_global_GECKO_documented_x_priorPCNL_moth_photo_guild_v1"
RESPONSE="recorded_moth_family_share_pp"

def require(ok,message):
    if not bool(ok):raise ValueError("ALAN_GECKO_MOTH_PREDATOR_SOURCE: "+str(message))

def load_original(joined,receipt):
    report=json.loads(Path(receipt).read_text())
    require(report.get("schema")==GECKO_COSUPPORT_SCHEMA and
            report.get("original_surveyed_moth_insect_geogrid_year_rows")==62688 and
            report.get("genuine_islands_in_source")==196 and
            sha256(joined)==report["original_source_joint_SHA256"],
            "Wrong original 2026-08 source-gecko and genuine historical PCNL join")
    d=pd.read_csv(joined,dtype={"island":str},low_memory=False)
    required={"island","year","latbin","lonbin","island_year",RESPONSE,
              "prior_log1p_PCNL","log_all_fourfamily_observations","centre_lat",
              "area_km2","prior_gecko_documented","prior_scink_documented",
              "prior_Hemidactylus_documented","prior_Phelsuma_documented",
              "prior_any_reptile_photo_documented",
              "source_gecko_same_cell_year","source_scink_same_cell_year"}
    require(required.issubset(d.columns) and len(d)==62688 and
            not d.duplicated(["island","year","latbin","lonbin"]).any() and
            d.year.between(2013,2025).all() and
            d[RESPONSE].between(0,100).all(),
            "Actual original documented-moth and predator source fields incomplete")
    for col in (
        "prior_gecko_documented","prior_scink_documented",
        "prior_Hemidactylus_documented","prior_Phelsuma_documented",
        "prior_any_reptile_photo_documented"):
        d[col]=d[col].astype(str).str.lower().eq("true")
    return d,report

def support_group(d, start=2015,end=2025,effort_min=5,
                  restrict_reptile_photos=False):
    x=d.loc[d.year.between(start,end)&
            d.prior_log1p_PCNL.notna()&
            d.log_all_fourfamily_observations.ge(np.log1p(effort_min))&
            d[RESPONSE].notna()&
            d.centre_lat.notna()&
            d.area_km2.gt(0)].copy()
    if restrict_reptile_photos:
        x=x.loc[x.prior_any_reptile_photo_documented].copy()
    candidates=x.groupby("island_year",observed=True).agg(
        sample_cells=("lonbin","size"),
        n_light_levels=("prior_log1p_PCNL","nunique"),
        span=("prior_log1p_PCNL",lambda s:float(s.max()-s.min())),
        island=("island","first"),
        gecko_prior=("prior_gecko_documented","first"),
        skink_prior=("prior_scink_documented","first"),
        prior_Hemidactylus=("prior_Hemidactylus_documented","first"),
        prior_Phelsuma=("prior_Phelsuma_documented","first"))
    for col in ("prior_gecko_documented","prior_scink_documented",
                "prior_Hemidactylus_documented","prior_Phelsuma_documented"):
        require(x.groupby("island_year")[col].nunique().le(1).all(),
                "Prior-year predator documentation must be constant within an island-year")
    keep=candidates.index[(candidates.sample_cells>=2)&
                          (candidates.n_light_levels>=2)&
                          (candidates.span>=.02)]
    z=x.loc[x.island_year.isin(keep)].copy()
    g=candidates.loc[keep]
    status={"observed_original_insect_photo_cells_selected":len(z),
         "original_islands_with_source_insect_sameyear_gradient":int(z.island.nunique()),
         "actual_island_year_light_gradient_groups":len(g),
         "source_prior_gecko_positive_groups":int(g.gecko_prior.sum()),
         "source_gecko_undocumented_groups":int((~g.gecko_prior).sum()),
         "prior_gecko_documented_true_islands":int(g.loc[g.gecko_prior,"island"].nunique()),
         "prior_gecko_undocumented_true_islands":int(g.loc[~g.gecko_prior,"island"].nunique()),
         "prior_scink_photo_positive_groups":int(g.skink_prior.sum()),
         "source_prior_Hemidactylus_photo_positive_groups":int(g.prior_Hemidactylus.sum()),
         "source_prior_Phelsuma_photo_positive_groups":int(g.prior_Phelsuma.sum()),
         "source_atleast_2_earlier_GBIF_calendar_years":start>=2015,
         "original_reptile_prior_documentation_restricted":restrict_reptile_photos,
         "no_reptile_GBIF_photo_is_not_predator_absence":True}
    return z,status

def fit(d,moderator,years=(2015,2025),only_reptile=False,
        include_lat_slope=False,include_area_slope=False,weight="island_equal",
        effort_min=5):
    x,diag=support_group(d,years[0],years[1],effort_min,only_reptile)
    if not len(x):
        return {"status":"no_eligible_original_photo_source","support":diag}
    observed_prior=x.groupby("island_year")[moderator].first()
    positive_islands=int(x.loc[x[moderator],"island"].nunique())
    comparison_islands=int(x.loc[~x[moderator],"island"].nunique())
    diag["moderator_documented"]=moderator
    diag["moderator_positive_source_islands"]=positive_islands
    diag["moderator_undocumented_source_islands"]=comparison_islands
    diag["source_positive_groups"]=int(observed_prior.sum())
    diag["source_undocumented_groups"]=int((~observed_prior).sum())
    if positive_islands<20 or comparison_islands<20 or x.island.nunique()<40:
        return {"status":"prior_documentation_interaction_unidentifiable_too_few_islands",
                "support":diag}
    if diag["source_positive_groups"]<40 or diag["source_undocumented_groups"]<40:
        return {"status":"too_few_true_island_year_strata",
                "support":diag}
    x["centre_lat_c"]=x.centre_lat-x.centre_lat.median()
    x["gecko_documentation_x_light"]=x.prior_log1p_PCNL*x[moderator].astype(float)
    cols=["prior_log1p_PCNL","gecko_documentation_x_light",
          "log_all_fourfamily_observations","centre_lat_c"]
    if include_lat_slope:
        x["lat_c_x_light"]=x.centre_lat_c*x.prior_log1p_PCNL
        cols.append("lat_c_x_light")
    if include_area_slope:
        x["area_c"]=np.log10(x.area_km2)-np.log10(x.area_km2).median()
        x["logarea_c_x_light"]=x.area_c*x.prior_log1p_PCNL
        cols.append("logarea_c_x_light")
    result=fit_source(x,RESPONSE,weight,cols)
    return {"source_support":diag,
            "source_estimand":"Difference in prior-light slopes between gecko documented vs gecko not documented in earlier GBIF photos, NOT gecko present/absent",
            "source_reference_not_verified_gecko_absent":True,
            "original_island_x_photo_year_fixed_effects":True,
            "island_original_source_cluster_uncertainty":True,
            "observational_detectability_and_invasion_status_unresolved":True,
            "fit":result}

def run(joined,receipt,out):
    d,source=load_original(joined,receipt)
    analyses={
       "gecko_prior_documented_main":fit(d,"prior_gecko_documented"),
       "gecko_prior_documented_with_latitude_light_slope":
              fit(d,"prior_gecko_documented",include_lat_slope=True),
       "gecko_prior_documented_with_island_area_light_slope":
              fit(d,"prior_gecko_documented",include_area_slope=True),
       "gecko_prior_documented_both_geographic_light_slope":
              fit(d,"prior_gecko_documented",include_area_slope=True,include_lat_slope=True),
       "gecko_prior_only_islands_with_any_prior_skink_or_gecko_photo":
              fit(d,"prior_gecko_documented",only_reptile=True),
       "scink_prior_documented_observer_negative_control":
              fit(d,"prior_scink_documented"),
       "Hemidactylus_prior_documented_specific_source":
              fit(d,"prior_Hemidactylus_documented"),
       "Phelsuma_prior_documented_day_gecko_source":
              fit(d,"prior_Phelsuma_documented"),
       "gecko_prior_2020plus":fit(d,"prior_gecko_documented",years=(2020,2025)),
       "gecko_prior_20plus_insect_source_records":
              fit(d,"prior_gecko_documented",effort_min=20),
       "gecko_prior_cell_year_equal_weight":
              fit(d,"prior_gecko_documented",weight="site_equal"),
    }
    report={
       "schema":SCHEMA,
       "GBIF_original_2026_08_frozen_source_DOI":"10.15468/dl.8yrbe7",
       "source_historical_corrected_PCNl_DOI":"10.5281/zenodo.17013414",
       "source_genuine_62688_GBIF_record_fraction_photo_cells":len(d),
       "source_genuine_196_moth_GSHHG_islands":int(d.island.nunique()),
       "source_joined_real_gecko_moth_source_CSV_SHA256":sha256(joined),
       "analyses":analyses,
       "original_gecko_source_absence_not_observed_not_biological_absence":True,
       "gecko_prior_photo_cannot_prove_same_site_moth_predation":True,
       "gecko_taxa_not_all_nocturnal_or_invasive":True,
       "no_true_moth_density_or_effective_pollen_flow_measured":True,
       "analyzed_recorded_photo_source_composition_not_ecological_predation":True,
       "postexploratory_multiple_sensitivity_model_unadjusted_nominal_p":True}
    out=Path(out)
    require(not out.exists(),"Never overwrite authentic gecko source association")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    safe={}
    for name,entry in analyses.items():
        v=entry.get("fit",{})
        safe[name]={"status":v.get("status",entry.get("status")),
                    "support":entry.get("source_support",entry.get("support")),
                    "exposure":v.get("terms",{}).get("prior_log1p_PCNL"),
                    "gecko_interaction":v.get("terms",{}).get("gecko_documentation_x_light")}
    print("ALAN_GENUINE_GECKO_PRIOR_SOURCE_VS_MOTH_GLOBAL_ALAN_MODEL",
          json.dumps(safe,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--joined",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.joined,a.receipt,a.out)
