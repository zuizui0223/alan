#!/usr/bin/env python3
"""Actual same-cell/year gecko-skink-moth GBIF documentation overlap gate.

Reuses the already SHA-verified original GBIF 2026-08 Squamata and 4-family
moth/butterfly sources, recorded at GSHHG 196 island / 62,688 true 2013–25
insect-bearing 0.05 degree cell-years, with year-before-original PCNL.

No unrecorded gecko, skink or moth is imputed biologically absent.
Photo of gecko and photo of moth in a ~5km cell in the same CALENDAR YEAR
does NOT establish a predation event, shared nocturnal hour, or lamp-level
co-occurrence. A negative reptile record is UNKNOWN occupancy. The purpose
is to gate whether source-co-documented sites can support a gecko-vs-skink
light-use hypothesis without the 7-island undocumented reptile comparator.
"""
from __future__ import annotations
import argparse, json, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA="ALAN_2026_11_gecko_moth_skink_true_colocated_photo_record_support_v1"
SOURCE_SCHEMA="ALAN_2026_11_genuine_global_gecko_moth_GBIF_prior_documentation_join_v1"
ORIGINAL_ROWS=62688
DOI="10.15468/dl.8yrbe7"

def require(condition,reason):
    if not bool(condition):raise ValueError("ALAN_TRUE_GECKO_COOBSERVED: "+str(reason))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def support(source,receipt):
    meta=json.loads(Path(receipt).read_text())
    require(meta["schema"]==SOURCE_SCHEMA and
            meta["original_surveyed_moth_insect_geogrid_year_rows"]==ORIGINAL_ROWS and
            meta["genuine_islands_in_source"]==196 and
            sha(source)==meta["original_source_joint_SHA256"],
            "Not exactly original source-verified gecko/moth evidence")
    d=pd.read_csv(source,dtype={"island":str},low_memory=False)
    keys=["island","year","latbin","lonbin"]
    require(len(d)==ORIGINAL_ROWS and not d.duplicated(keys).any(),
            "Not one original geographically co-observed moth source-cell/year")
    required={"source_gecko_same_cell_year","source_scink_same_cell_year",
              "source_gecko_same_island_year","source_scink_same_island_year",
              "prior_log1p_PCNL","recorded_moth_family_share_pp",
              "moth_family_records","butterfly_family_records"}
    require(required.issubset(d.columns),"Original gecko source columns absent")
    for name in ("source_gecko_same_cell_year","source_scink_same_cell_year",
                 "source_gecko_same_island_year","source_scink_same_island_year"):
        require(d[name].notna().all() and d[name].ge(0).all() and
                (d[name]%1).eq(0).all(),
                "Original photographed reptile records non-integer or missing")
    require(d.prior_log1p_PCNL.notna().all() and d.year.between(2013,2025).all(),
            "Original year-before ALAN exposure source missing")
    d["gecko_doc"]=d.source_gecko_same_cell_year.gt(0)
    d["skink_doc"]=d.source_scink_same_cell_year.gt(0)
    d["both_lizard_guilds_doc"]=d.gecko_doc&d.skink_doc
    d["one_lizard_guild_doc"]=d.gecko_doc|d.skink_doc
    d["at_least_two_reptile_photo_events"]=(
        d.source_gecko_same_cell_year+d.source_scink_same_cell_year).ge(2)
    d["at_least_five_reptile_photo_events"]=(
        d.source_gecko_same_cell_year+d.source_scink_same_cell_year).ge(5)

    gecko=d.loc[d.gecko_doc].copy()
    skink=d.loc[d.skink_doc].copy()
    together=d.loc[d.both_lizard_guilds_doc].copy()
    either=d.loc[d.one_lizard_guild_doc].copy()
    group=either.groupby(["island","year"],sort=False).agg(
        original_real_cell_years=("lonbin","size"),
        gecko_doc_cells=("gecko_doc","sum"),
        skink_doc_cells=("skink_doc","sum"),
        both_doc_cells=("both_lizard_guilds_doc","sum"),
        pcNl_min=("prior_log1p_PCNL","min"),
        pcNl_max=("prior_log1p_PCNL","max"),
        distinct_PCNL_levels=("prior_log1p_PCNL","nunique"))
    support_both=group.loc[group.gecko_doc_cells.gt(0)&group.skink_doc_cells.gt(0)]
    gradients=support_both.loc[
       support_both.original_real_cell_years.ge(2) &
       support_both.distinct_PCNL_levels.ge(2) &
       (support_both.pcNl_max-support_both.pcNl_min).ge(.02)]
    inside=group.loc[group.gecko_doc_cells.ge(2)&
                     group.skink_doc_cells.ge(2)&
                     group.original_real_cell_years.ge(4)&
                     group.distinct_PCNL_levels.ge(2)&
                     (group.pcNl_max-group.pcNl_min).ge(.02)]
    true_both=together.groupby(["island","year"]).size()
    per_year=[]
    for y,g in d.groupby("year",sort=True):
        per_year.append({"observation_year":int(y),"true_moth_plantphoto_cells":len(g),
            "source_documented_gecko_same_cell_years":int(g.gecko_doc.sum()),
            "source_documented_skink_same_cell_years":int(g.skink_doc.sum()),
            "same_source_gecko_and_skink_cell_years":int(g.both_lizard_guilds_doc.sum()),
            "source_islands_with_photographed_geckos":int(g.loc[g.gecko_doc,"island"].nunique())})
    report={
       "schema":SCHEMA,"original_gbif_frozen_reptile_and_insect_DOI":DOI,
       "original_source_join_sha256":meta["original_source_joint_SHA256"],
       "total_authenticated_global_moth_source_cell_years":len(d),
       "total_source_independent_original_islands":int(d.island.nunique()),
       "same_year_gecko_photographed_moth_source_cell_years":len(gecko),
       "same_year_gecko_photographed_moth_source_islands":int(gecko.island.nunique()),
       "same_year_skink_photographed_moth_source_cell_years":len(skink),
       "same_year_skink_photographed_moth_source_islands":int(skink.island.nunique()),
       "same_cell_year_Gecko_and_Scincidae_photographed_together":len(together),
       "same_cell_year_both_lizard_guilds_original_islands":int(together.island.nunique()),
       "same_island_year_both_lizard_guilds_documented_somewhere":len(support_both),
       "same_island_year_both_lizard_guilds_documented_independent_original_islands":
           int(support_both.reset_index().island.nunique()),
       "source_eligible_island_year_groups_with_both_lizard_guilds_and_actual_light_gradient":
           len(gradients),
       "source_eligible_independent_islands_with_both_lizard_guilds_and_light_gradient":
           int(gradients.reset_index().island.nunique()),
       "stronger_same_island_year_2plus_gecko_photo_cells_and_2plus_skink_photo_cells_with_light":
           len(inside),
       "stronger_same_island_year_eligible_original_independent_islands":
           int(inside.reset_index().island.nunique()),
       "same_photo_cell_year_gecko_and_skink_documented_and_at_least_2_reptile_photo_events":
           int(together.at_least_two_reptile_photo_events.sum()),
       "same_photo_cell_year_gecko_and_skink_documented_and_at_least_5_reptile_photo_events":
           int(together.at_least_five_reptile_photo_events.sum()),
       "source_by_year":per_year,
       "original_0p05deg_cell_nearby_not_moth_predation":True,
       "no_gecko_photograph_not_biological_absence":True,
       "source_lizard_observations_not_standardized_light_census":True,
       "source_moth_GBIF_fraction_not_true_insect_population_abundance":True,
       "real_source_yearly_ALAN_pixel_not_local_streetlight_or_lux":True,
       "no_observed_direct_gecko_predation_moth_diet_or_effective_pollination":True,
       "no_gecko_light_selectivity_coefficient_fitted_yet":True}
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    out=a.out
    require(not out.exists(),"Do not overwrite original source eligibility")
    result=support(a.source,a.receipt)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ALAN_TRUE_GECKO_SKINK_MOTH_COLLOCATED_SOURCE_SUPPORT",
          json.dumps({k:v for k,v in result.items() if k not in
                     ("source_by_year","original_source_join_sha256")},
                     sort_keys=True),flush=True)
if __name__=="__main__":main()
