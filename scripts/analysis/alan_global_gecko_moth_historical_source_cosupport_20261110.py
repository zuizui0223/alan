#!/usr/bin/env python3
"""Join actually documented prior-year gecko & skink records to global moth photos.

Source: original GBIF 2026-08 full Squamata observations on GSHHG islands,
separate from original 2013–25 photographed moth/butterfly four-family source
and year-minus-one PCNL (62,688 primary eligible photo-grid site-years;
196 independent source islands). A gecko NOT photographed is "undocumented",
NOT absent; a gecko photographed does NOT prove eating a moth.

A PRIOR gecko record (2013..year-1) is a temporally ordered archival
source-evidence moderator, not the prior-year biological occupancy truth.
Require 2015+ for >=2 source-gecko documentation years and quantify actual
reptile-photo effort co-support, Scincidae observation comparator, site-year
and year-level overlap. No causal ALAN × predator conclusions in this stage.
"""
from __future__ import annotations
import argparse,collections,hashlib,json
from pathlib import Path
import pandas as pd
import numpy as np

from scripts.analysis.alan_global_trueGBIF_PCNl_13year_darkpatch_20261110 import source_eligible
from scripts.analysis.alan_global_GBIF_gecko_island_original_stream_20261110 import (
    SCHEMA as REPTILE_SCHEMA,GECKO_FAMILIES,sha256 as file_sha)

SCHEMA="ALAN_2026_11_genuine_global_gecko_moth_GBIF_prior_documentation_join_v1"
BASELINE_SHA="64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852"
START=2015

def require(ok,msg):
    if not bool(ok):raise ValueError("ALAN_GECKO_GBIF_JOIN: "+str(msg))

def prepare(gecko,genus,receipt):
    r=json.loads(Path(receipt).read_text())
    require(r["schema"]==REPTILE_SCHEMA and
            r["full_original_Squamata_all_year_and_basis_rows"]==11634053 and
            file_sha(gecko)==r["source_gecko_family_cell_year_CSV_SHA256"] and
            file_sha(genus)==r["source_gecko_genus_cell_year_CSV_SHA256"],
            "Gecko source is not the complete original 2026-08 Squamata file")
    g=pd.read_csv(gecko,dtype={"island":str})
    h=pd.read_csv(genus,dtype={"island":str})
    keys=["island","year","latbin","lonbin"]
    require(not g.duplicated(keys+["actual_source_taxon_family"]).any() and
            g.year.between(2013,2025).all() and
            set(g.actual_source_taxon_family).issubset(set(GECKO_FAMILIES)|{"Scincidae"}),
            "Actual original gecko and Scincidae source unexpected")
    require(not h.duplicated(keys+["original_genus"]).any() and
            h.year.between(2013,2025).all(),
            "Named gecko genus original time/location not clean")
    g["source_gecko_record_count"]=np.where(g.actual_source_taxon_family.isin(
            GECKO_FAMILIES),g.actual_GBIF_human_observation_records,0)
    g["source_scink_record_count"]=np.where(g.actual_source_taxon_family.eq(
            "Scincidae"),g.actual_GBIF_human_observation_records,0)
    cell=g.groupby(keys,sort=False)[[
        "source_gecko_record_count","source_scink_record_count"]].sum().reset_index()
    year=g.groupby(["island","year"],sort=False)[[
        "source_gecko_record_count","source_scink_record_count"]].sum().reset_index()
    h["source_Hemidactylus"]=np.where(h.original_genus.eq("Hemidactylus"),
                          h.original_source_records,0)
    h["source_Phelsuma"]=np.where(h.original_genus.eq("Phelsuma"),
                          h.original_source_records,0)
    genus_year=h.groupby(["island","year"],sort=False)[[
        "source_Hemidactylus","source_Phelsuma"]].sum().reset_index()
    year=year.merge(genus_year,on=["island","year"],how="left",validate="one_to_one")
    year[["source_Hemidactylus","source_Phelsuma"]]=year[[
        "source_Hemidactylus","source_Phelsuma"]].fillna(0)
    require((year.source_gecko_record_count+year.source_scink_record_count>0).all(),
            "Empty source gecko/scink island-year in original")
    require(int(year.source_gecko_record_count.sum())==
            r["source_gecko_photo_observation_records"] and
            int(year.source_scink_record_count.sum())==
            r["source_Scincidae_photo_observation_records"],
            "Reptile source photo effort not reconciled")
    return cell,year,r

def link(base,cell,year):
    k=["island","year","latbin","lonbin"]
    require(len(base)==62688 and base.island.nunique()==196,
            "Not original genuine global 62,688 moth source eligible sites")
    out=base.copy()
    same=cell.rename(columns={
        "source_gecko_record_count":"source_gecko_same_cell_year",
        "source_scink_record_count":"source_scink_same_cell_year"})
    out=out.merge(same,on=k,how="left",validate="one_to_one")
    out[["source_gecko_same_cell_year","source_scink_same_cell_year"]]=out[[
        "source_gecko_same_cell_year","source_scink_same_cell_year"]].fillna(0)
    at_year=year.rename(columns={
        "source_gecko_record_count":"source_gecko_same_island_year",
        "source_scink_record_count":"source_scink_same_island_year",
        "source_Hemidactylus":"source_Hemidactylus_same_island_year",
        "source_Phelsuma":"source_Phelsuma_same_island_year"})
    out=out.merge(at_year,on=["island","year"],how="left",validate="many_to_one")
    for name in at_year.columns:
        if name not in ("island","year"):out[name]=out[name].fillna(0)
    priors=[]
    for observation_year in range(START,2026):
        prior=year.loc[year.year.lt(observation_year)].groupby("island")[[
            "source_gecko_record_count","source_scink_record_count",
            "source_Hemidactylus","source_Phelsuma"]].sum()
        grp=out.loc[out.year.eq(observation_year),[
            "island","year","latbin","lonbin"]].copy()
        if len(grp)==0:continue
        for field in ("source_gecko_record_count","source_scink_record_count",
                      "source_Hemidactylus","source_Phelsuma"):
            grp["prior_"+field]=grp.island.map(prior[field]).fillna(0)
        priors.append(grp)
    prev=pd.concat(priors,ignore_index=True)
    out=out.merge(prev,on=k,how="left",validate="one_to_one")
    out["prior_gecko_documented"]=out.prior_source_gecko_record_count.gt(0)
    out["prior_scink_documented"]=out.prior_source_scink_record_count.gt(0)
    out["prior_any_reptile_photo_documented"]=(out.prior_gecko_documented |
                                             out.prior_scink_documented)
    out["prior_Hemidactylus_documented"]=out.prior_source_Hemidactylus.gt(0)
    out["prior_Phelsuma_documented"]=out.prior_source_Phelsuma.gt(0)
    out["original_gecko_window_available_ge2years"]=out.year.ge(START)
    # The source gapped years are 'no source records' not absence.
    require(len(out)==len(base) and
            not out.duplicated(k).any() and
            (out.prior_source_gecko_record_count.loc[
                out.year.eq(2013)].isna()).all(),
            "Original observation year or source-gecko prior-year mismatch")
    return out

def report_coverage(d):
    w=d.loc[d.original_gecko_window_available_ge2years].copy()
    grouped=w.groupby(["island","year"],sort=False).agg(
         photo_cells=("lonbin","size"),
         true_light_levels=("prior_log1p_PCNL","nunique"),
         min_PCNL=("prior_log1p_PCNL","min"),
         max_PCNL=("prior_log1p_PCNL","max"),
         prior_gecko=("prior_gecko_documented","first"),
         prior_scink=("prior_scink_documented","first"),
         prior_Hemidactylus=("prior_Hemidactylus_documented","first"),
         prior_Phelsuma=("prior_Phelsuma_documented","first")).reset_index()
    groups=grouped.loc[grouped.photo_cells.ge(2)&grouped.true_light_levels.ge(2)&(
                     grouped.max_PCNL-grouped.min_PCNL>=.02)]
    diagnostic={"eligible_2015_2025_moth_site_years":len(w),
        "eligible_original_independent_moth_islands":int(w.island.nunique()),
        "year_2015to2025_true_within_island_light_groups":len(groups),
        "source_prior_gecko_documented_islands":int(
            w.loc[w.prior_gecko_documented,"island"].nunique()),
        "source_prior_scink_documented_islands":int(
            w.loc[w.prior_scink_documented,"island"].nunique()),
        "source_prior_any_reptile_documented_islands":int(
            w.loc[w.prior_any_reptile_photo_documented,"island"].nunique()),
        "true_prior_gecko_positive_island_year_groups_with_light":
           int(groups.prior_gecko.sum()),
        "true_prior_gecko_not_documented_island_year_groups_with_light":
           int((~groups.prior_gecko).sum()),
        "true_prior_skink_positive_island_year_groups_with_light":
           int(groups.prior_scink.sum()),
        "true_prior_Hemidactylus_photo_positive_groups":
           int(groups.prior_Hemidactylus.sum()),
        "true_prior_Phelsuma_photo_positive_groups":
           int(groups.prior_Phelsuma.sum()),
        "source_moth_sites_with_gecko_same_0p05_photo_cell_year":
           int(w.source_gecko_same_cell_year.gt(0).sum()),
        "source_moth_sites_with_skink_same_0p05_photo_cell_year":
           int(w.source_scink_same_cell_year.gt(0).sum()),
        "source_moth_islands_with_any_same_year_gecko_source":
           int(w.loc[w.source_gecko_same_island_year.gt(0),"island"].nunique()),
        "source_moth_islands_with_any_earlier_gecko_record":
           int(w.loc[w.prior_gecko_documented,"island"].nunique()),
        "conditional_true_reptile_photo_background_islands":int(
           w.loc[w.prior_any_reptile_photo_documented,"island"].nunique()),
        "NO_prior_gecko_record_is_NOT_gecko_free_island":True,
        "same_cell_gecko_photo_not_direct_observed_predation":True,
        "no_2026_or_future_gecko_detection_leaked_into_older_insects":True}
    return diagnostic

def run(source,source_receipt,gecko,genus,gecko_receipt,out):
    require(file_sha(source)=="64b269dd796b9095853c371909a0dd760921df2e396636a83fc59a2d659e4852",
            "Historical globally original moth light photo-cell source differs")
    base=source_eligible(source,source_receipt)
    cell,year,gecko_info=prepare(gecko,genus,gecko_receipt)
    joint=link(base,cell,year)
    stats=report_coverage(joint)
    out=Path(out);require(not out.exists(),"Never overwrite actual gecko moth source")
    out.mkdir(parents=True)
    csv=out/"ALAN_real_original_62688_moth_source_gecko_prior_photo_documentation.csv.gz"
    joint.to_csv(csv,index=False,compression="gzip")
    receipt={
        "schema":SCHEMA,
        "original_GBIF_2026_08_source_DOI":"10.15468/dl.8yrbe7",
        "original_surveyed_moth_insect_geogrid_year_rows":62688,
        "genuine_islands_in_source":196,
        "real_source_gecko_and_scink_coverage":stats,
        "original_source_joint_SHA256":file_sha(csv),
        "original_source_gecko_family_SHA256":gecko_info[
            "source_gecko_family_cell_year_CSV_SHA256"],
        "source_classification_gecko_documented_not_population_occupancy":True,
        "source_earlier_gecko_photos_not_proof_present_at_moth_visit_date":True,
        "source_Hemidactylus_photo_not_automatically_invasive":True,
        "source_Phelsuma_day_geckos_can_shift_night_under_ALAN":True,
        "source_skink_records_are_observer_activity_proxy_not_true_predator_absence":True,
        "no_moth_pollen_or_gecko_predation_interaction_fitted":True}
    (out/"ALAN_2026_11_original_global_gecko_moth_photo_source_join_receipt.json").write_text(
         json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("ALAN_ACTUAL_GLOBAL_GECKO_MOTH_PREYEAR_SOURCE_COSUPPORT",
          json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for field in ("source","source_receipt","gecko","genus","gecko_receipt","out"):
        p.add_argument("--"+field.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.source_receipt,a.gecko,a.genus,a.gecko_receipt,a.out)
