#!/usr/bin/env python3
"""Collect complete GBIF source-verified 4-family 2013-2025 island site-years.

A source-only observation and sampling coverage result. Zero moth records
in a photographed area is not evidence of true local absence. Plant
pollination, insect density and causal artificial nightlight are UNMEASURED.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import pandas as pd
import numpy as np

FAMILIES=("Noctuidae","Geometridae","Nymphalidae","Pieridae")
MOTH=("Noctuidae","Geometridae")
BUTTERFLY=("Nymphalidae","Pieridae")
SCHEMA="ALAN_GBIF_2026_08_4family_global_GSHHG_island_yearcell_real_counts_v1"
SOURCE_SCHEMA="ALAN_GBIF_202608_true_fourfamily_GSHHG_island_cell_year_occurrences_v1"

def require(v,m):
    if not bool(v):raise ValueError("ALAN_GBIF_TRUE_GLOBAL_COLLECT: "+str(m))

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):h.update(part)
    return h.hexdigest()

def run(root,out):
    root=Path(root);out=Path(out)
    require(not out.exists(),"No overwrite authentic full global observed support")
    frames=[];audits={}
    for family in FAMILIES:
        matches=list(root.rglob(f"{family}_GBIF_202608_real_island_source_audit.json"))
        require(len(matches)==1,
                "Exactly one verified FULL original source per family needed "+family)
        a=json.loads(matches[0].read_text())
        require(a.get("schema")==SOURCE_SCHEMA and
                a.get("source_original_family")==family and
                a.get("source_scan_complete_all_files_all_row_groups") is True and
                a.get("source_scan_max_row_groups_only") is None and
                a.get("source_expected_family_original_all_year_all_basis_records")==
                  a.get("source_actual_scan",{}).get("source_original_records_read") and
                a.get("GBIF_source_2026_08_DOI")=="10.15468/dl.8yrbe7",
                "A source family has incomplete real global Parquet scan "+family)
        c=list(root.rglob(f"{family}_GBIF_202608_original_island_photo_cell_year.csv.gz"))
        require(len(c)==1 and sha(c[0])==a["source_aggregated_cell_CSV_SHA256"],
                "Source family genuine aggregated grid values changed "+family)
        df=pd.read_csv(c[0],dtype={"island":str})
        require(not df.duplicated(["island","year","latbin","lonbin","family"]).any() and
                set(df.family)=={family} and
                df.year.between(2013,2025).all(),
                "Family grid/year/source keys not unique or observed")
        require(len(df)==a["actual_source_island_cell_year_family_rows"] and
                int(df.original_GBIF_human_occurrence_records.sum())==
                  a["source_island_level_photo_observation_record_total"],
                "Source family row count did not reconcile original records")
        frames.append(df)
        audits[family]={"original_source_all_years_all_basis_rows":
                         a["source_expected_family_original_all_year_all_basis_records"],
                        "original_source_human_observation_2013to2025_on_islands":
                         a["source_island_level_photo_observation_record_total"],
                        "actual_island_0p05_yearcells":len(df),
                        "actual_islands_with_family_records":int(df.island.nunique()),
                        "source_aggregated_sha256":a["source_aggregated_cell_CSV_SHA256"]}
    all_records=pd.concat(frames,ignore_index=True)
    require(not all_records.duplicated(["island","year","latbin","lonbin","family"]).any(),
            "Duplicate actual source island/year/cell/family")
    cell=all_records.pivot_table(
        index=["island","year","latbin","lonbin"],
        columns="family",values="original_GBIF_human_occurrence_records",
        fill_value=0,aggfunc="sum").reset_index()
    for fam in FAMILIES:
        if fam not in cell:cell[fam]=0
        cell[fam]=cell[fam].astype("int64")
    cell["moth_family_records"]=cell[list(MOTH)].sum(axis=1)
    cell["butterfly_family_records"]=cell[list(BUTTERFLY)].sum(axis=1)
    cell["all_four_family_records"]=cell.moth_family_records+cell.butterfly_family_records
    cell["relative_moth_family_record_fraction"]=(
       cell.moth_family_records/cell.all_four_family_records.where(
           cell.all_four_family_records.gt(0),np.nan))
    cell["has_at_least_5_four_family_records"]=cell.all_four_family_records.ge(5)
    cell["has_both_moth_and_butterfly_family_records"]=(
       cell.moth_family_records.gt(0)&cell.butterfly_family_records.gt(0))
    eligible=cell.loc[cell.has_at_least_5_four_family_records &
                      cell.has_both_moth_and_butterfly_family_records]
    out.mkdir(parents=True)
    outcsv=out/"GBIF_202608_actual_global_island_0p05deg_year_4family_occurrence_counts.csv.gz"
    cell.to_csv(outcsv,index=False,compression="gzip")
    receipt={
        "schema":SCHEMA,"global_original_snapshot_DOI":"10.15468/dl.8yrbe7",
        "explicit_family_source_count":len(FAMILIES),
        "original_source_family_and_SHA_receipts":audits,
        "source_true_GSHHG_20km2_landmass_universe":"GSHHG 2.3.7 L1, not 4663-island GID alpha universe",
        "original_fourfamily_record_count_on_actual_islands":
           int(cell.all_four_family_records.sum()),
        "observed_original_island_0p05_cell_years":len(cell),
        "observed_actual_distinct_islands":int(cell.island.nunique()),
        "original_source_observation_years":sorted(cell.year.unique().astype(int).tolist()),
        "positive_moth_family_islands":int(cell.loc[
            cell.moth_family_records.gt(0),"island"].nunique()),
        "positive_butterfly_family_islands":int(cell.loc[
            cell.butterfly_family_records.gt(0),"island"].nunique()),
        "same_cell_year_moth_and_butterfly_source_present":int(
            cell.has_both_moth_and_butterfly_family_records.sum()),
        "at_least_5_records_and_both_guilds_cell_year":len(eligible),
        "at_least_5_records_and_both_guilds_islands":int(eligible.island.nunique()),
        "islands_with_ge_two_distinct_supported_source_cells_same_observation_year":int(
            (eligible.groupby(["island","year"]).size().ge(2)).sum()),
        "source_grid_csv_SHA256":sha(outcsv),
        "record_fraction_is_documentation_and_observation_effort_not_capture_abundance":True,
        "source_zeros_in_family_record_category_not_biological_nonoccurrence":True,
        "no_actual_moth_visitation_pollen_deposition_or_population_density_measured":True,
        "GBIF_2026_08_DWCA_is_preview_and_doi_provenance_preserved":True,
        "no_satellite_ALAN_ecological_effect_fitted":True
    }
    (out/"ALAN_2026_08_real_global_fourfamily_island_counts_receipt.json").write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("ALAN_GBIF_TRUE_COMPLETE_GLOBAL_ISLAND_FOURFAMILY_SUPPORT",
          json.dumps({k:receipt[k] for k in (
          "original_fourfamily_record_count_on_actual_islands",
          "observed_original_island_0p05_cell_years",
          "observed_actual_distinct_islands",
          "same_cell_year_moth_and_butterfly_source_present",
          "at_least_5_records_and_both_guilds_cell_year",
          "at_least_5_records_and_both_guilds_islands")},
          sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.root,a.out)
