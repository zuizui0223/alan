#!/usr/bin/env python3
"""Outcome-blind preflight for historical ALAN at GloPL experimental sites.

Only original SOURCE experimental year_of_study, location, species and
publication identity are included in the signed cohort before any
pollen-limitation effect is inspected for a model. The GloPL 'Year'
column is bibliographic; never substitute it for experimental year.
Published GloPL geo-coordinates may be approximate, not field lux.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

SHA="264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de"
SCHEMA="ALAN_GloPL_true_experiment_year_1993_2015_prePCNL_source_gate_v1"

def require(v,msg):
    if not v:raise ValueError("ALAN_GLOPL_YEAR_PREFLIGHT: "+msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def run(src,out):
    require(sha(src)==SHA,"Original GloPL Dryad-identical file SHA changed")
    raw=pd.read_csv(src,encoding="latin-1",dtype=str,keep_default_na=False)
    require(len(raw)==2969,"Original published experimental case count changed")
    cols=["Latitude","Longitude","year_of_study","Species_accepted_names",
          "DOI","Author","Year","PL_Effect_Size","PL_Effect_Size_Type1",
          "PL_Effect_Size_Type2","Constant_added","Level_of_Supplementation"]
    require(set(cols).issubset(raw.columns),"Original experiment-year or outcome schema changed")
    d=raw[cols].copy()
    d.insert(0,"original_row_id",range(len(d)))
    d["year_exp"]=pd.to_numeric(d.year_of_study,errors="coerce")
    d["lat"]=pd.to_numeric(d.Latitude,errors="coerce")
    d["lon"]=pd.to_numeric(d.Longitude,errors="coerce")
    year=(d.year_exp.between(1993,2015)&d.year_exp.eq(np.floor(d.year_exp)))
    valid=year & d.lat.between(-90,90) & d.lon.between(-180,180)
    # Publication source identity is a cluster, not a locally controlled plot.
    DOI=d.DOI.astype(str).str.strip().str.casefold()
    fallback=(d.Author.astype(str).str.strip()+"|"+
              d.Year.astype(str).str.strip())
    d["study_key"]=np.where(DOI.ne(""),DOI,fallback)
    d=d.loc[valid & d.study_key.ne("|")].copy()
    d["year_exp"]=d.year_exp.astype(int)
    d["species"]=d.Species_accepted_names.str.replace("_"," ",regex=False).str.strip()
    d["years_of_pre_experiment_ALAN"]=d.year_exp-1992
    require((d.years_of_pre_experiment_ALAN>=1).all(),
            "Historical ALAN includes contemporaneous or future experiment year")
    require(not d.original_row_id.duplicated().any(),"Original GloPL row identity duplicated")
    # Keep positive and negative PL effect (log treatment vs natural) without
    # outcome-based filtering. Only finite outcome test occurs in model stage.
    out=Path(out)
    require(not out.exists(),"No overwrite of frozen original GloPL population")
    out.mkdir(parents=True)
    cohort=out/"original_GloPL_experiment_year_georeferenced_1993_2015.csv.gz"
    d.to_csv(cohort,index=False,compression="gzip")
    audit={
      "schema":SCHEMA,"GloPL_source_SHA256":SHA,
      "original_GloPL_rows":len(raw),"selected_original_year_rows":len(d),
      "experimental_year_range":[int(d.year_exp.min()),int(d.year_exp.max())],
      "selected_unique_GloPL_sites":int(d[["lat","lon"]].drop_duplicates().shape[0]),
      "selected_unique_GloPL_site_years":int(d[["lat","lon","year_exp"]].drop_duplicates().shape[0]),
      "selected_unique_published_studies":int(d.study_key.nunique()),
      "selected_unique_accepted_names":int(d.species.nunique()),
      "selected_island_status":"unknown_pending_GSHHG_experiment_point_overlay",
      "GloPL_year_of_study_used_not_bibliographic_publication_year":True,
      "effect_size_not_inspected_for_source_eligibility":True,
      "GloPL_no_moth_specific_pollinator_identity":True,
      "source_frozen_csv_SHA256":sha(cohort),
      "historical_PCNLLight_exposure_not_yet_measured":True
    }
    (out/"glopl_true_experiment_year_prior_light_source_receipt.json").write_text(
        json.dumps(audit,sort_keys=True,indent=2)+"\n")
    print("ALAN_GLOPL_ACTUAL_STUDY_YEAR_PRE_ALAN_GEO_COHORT",
          json.dumps(audit,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
