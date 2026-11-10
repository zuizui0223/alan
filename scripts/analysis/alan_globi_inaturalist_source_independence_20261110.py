#!/usr/bin/env python3
"""Reclassify real photo-recorded accepted plants after removing iNaturalist-sourced GloBI claims.

The 61,967 GloBI moth/butterfly source claims include a published
iNaturalist-derived interaction archive. This is a SOURCE-INDEPENDENCE
sensitivity (not local nocturnal pollination verification).
Rebuild species×island×0.05°photo cell×year from the 16 genuine disjoint
93m original photo-source SQLite groups; collapse duplicate source group
photo records to ONE accepted taxon per cell-year before comparing to PCNL.
Never relabel 'not in GloBI' as no moth visitation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import alan_full93m_moth_PCNLPrior_two_way_FE_20261109 as model

GLOBI_SHA="a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994"
PCNL_SHA="440a8313d71168e0a4dfd2f41a7bef45623abb068b07e3282026f1c35dccbd2c"
SOURCE_GROUPS=16
SOURCE_ROWS=85318
FIELDS=("moth_external","butterfly_external","any_external",
        "moth_inat_only","butterfly_inat_only","both_external",
        "external_moth_no_external_butterfly")


def require(ok,message):
    if not ok: raise ValueError("ALAN_GLOBI_INAT_SOURCE: "+message)


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()


def source_labels(globi):
    require(sha256(globi)==GLOBI_SHA,"GloBI source study snapshot changed")
    g=pd.read_csv(globi,compression="gzip",dtype=str).fillna("")
    require(len(g)==61967 and {"plant_name","guild","relation","study_title"}<=set(g),
            "Not the pinned 61,967 GloBI study-claim cache")
    is_inat=g.study_title.str.contains("iNaturalist",case=False,regex=False)
    is_moth=g.guild.eq("moth_family_visit")
    is_butterfly=g.guild.eq("butterfly_family_visit")
    require(set(g.guild)=={"moth_family_visit","butterfly_family_visit"},
            "Unexpected GloBI guild; all-guild information not available")
    g["moth_global"]=is_moth
    g["butterfly_global"]=is_butterfly
    g["moth_external"]=is_moth & ~is_inat
    g["butterfly_external"]=is_butterfly & ~is_inat
    g["moth_inat"]=is_moth & is_inat
    g["butterfly_inat"]=is_butterfly & is_inat
    f=g.groupby("plant_name",observed=True)[
        ["moth_global","butterfly_global","moth_external","butterfly_external",
         "moth_inat","butterfly_inat"]].any().reset_index()
    f=f.rename(columns={"plant_name":"species"})
    f["any_external"]=f.moth_external | f.butterfly_external
    f["both_external"]=f.moth_external & f.butterfly_external
    f["moth_inat_only"]=f.moth_global & ~f.moth_external
    f["butterfly_inat_only"]=f.butterfly_global & ~f.butterfly_external
    f["external_moth_no_external_butterfly"]=f.moth_external & ~f.butterfly_external
    source_summary={
       "GloBI_total_study_claims":len(g),
       "GloBI_iNaturalist_title_claims":int(is_inat.sum()),
       "GloBI_moth_plant_names":int(f.moth_global.sum()),
       "GloBI_moth_only_iNaturalist_plant_names":int(f.moth_inat_only.sum()),
       "GloBI_moth_non_iNaturalist_plant_names":int(f.moth_external.sum()),
       "GloBI_butterfly_only_iNaturalist_plant_names":int(f.butterfly_inat_only.sum())}
    return f,source_summary


def original_species(groups,sites,study):
    keys=["island","latbin","lonbin","year"]
    source_keys=sites[keys].copy()
    source_keys["island"]=source_keys.island.astype(str)
    require(source_keys.duplicated(keys).sum()==0,"Photo-site×year source duplicated")
    blocks=[]
    for group in range(SOURCE_GROUPS):
        path=Path(groups)/f"allflora-island-group-{group:02d}.sqlite"
        require(path.is_file(),f"Missing original disjoint source group {group}")
        with sqlite3.connect(path) as conn:
            d=pd.read_sql_query("""SELECT island,latbin,lonbin,year,species,
                 moth_positive_photo_ids,butterfly_positive_photo_ids
                 FROM flora""",conn)
        d["island"]=d.island.astype(str)
        d=d.merge(source_keys,on=keys,how="inner",validate="many_to_one")
        d["original_moth"]=d.moth_positive_photo_ids.gt(0)
        d["original_butterfly"]=d.butterfly_positive_photo_ids.gt(0)
        blocks.append(d[keys+["species","original_moth","original_butterfly"]])
    d=pd.concat(blocks,ignore_index=True)
    key=keys+["species"]
    repeated=d.loc[d.duplicated(key,keep=False)]
    if len(repeated):
        conflicts=repeated.groupby(key,observed=True)[
            ["original_moth","original_butterfly"]].nunique()
        require(int(conflicts.to_numpy().max())==1,
                "Disjoint photo original groups disagree on source guild identity")
    d=d.drop_duplicates(key)
    d=d.merge(study,on="species",how="left",validate="many_to_one")
    guild=[c for c in study.columns if c!="species"]
    d[guild]=d[guild].fillna(False).astype(bool)
    require(d.original_moth.eq(d.moth_global).all() and
            d.original_butterfly.eq(d.butterfly_global).all(),
            "Reclassified GloBI accepted names differ from original GSHHG source labels")
    return d


def run(pcnl,globi,groups,out):
    require(sha256(pcnl)==PCNL_SHA,"Original 33-year PCNL source changed")
    site=pd.read_csv(pcnl)
    require(len(site)==SOURCE_ROWS,"Genuine same-source accepted plant denominator changed")
    site["island"]=site.island.astype(str)
    study,source_receipt=source_labels(globi)
    species=original_species(groups,site,study)
    keys=["island","latbin","lonbin","year"]
    counts=species.groupby(keys,observed=True).agg(
        accepted_species=("species","size"),
        moth_global=("moth_global","sum"),
        butterfly_global=("butterfly_global","sum"),
        **{k:(k,"sum") for k in FIELDS}).reset_index()
    joined=site.merge(counts,on=keys,how="left",validate="one_to_one",indicator=True)
    require(joined._merge.eq("both").all(),
            "Original same-source actual photographed-site rows missing")
    require(joined.accepted_species.eq(joined.all_accepted_recorded_species).all(),
            "Accepted species per original photo site did not conserve denominator")
    require(joined.moth_global.eq(joined.moth_visitor_positive_recorded_species).all(),
            "Known original moth-positive source numerator did not reproduce")
    joined=joined.drop(columns="_merge")
    current=model.prepare(joined)
    require(len(current)==77647 and current.island.nunique()==562,
            "Original within-island light-support cohort unexpectedly changed")
    results={}
    for k in ["moth_global","butterfly_global",*FIELDS]:
        q=current.copy()
        q["frac_pp"]=100.0*q[k].astype(float)/q.all_accepted_recorded_species
        variants={
           "island_equal":model.fit(q,"island_equal"),
           "island_equal_species_adjusted":model.fit(q,"island_equal",True),
        }
        require(all(x["status"]=="estimated_exploratory_observational"
                    for x in variants.values()),"Could not identify original source beta")
        results[k]={"n_species_site_records":int(current[k].sum()),
                    "models":{name:x["terms"][model.LIGHT] for name,x in variants.items()}}
    get=lambda k:results[k]["models"]["island_equal"]["beta_pp"]
    require(abs(get("moth_global")-.44933020729608203)<1.e-6,
            "Original empirical moth positive beta was not reproduced")
    require(abs(get("moth_global")-get("moth_external")-
                get("moth_inat_only"))<1.e-8,
            "Published source independence components do not reconstruct regression")
    report={
       "schema":"ALAN_real93m_GloBI_iNaturalist_source_independence_v1",
       "source_GloBI":source_receipt,
       "original_93m_photo_IDs":93209351,
       "original_source_site_years":SOURCE_ROWS,
       "original_within_island_site_years":len(current),
       "original_within_island_islands":int(current.island.nunique()),
       "unique_accepted_species_site_years":len(species),
       "estimates":results,
       "iNaturalist_origin_removed_from_labels_not_from_plant_photo_observations":True,
       "non_iNaturalist_GloBI_source_not_necessarily_independent_of_other_citizen_science":True,
       "non_iNaturalist_GloBI_positive_not_proven_effective_moth_pollination":True,
       "GloBI_nonrecord_not_biological_absence":True,
       "native_or_uncultivated_plant_status_not_in_original_93m_source":True}
    path=Path(out)
    require(not path.exists(),"Refuse output overwrite")
    path.mkdir(parents=True)
    (path/"ALAN_original_GloBI_iNaturalist_source_removed_20261110.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("ALAN_GLOBI_NON_INAT_SOURCE_REAL_562ISLANDS",
          json.dumps({k:x["models"]["island_equal"] for k,x in results.items()},
                     sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pcnl",required=True,type=Path)
    p.add_argument("--globi",required=True,type=Path)
    p.add_argument("--groups",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.pcnl,a.globi,a.groups,a.out)


if __name__=="__main__": main()
