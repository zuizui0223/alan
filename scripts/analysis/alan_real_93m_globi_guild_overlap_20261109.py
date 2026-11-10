#!/usr/bin/env python3
"""Reconstruct moth/butterfly *documentation* overlap from real 93m photo source.

Sixteen original GitHub Actions group SQLite files (run 42616279235) are
biological source artifacts, not disjoint accepted-species lists. Unify each
island×photo-cell×year×WCVP accepted species before counting any group;
duplicate species recorded in multiple original hash groups count ONCE.
Join to the SHA-pinned 85,318 eligible site-years independently selected
before ALAN or colour. GloBI no record never means no actual pollination.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA="e7800476ae6d35acb9ef71525f32521596f884733e2e8ee9f87fda7acb8e9eab"
CATEGORIES=("moth","butterfly","both","moth_record_no_butterfly",
            "butterfly_record_no_moth","neither_documented",
            "moth_pollinates","night_pollen")
CREATE="""CREATE TABLE unit (
 island TEXT,latbin INTEGER,lonbin INTEGER,year INTEGER,species TEXT,
 moth INTEGER,butterfly INTEGER,pollinates INTEGER,night_pollen INTEGER,
 PRIMARY KEY(island,latbin,lonbin,year,species))"""

def need(ok,msg):
    if not ok: raise ValueError("ALAN_GUILDS_93M: "+msg)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def combine(groups,source,out):
    out=Path(out)
    need(not out.exists(),"Refuse replacing previous empirical evidence")
    need(sha(source)==SOURCE_SHA,"Precolour source coverage changed")
    sources=[Path(groups)/f"allflora-island-group-{g:02d}.sqlite" for g in range(16)]
    need(all(x.is_file() for x in sources),"Full 16 group signed real GSHHG source required")
    out.mkdir(parents=True)
    conn=sqlite3.connect(out/"temp_join.sqlite")
    conn.execute(CREATE)
    for g,path in enumerate(sources):
        conn.execute("ATTACH DATABASE ? AS original",(str(path),))
        conn.execute("""INSERT INTO unit
          SELECT island,latbin,lonbin,year,species,
            (moth_positive_photo_ids>0),(butterfly_positive_photo_ids>0),
            (moth_pollinates_photo_ids>0),(night_pollen_photo_ids>0)
          FROM original.flora WHERE 1
          ON CONFLICT(island,latbin,lonbin,year,species)
          DO UPDATE SET
             moth=MAX(unit.moth,excluded.moth),
             butterfly=MAX(unit.butterfly,excluded.butterfly),
             pollinates=MAX(unit.pollinates,excluded.pollinates),
             night_pollen=MAX(unit.night_pollen,excluded.night_pollen)""")
        conn.commit()
        conn.execute("DETACH DATABASE original")
    q="""SELECT island,latbin,lonbin,year,COUNT(*) AS all_species,
        SUM(moth) AS moth,SUM(butterfly) AS butterfly,
        SUM(moth*butterfly) AS both,
        SUM(moth*(1-butterfly)) AS moth_record_no_butterfly,
        SUM((1-moth)*butterfly) AS butterfly_record_no_moth,
        SUM((1-moth)*(1-butterfly)) AS neither_documented,
        SUM(pollinates) AS moth_pollinates,
        SUM(night_pollen) AS night_pollen
        FROM unit GROUP BY island,latbin,lonbin,year"""
    cat=pd.read_sql_query(q,conn)
    conn.close()
    site=pd.read_csv(source)
    cat["island"]=cat.island.astype(str)
    site["island"]=site.island.astype(str)
    d=site.merge(cat,on=["island","latbin","lonbin","year"],
                 how="left",validate="one_to_one",indicator=True)
    need(len(d)==85318 and d._merge.eq("both").all(),
         "Source island-year group reconciliation is incomplete")
    need(d.all_species.eq(d.all_accepted_recorded_species).all()
         and d.moth.eq(d.moth_visitor_positive_recorded_species).all(),
         "Independent original documented guild counts disagree with signed source")
    need((d.both + d.moth_record_no_butterfly +
          d.butterfly_record_no_moth + d.neither_documented == d.all_species).all(),
         "Four evidence categories do not close to source denominator")
    need(d.moth_pollinates.le(d.moth).all() and d.night_pollen.le(d.moth).all(),
         "Nocturnal pollen / GloBI pollinates must have a positive source moth claim")
    result={"schema":"ALAN_full93m_original_85318_GloBI_moth_butterfly_accepted_species_source_overlap_v1",
            "GSHHG_source_photo_run":42616279235,
            "signed_site_source_SHA256":SOURCE_SHA,
            "n_real_site_years":len(d),
            "n_real_islands":int(d.island.nunique()),
            "source_species_site_units":int(d.all_species.sum()),
            "evidence_categories":{}}
    for key in CATEGORIES:
        ratio=d[key]/d.all_species
        result["evidence_categories"][key]={
            "species_site_units":int(d[key].sum()),
            "sites_with_at_least_one_recorded_species":int(d[key].gt(0).sum()),
            "site_equal_mean_fraction":float(ratio.mean()),
            "island_equal_mean_fraction":float(ratio.groupby(d.island).mean().mean())}
    result["moth_positive_also_butterfly_positive_fraction"]=(
        result["evidence_categories"]["both"]["species_site_units"]/
        result["evidence_categories"]["moth"]["species_site_units"])
    result["moth_record_without_recorded_butterfly_is_NOT_proven_moth_exclusive"]=True
    result["both_recorded_does_NOT_prove_two_effective_pollination_services"]=True
    result["guild_no_record_is_NOT_biological_absence"]=True
    result["all_source_photo_UUID_not_botanical_native_population"]=True
    (out/"source_guild_overlap_result.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    d.drop(columns="_merge").to_csv(out/"source_guild_overlap_site_year.csv.gz",
                                      index=False,compression="gzip")
    (out/"temp_join.sqlite").unlink(missing_ok=True)
    print("ALAN_93M_GLOBI_SOURCE_EVIDENCE_OVERLAP",json.dumps(result,sort_keys=True))
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--groups",type=Path,required=True)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();combine(a.groups,a.source,a.out)
if __name__=="__main__": main()
