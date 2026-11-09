#!/usr/bin/env python3
"""Count BOTH WCVP accepted flora and GloBI moth-positive flora by real island sites.

Single pass over 93,209,351 SHA-pinned ORIGINAL public photo metadata IDs.
No RGB downloads or inferred floral white colour; GloBI visit-positive is NOT
effective nocturnal pollination, and no GloBI record means UNKNOWN (not absence).
The SAME observation opportunity defines numerator and denominator; photo IDs
and true observation UUIDs are counted independently. Exactly one original
observation UUID contributes to observational effort (despite sibling photos).
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import tarfile

import numpy as np
import pandas as pd

from scripts.analysis import alan_moth_species_historical_pcnl_inventory_20261109 as original
from scripts.analysis import alan_full93m_globi_moth_original_photo_census_20261109 as source

YEARS = (2019, 2021, 2023, 2024, 2025)
FIELDS = source.COLUMNS
TABLE = """CREATE TABLE IF NOT EXISTS flora (
 island TEXT NOT NULL, latbin INTEGER NOT NULL, lonbin INTEGER NOT NULL,
 year INTEGER NOT NULL, species TEXT NOT NULL,
 photo_ids INTEGER NOT NULL, observation_uuids INTEGER NOT NULL,
 moth_positive_photo_ids INTEGER NOT NULL, moth_positive_observation_uuids INTEGER NOT NULL,
 moth_pollinates_photo_ids INTEGER NOT NULL, night_pollen_photo_ids INTEGER NOT NULL,
 butterfly_positive_photo_ids INTEGER NOT NULL,
 lon REAL NOT NULL, lat REAL NOT NULL,
 PRIMARY KEY (island,latbin,lonbin,year,species)
)"""
UPSERT = """INSERT INTO flora VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
ON CONFLICT(island,latbin,lonbin,year,species) DO UPDATE SET
 photo_ids=photo_ids+excluded.photo_ids,
 observation_uuids=observation_uuids+excluded.observation_uuids,
 moth_positive_photo_ids=moth_positive_photo_ids+excluded.moth_positive_photo_ids,
 moth_positive_observation_uuids=moth_positive_observation_uuids+excluded.moth_positive_observation_uuids,
 moth_pollinates_photo_ids=moth_pollinates_photo_ids+excluded.moth_pollinates_photo_ids,
 night_pollen_photo_ids=night_pollen_photo_ids+excluded.night_pollen_photo_ids,
 butterfly_positive_photo_ids=butterfly_positive_photo_ids+excluded.butterfly_positive_photo_ids"""


def require(ok, reason):
    original.require(ok, "full93m full-flora denominator: " + str(reason))


def normalize_rows(data, accepted_map, evidence):
    require(tuple(data.columns) == FIELDS, "Original complete 93m row format changed")
    require(data.photo_id.str.fullmatch(r"[0-9]+").all()
            and data.observation_uuid.ne("").all()
            and data.inat_taxon_id.ne("").all(), "Original photo/UUID/taxon identity malformed")
    x = data.loc[data.valid_observed_date.eq("true") &
                 data.valid_geopoint.eq("true")].copy()
    x["year"] = pd.to_numeric(x.observed_on.str.slice(0,4),errors="coerce")
    x = x.loc[x.year.isin(YEARS)].copy()
    x["species"] = x.inat_taxon_id.map(accepted_map)
    x = x.loc[x.species.notna()].copy()
    x["latitude"] = pd.to_numeric(x.latitude, errors="coerce")
    x["longitude"] = pd.to_numeric(x.longitude, errors="coerce")
    x = x.loc[x.latitude.between(-90,90) & x.longitude.between(-180,180)].copy()
    x["year"] = x.year.astype(int)
    x["moth_mask"] = x.inat_taxon_id.map(evidence).fillna(0).astype(int)
    return x


def grouped_island_rows(x, islands, seen_uuids):
    if x.empty:
        return [],0,0
    within, ambiguous = original.join_geography(x, islands)
    if within.empty:
        return [],0,ambiguous
    within["latbin"] = np.floor((within.latitude+90)*20).astype(int)
    within["lonbin"] = np.floor((within.longitude+180)*20).astype(int)
    new = ~within.observation_uuid.isin(seen_uuids)
    new &= ~within.observation_uuid.duplicated(keep="first")
    seen_uuids.update(within.observation_uuid.tolist())
    within["obs"] = new.astype(int)
    within["moth_photo"] = ((within.moth_mask & 1)!=0).astype(int)
    within["moth_obs"] = within.moth_photo * within.obs
    within["pollinates"] = ((within.moth_mask & 2)!=0).astype(int)
    within["night_pollen"] = ((within.moth_mask & 4)!=0).astype(int)
    within["butterfly"] = ((within.moth_mask & 8)!=0).astype(int)
    g = within.groupby(["island","latbin","lonbin","year","species"],
                       sort=False,observed=True).agg(
        photo_ids=("photo_id","size"),
        observation_uuids=("obs","sum"),
        moth_positive_photo_ids=("moth_photo","sum"),
        moth_positive_observation_uuids=("moth_obs","sum"),
        moth_pollinates_photo_ids=("pollinates","sum"),
        night_pollen_photo_ids=("night_pollen","sum"),
        butterfly_positive_photo_ids=("butterfly","sum"),
        lon=("longitude","first"),lat=("latitude","first")).reset_index()
    return list(g.itertuples(index=False,name=None)),len(within),ambiguous


def count_group(args):
    require(args.group in range(16), "Original source group must be 0..15")
    require(not args.db.exists(), "Refuse to overwrite existing source group result")
    receipt=json.loads(args.receipt.read_text())
    registry=source.source_registry(args.registry,args.group)
    require(receipt.get("schema")=="mois_onepass_full93m_disjoint_photo_target_tiers_v1"
           and receipt.get("source_snapshot_SHA256")==registry["full_source_manifest_sha256"]
           and len(receipt.get("parts",[]))==8
           and set(int(x["part"]) for x in receipt["parts"]) ==
                 set(range(8*args.group,8*(args.group+1))),
           "Original 93m group provenance receipt differs")
    accepted=original.accepted_taxa(args.wcvp,args.names)
    evidence,_=source.evidence_map(args.wcvp,args.names,args.globi)
    islands=original.load_islands(args.coast)
    db=sqlite3.connect(args.db)
    db.execute(TABLE)
    seen=set()
    counts=Counter()
    file_shas={}
    try:
        for tier in range(14):
            filename=f"tier-{tier:02d}-group-{args.group:02d}.tar"
            path=args.archives/filename
            digest=original.sha256(path)
            require(digest==registry["published_tier_archive_SHA256"][filename],
                    "Original complete 93m source tier SHA changed")
            file_shas[filename]=digest
            tier_n=0
            with tarfile.open(path,mode="r:") as tar:
                items={v.name:v for v in tar.getmembers() if v.isfile()}
                expected={f"tier-{tier:02d}/part-{p:03d}.csv.gz"
                    for p in range(args.group*8,(args.group+1)*8)}
                require(set(items)==expected, "Original full-tier archive parts incomplete")
                for item in sorted(expected):
                    with gzip.GzipFile(fileobj=tar.extractfile(items[item])) as stream:
                        reader=pd.read_csv(stream,dtype=str,keep_default_na=False,
                                           chunksize=args.chunksize)
                        for raw in reader:
                            tier_n+=len(raw)
                            counts["all_93209351_source_photos_examined"]+=len(raw)
                            counts["all_source_photos_with_GloBI_moth_visitor_positive"]+=sum(
                                bool(evidence.get(k,0)&1) for k in raw.inat_taxon_id)
                            x=normalize_rows(raw,accepted,evidence)
                            counts["all_fiveyear_accepted_source_photos_valid_geo"]+=len(x)
                            rows,n,ambiguous=grouped_island_rows(x,islands,seen)
                            counts["all_fiveyear_accepted_source_photos_exact_GSHHG_island"]+=n
                            counts["ambiguous_island_polygon_photo_ids"]+=ambiguous
                            db.executemany(UPSERT,rows)
            require(tier_n==int(receipt["counts_by_disjoint_tier"][str(tier)]),
                    "Original tier counts did not conserve source")
            counts["completed_verified_original_tiers"]+=1
            db.commit()
        require(counts["all_93209351_source_photos_examined"] ==
                int(receipt["original_photo_IDs_in_this_group"]),
                "Original full93m group photo counts mismatch")
        counts["distinct_accepted_island_observation_UUIDs"]=len(seen)
    finally:
        db.close()
    args.db.with_suffix(".json").write_text(json.dumps({
        "schema":"ALAN_full93m_genuine_WCVP_allflora_moth_evidence_island_source_group_v1",
        "group":args.group,"counts":dict(counts),
        "source_tier_sha256":file_shas,
        "registry_sha256":source.REGISTRY_SHA,
        "wcvp_sha256":source.NAMES_SHA,
        "accepted_sha256":source.ACCEPTED_SHA,
        "GloBI_sha256":source.GLOBI_SHA,
        "source_group_receipt_sha256":original.sha256(args.receipt),
        "no_GloBI_record_means_UNKNOWN_not_no_moth":True,
        "photo_not_plant_abundance":True,
        "source_RGB_or_ROI9_measured":0},sort_keys=True,indent=2)+"\n")


def export_query(db,query,path):
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        cur=db.execute(query)
        writer=csv.writer(f)
        writer.writerow([col[0] for col in cur.description])
        for item in cur:
            writer.writerow(item)


def collect(args):
    require(not args.out.exists(), "Do not overwrite pooled 93m source census")
    files=[args.groups/f"allflora-island-group-{i:02d}.sqlite" for i in range(16)]
    require(all(p.is_file() and p.with_suffix(".json").is_file() for p in files),
            "Missing one of 16 original signed full-flora source groups")
    args.out.mkdir(parents=True)
    dest=sqlite3.connect(args.out/"allflora_moth_island.sqlite")
    dest.execute(TABLE)
    totals=Counter()
    for i,path in enumerate(files):
        receipt=json.loads(path.with_suffix(".json").read_text())
        require(receipt.get("schema") ==
                "ALAN_full93m_genuine_WCVP_allflora_moth_evidence_island_source_group_v1"
                and receipt.get("group")==i
                and receipt.get("counts",{}).get("completed_verified_original_tiers")==14,
                "Incomplete full-flora tier source group")
        totals.update(receipt["counts"])
        with sqlite3.connect(path) as db:
            cur=db.execute("SELECT * FROM flora")
            while True:
                batch=cur.fetchmany(5000)
                if not batch:break
                dest.executemany(UPSERT,batch)
            dest.commit()
    require(totals["all_93209351_source_photos_examined"]==source.PHOTO_TOTAL,
            "No claim of full 93,209,351 photographs without full scan")
    require(totals["all_source_photos_with_GloBI_moth_visitor_positive"]==34_346_428,
            "GloBI positive numerator not equal original independently verified source")
    export_query(dest, """SELECT island,latbin,lonbin,year,
        COUNT(DISTINCT species) AS all_accepted_recorded_species,
        SUM(CASE WHEN moth_positive_photo_ids>0 THEN 1 ELSE 0 END) AS
            moth_visitor_positive_recorded_species,
        SUM(photo_ids) AS all_accepted_photo_ids,
        SUM(moth_positive_photo_ids) AS moth_positive_photo_ids,
        SUM(observation_uuids) AS all_accepted_observation_uuids,
        SUM(moth_positive_observation_uuids) AS moth_positive_observation_uuids,
        CAST(SUM(CASE WHEN moth_positive_photo_ids>0 THEN 1 ELSE 0 END) AS REAL)
            / COUNT(DISTINCT species) AS moth_positive_fraction_OF_RECORDED_WCVP_SPECIES,
        MIN(lon) AS lon,MIN(lat) AS lat
        FROM flora GROUP BY island,latbin,lonbin,year
        ORDER BY island,year,latbin,lonbin""",
        args.out/"allflora_moth_island_site_year.csv")
    export_query(dest, """SELECT island,year,
        COUNT(DISTINCT species) AS all_accepted_recorded_species,
        COUNT(DISTINCT CASE WHEN moth_positive_photo_ids>0 THEN species END) AS
             moth_positive_recorded_species,
        SUM(photo_ids) AS all_accepted_photo_ids,
        SUM(moth_positive_photo_ids) AS moth_positive_photo_ids,
        SUM(observation_uuids) AS all_accepted_observation_uuids,
        SUM(moth_positive_observation_uuids) AS moth_positive_observation_uuids,
        CAST(COUNT(DISTINCT CASE WHEN moth_positive_photo_ids>0 THEN species END)
             AS REAL)/COUNT(DISTINCT species) AS moth_positive_fraction_OF_RECORDED_WCVP_SPECIES
        FROM flora GROUP BY island,year ORDER BY island,year""",
        args.out/"allflora_moth_island_year.csv")
    export_query(dest, """SELECT island,
        COUNT(DISTINCT species) AS all_accepted_recorded_species,
        COUNT(DISTINCT CASE WHEN moth_positive_photo_ids>0 THEN species END) AS moth_positive_recorded_species,
        SUM(photo_ids) AS all_accepted_photo_ids,
        SUM(moth_positive_photo_ids) AS moth_positive_photo_ids,
        SUM(observation_uuids) AS all_accepted_observation_uuids
        FROM flora GROUP BY island ORDER BY island""",
        args.out/"allflora_moth_island_all_years.csv")
    dest.close()
    (args.out/"source_receipt.json").write_text(json.dumps({
        "schema":"ALAN_full93m_WCVP_allflora_moth_visitor_positive_same_source_denominator_v1",
        "source_totals":dict(totals),
        "same_source_moth_positive_species_fraction_IS_NOT_obligate_nocturnal_dependence":True,
        "missing_GloBI_record_is_NOT_moth_absence":True,
        "nonrandom_iNaturalist_species_composition_not_true_flora_prevalence":True,
        "survey_effort_and_species_detection_not_corrected":True,
        "historic_light_slope_not_fitted":True},sort_keys=True,indent=2)+"\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="mode",required=True)
    a=sub.add_parser("group")
    for k in ("archives","registry","receipt","wcvp","names","globi","coast","db"):
        a.add_argument("--"+k,type=Path,required=True)
    a.add_argument("--group",type=int,required=True)
    a.add_argument("--chunksize",type=int,default=40000)
    b=sub.add_parser("collect")
    b.add_argument("--groups",type=Path,required=True)
    b.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if args.mode=="group":
        require(args.chunksize>0,"Chunksize must be positive")
        count_group(args)
    else:
        collect(args)


if __name__=="__main__":
    main()
