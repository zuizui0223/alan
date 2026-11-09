#!/usr/bin/env python3
"""Extract EVERY signed original photo ID for published GloBI moth-visited plants.

Input is original public MOIS 93,209,351 distinct-photo SHA registry,
14 independent observation-hash tiers x 16 geographical-independent source
groups, each containing 8 original CSV.gz photo+UUID partitions. Source GloBI
and WCVP are original signed evidence; this program never fetches RGB images,
declares observed flower petals, infers moth visitation from petal features,
or interprets missing GloBI edges as absence of moth pollination.

Output: one gzip CSV per tier per original source group containing just the
documented moth-linked original photo IDs, true observation UUID, source plant
taxon ID, observed date/location, license/format and positive evidence flags.
This emits ALL matched original photo IDs, not a probabilistic subsample.
Exact expected 34,346,428 IDs is asserted only by the separate 16-group
SHA-verified ledger collector after every original archive is processed.
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
import tarfile
import pandas as pd

from scripts.analysis import alan_full93m_globi_moth_original_photo_census_20261109 as global_source

GROUPS=16
TIERS=14
N_PHOTOS=93209351
EXPECTED_MOTH_COUNT=34346428
SOURCE_TIER_FIELDS=global_source.COLUMNS
OUTPUT_FIELDS=("photo_id","observation_uuid","inat_taxon_id","observed_on",
    "latitude","longitude","photo_license","image_extension",
    "frozen_loader_supported","valid_observed_date","valid_geopoint",
    "tier","moth_database_pollinates","walton2020_nocturnal_pollen_carriage",
    "butterfly_record_positive")
SCHEMA="ALAN_93m_original_GloBI_moth_positive_complete_source_photo_IDs_v1"

def require(ok,message):
    if not bool(ok):
        raise ValueError("ALAN_FULL93M_MOTH_PHOTO_ID_EXTRACTION: "+message)

def sha(path):
    return global_source.sha(path)

def extract_rows(rows, taxon_map, writers, source_tier,
                 stats, distinct_obs, positive_counts):
    """Stream source CSV rows; no 93m-ID in-memory state."""
    n=0
    for r in rows:
        n+=1
        require(r["tier"]==str(source_tier) and r["photo_id"].isdecimal()
                and bool(r["observation_uuid"]) and bool(r["inat_taxon_id"]),
                "One original source row is not an immutable photo/UUID/taxon")
        evidence_mask=taxon_map.get(r["inat_taxon_id"],0)
        if not evidence_mask&1:
            continue
        row={field:r[field] for field in OUTPUT_FIELDS[:-3]}
        row["moth_database_pollinates"]="true" if evidence_mask&2 else "false"
        row["walton2020_nocturnal_pollen_carriage"]="true" if evidence_mask&4 else "false"
        row["butterfly_record_positive"]="true" if evidence_mask&8 else "false"
        writers.writerow(row)
        stats["moth_positive_original_photo_IDs"]+=1
        if evidence_mask&2:stats["moth_database_pollinates_photo_IDs"]+=1
        if evidence_mask&4:stats["walton_nocturnal_pollen_carriage_photo_IDs"]+=1
        if evidence_mask&8:stats["also_butterfly_positive_photo_IDs"]+=1
        if (r["valid_observed_date"]=="true"
            and r["valid_geopoint"]=="true"
            and r["frozen_loader_supported"]=="true"):
            stats["true_source_RGB_eligible_moth_photo_IDs"]+=1
        if r["observed_on"][:4] in ("2019","2021","2023","2024","2025"):
            stats["same_5year_moth_original_photo_IDs"]+=1
        # Within each disjoint source tier, multiple photos may legitimately
        # belong to the same original observation: do NOT count these twice.
        distinct_obs.add(r["observation_uuid"])
        positive_counts[r["inat_taxon_id"]]+=1
    return n

def run(archives,registry,receipt,wcvp,names,globi,group,out):
    archives,registry,receipt,out=map(Path,(archives,registry,receipt,out))
    require(not out.exists(),"Refuse overwrite of original moth photo-ID extraction")
    require(0<=group<GROUPS,"Original 16 source group range")
    original=global_source.source_registry(registry,group)
    group_receipt=json.loads(receipt.read_text())
    require(group_receipt.get("schema")==
                "mois_onepass_full93m_disjoint_photo_target_tiers_v1" and
            group_receipt.get("source_snapshot_SHA256")==
                original["full_source_manifest_sha256"] and
            len(group_receipt.get("parts",[]))==8 and
            set(p["part"] for p in group_receipt["parts"])==
                set(range(8*group,8*(group+1))),
            "Original MOIS 93m group receipt or 128 source partition verification drift")
    taxon_map,catalog=global_source.evidence_map(wcvp,names,globi)
    (out/"moth_original_photo_IDs").mkdir(parents=True)
    grouped=[]
    source_total=0
    for tier in range(TIERS):
        filename=f"tier-{tier:02d}-group-{group:02d}.tar"
        file=archives/filename
        require(file.is_file() and
                sha(file)==original["published_tier_archive_SHA256"][filename],
                "Original photo ID target tier tar source SHA differs")
        export=out/"moth_original_photo_IDs"/(
            f"moth-positive-original-tier-{tier:02d}-group-{group:02d}.csv.gz")
        stats=Counter()
        seen=set()
        taxon_counts=Counter()
        n=0
        with gzip.open(export,"wt",encoding="utf-8",newline="",compresslevel=5) as dest:
            writer=csv.DictWriter(dest,fieldnames=OUTPUT_FIELDS,
                                  lineterminator="\n")
            writer.writeheader()
            with tarfile.open(file,"r:") as tar:
                members={m.name:m for m in tar.getmembers() if m.isfile()}
                expected={f"tier-{tier:02d}/part-{p:03d}.csv.gz"
                           for p in range(group*8,(group+1)*8)}
                require(set(members)==expected,
                        "Missing or repeated source part in 93m signed tier archive")
                for name in sorted(expected):
                    with gzip.GzipFile(fileobj=tar.extractfile(members[name])) as z:
                        reader=csv.DictReader(io.TextIOWrapper(
                            z,encoding="utf-8",newline=""))
                        require(tuple(reader.fieldnames or ())==SOURCE_TIER_FIELDS,
                                "Original 93m source photo/observation fields changed")
                        n+=extract_rows(reader,taxon_map,writer,tier,stats,seen,
                                        taxon_counts)
        require(n==int(group_receipt["counts_by_disjoint_tier"][str(tier)]),
                "Original source photo ID count not fully scanned in signed tier")
        info={"tier":tier,"group":group,
              "all_original_photo_IDs_examined":n,
              "source_moth_positive_distinct_observation_UUIDs":len(seen),
              "source_moth_positive_distinct_original_iNat_taxa":len(taxon_counts),
              "moth_original_photo_file_SHA256":sha(export),
              "source_original_tier_sha256":sha(file),
              **dict(stats)}
        grouped.append(info)
        source_total+=n
        print("ALAN_FULL93M_MOTH_REAL_SOURCE_ID_TIER_EXTRACTED",
              json.dumps({k:info[k] for k in (
                "tier","group","all_original_photo_IDs_examined",
                "moth_positive_original_photo_IDs",
                "source_moth_positive_distinct_observation_UUIDs")},
                sort_keys=True),flush=True)
    require(source_total==int(group_receipt["original_photo_IDs_in_this_group"]),
            "Source 93m original photo IDs missing after 14 disjoint tiers")
    result={
       "schema":SCHEMA,
       "source_group":group,
       "all_original_photo_IDs_scanned":source_total,
       "GloBI_moth_positive_original_photo_IDs_extracted":sum(
           d.get("moth_positive_original_photo_IDs",0) for d in grouped),
       "GloBI_moth_positive_original_photo_IDs_with_supported_RGB_format_geopoint_date":sum(
           d.get("true_source_RGB_eligible_moth_photo_IDs",0) for d in grouped),
       "published_positive_species_source_catalog":catalog,
       "source_tier_rows":grouped,
       "global_all_original_photo_IDs":N_PHOTOS,
       "globally_expected_moth_original_photo_IDs":EXPECTED_MOTH_COUNT,
       "preserved_original_observation_UUIDs_in_output":True,
       "original_photo_IDs_not_repeated_across_signed_14x16_tiers":True,
       "missing_GloBI_interaction_not_moth_absence":True,
       "documented_moth_visitors_not_exclusive_moth_pollinators":True,
       "original_source_RGB_images_downloaded":0,
       "original_flower_ROI9_pixels_inferred":0,
       "no_native_flora_status_inferred":True,
       "original_source_shas":{
          "registry":global_source.REGISTRY_SHA,
          "WCVP_taxonomy":global_source.NAMES_SHA,
          "accepted_WCVP_names":global_source.ACCEPTED_SHA,
          "signed_GloBI_positive_interactions":global_source.GLOBI_SHA,
          "source_group_receipt":sha(receipt)}
    }
    (out/f"original_moth_photoID_group{group:02d}_verified_receipt.json"
      ).write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_FULL93M_MOTH_COMPLETE_GROUP_PHOTO_IDS_VERIFIED",
          json.dumps({"group":group,
                      "original_photos_scanned":source_total,
                      "moth_positive_original_photo_IDs":result[
                        "GloBI_moth_positive_original_photo_IDs_extracted"]},
                        sort_keys=True),flush=True)
    return result

def collect(receipts,out):
    out=Path(out)
    require(not out.exists(),"Cannot replace previously verified full moth source count")
    all_=[]
    for p in receipts:
        r=json.loads(Path(p).read_text())
        require(r["schema"]==SCHEMA and r["source_group"] in range(GROUPS),
                "Invalid completed original MOIS source photo ID group")
        all_.append(r)
    require(len(all_)==GROUPS and set(x["source_group"] for x in all_)==set(range(GROUPS)),
            "Cannot declare all 93m original photo IDs enumerated from partial source")
    total=sum(x["all_original_photo_IDs_scanned"] for x in all_)
    moth=sum(x["GloBI_moth_positive_original_photo_IDs_extracted"] for x in all_)
    eligible=sum(x["GloBI_moth_positive_original_photo_IDs_with_supported_RGB_format_geopoint_date"]
                 for x in all_)
    require(total==N_PHOTOS and moth==EXPECTED_MOTH_COUNT,
            "Entire original 93m photo-ID set was NOT exactly covered")
    tiers=[]
    for source in all_:
        tiers.extend(source["source_tier_rows"])
    require(len(tiers)==GROUPS*TIERS,
            "All 224 original source photo target tiers were not accounted for")
    d={
      "schema":"ALAN_full93209351_original_photo_moth_positive_exhaustive_ID_manifest_v1",
      "original_all_photo_IDs_exactly_inspected":total,
      "moth_positive_all_original_photo_IDs_verified_and_materialized":moth,
      "moth_positive_original_RGB_format_geopoint_date_eligible":eligible,
      "source_group_receipts_completed":GROUPS,
      "source_photo_tiers_verified":len(tiers),
      "source_GloBI_positive_species_not_moth_effectiveness":True,
      "all_filtered_original_photo_IDs_in_group_output_artifacts_not_one_big_photo_file":True,
      "raw_original_RGB_images_downloaded":0,
      "original_colour_ROI9_measured":0,
      "group_receipt_SHA256":{
         str(i["source_group"]):sha(receipts[all_.index(i)])
         for i in all_},
    }
    out.mkdir(parents=True)
    (out/"ALAN_full93m_moth_positive_all_original_photo_IDs_recovered_manifest.json"
     ).write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    pd.DataFrame(tiers).to_csv(
        out/"full93m_moth_original_photo_IDs_by_14_tiers_16_groups.csv",
        index=False)
    print("ALAN_FULL93209351_MOTH_POSITIVE_34346428_ALL_ORIGINAL_PHOTO_IDS_EXTRACTED",
          json.dumps(d,sort_keys=True),flush=True)
    return d

def main():
    p=argparse.ArgumentParser(description=__doc__)
    subs=p.add_subparsers(dest="mode",required=True)
    q=subs.add_parser("group")
    for a in ("archives","registry","receipt","wcvp","names","globi","out"):
        q.add_argument("--"+a,type=Path,required=True)
    q.add_argument("--group",type=int,required=True)
    q=subs.add_parser("collect")
    q.add_argument("--receipts",type=Path,nargs="+",required=True)
    q.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.mode=="group":
        run(a.archives,a.registry,a.receipt,a.wcvp,a.names,a.globi,a.group,a.out)
    else:
        collect(a.receipts,a.out)
if __name__=="__main__":main()
