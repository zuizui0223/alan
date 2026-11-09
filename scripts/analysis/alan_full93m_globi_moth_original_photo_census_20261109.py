#!/usr/bin/env python3
"""Census all 93,209,351 ORIGINAL iNaturalist photo IDs by moth-guild evidence.

Source 14 disjoint SHA-signed observation-hash photo tiers x 16 groups x 8
subshards (224 archives); every original photo ID is enumerated EXACTLY once.
Taxonomy uses 2026-10-19 frozen iNat->WCVP accepted-plant match, followed by
positive *study-associated* GloBI records from its signed 2026-11-05 cache:
 (1) source reports a moth-family flower visit OR pollinates relation;
 (2) GloBI "pollinates" relation involving moth, not actual verified service;
 (3) published Walton 2020 exact DOI explicitly sampled nocturnal moth pollen
     carriage, NOT effective stigma deposition or local moth dependence.
Butterfly records are a positive ancillary flag, NOT proven moth absence.

Read 5.8m original photo metadata records per group in bounded gzip streams;
never download any JPEG, modify source samples, or create white/colour data.
Compute exact original eligible/photo ID opportunities and selected observation
UUID counts, grouped by sampling tier, recorded year and available file/coords.
Original source 93m indexes 42,111 nominal plant taxa and does not cover all
biota; unmatched accepted taxa have UNKNOWN moth status, never "not moth".
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

REGISTRY_SCHEMA="mois_complete_full93m_original_photo_target_registry_v1"
SCHEMA="ALAN_full93209351_original_photo_GloBI_positive_moth_evidence_exact_census_v1"
REGISTRY_SHA="4ed5139e510ba8442d4cdb03f2b769d94fac942832b1d544d7a7925bc431255b"
NAMES_SHA="dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2"
ACCEPTED_SHA="9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3"
GLOBI_SHA="a6b97aa0941215d2f645e4c4fbf19490492b93680a7d7c0502814c7a66902994"
PHOTO_TOTAL=93_209_351
SOURCE_TAXA=42_111
N_TIERS=14
N_GROUPS=16
PARTS_PER_GROUP=8
COLUMNS=("photo_id","observation_uuid","inat_taxon_id","observer_id",
         "photographer_id","observed_on","latitude","longitude",
         "positional_accuracy","photo_license","image_extension",
         "frozen_loader_supported","valid_observed_date","valid_geopoint",
         "observation_inclusion_score","tier")
FLAGS=("moth_visit_positive","moth_database_pollinates_positive",
       "walton2020_nocturnal_pollen_transport_positive",
       "butterfly_visit_positive")
MASKS=(1,2,4,8)
SUPPORTED_YEAR={"2019","2021","2023","2024","2025"}

def require(ok,msg):
    if not bool(ok):
        raise ValueError("ALAN_FULL93M_MOTH_SOURCE: "+str(msg))
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):h.update(part)
    return h.hexdigest()

def evidence_map(wcvp,names,globi):
    require(sha(wcvp)==NAMES_SHA and sha(names)==ACCEPTED_SHA and
            sha(globi)==GLOBI_SHA,
            "Signed original WCVP accepted names or GloBI source changed")
    mapping=pd.read_csv(wcvp,compression="gzip",dtype=str,
        usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"],
        keep_default_na=False)
    assert mapping.inat_taxon_id.is_unique
    mapping=mapping.loc[mapping.wcvp_match_state.eq(
        "unique_exact_name_to_accepted_id")].copy()
    accepted=pd.read_csv(names,compression="gzip",dtype=str,
        usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"],
        keep_default_na=False).drop_duplicates()
    require(accepted.wcvp_accepted_plant_name_id.is_unique,
            "One accepted WCVP ID has multiple botanical names")
    mapping=mapping.merge(accepted,on="wcvp_accepted_plant_name_id",
                          how="left",validate="many_to_one")
    require(mapping.inat_taxon_id.is_unique,
            "Inat taxon ID repeated after botanical accepted-name mapping")
    g=pd.read_csv(globi,compression="gzip",dtype=str,keep_default_na=False)
    require(len(g)==61967 and
            set(["plant_name","guild","relation","study_title"])<=set(g) and
            set(g.guild)=={"moth_family_visit","butterfly_family_visit"} and
            set(g.relation)=={"visitsFlowersOf","pollinates"},
            "Source original study-linked moth versus butterfly record fields changed")
    groups={}
    for name,f in g.groupby("plant_name",observed=True):
        moth=f.loc[f.guild.eq("moth_family_visit")]
        butterfly=f.loc[f.guild.eq("butterfly_family_visit")]
        mask=(1 if len(moth)>0 else 0)
        if moth.relation.eq("pollinates").any():mask|=2
        if moth.study_title.str.contains("rsbl.2019.0877",regex=False,
                                           case=False,na=False).any():mask|=4
        if len(butterfly)>0:mask|=8
        groups[str(name)]=mask
    by_taxon={}
    botanical_unmatched=[]
    for row in mapping.itertuples(index=False):
        if not row.wcvp_accepted_name or pd.isna(row.wcvp_accepted_name):
            botanical_unmatched.append(str(row.inat_taxon_id))
            continue
        by_taxon[str(row.inat_taxon_id)]=groups.get(str(row.wcvp_accepted_name),0)
    catalog={
       "source_GloBI_records":len(g),
       "source_positive_moth_plant_names":sum(bool(v&1) for v in groups.values()),
       "source_moth_pollinates_plant_names":sum(bool(v&2) for v in groups.values()),
       "source_Walton2020_pollen_transport_plant_names":sum(bool(v&4) for v in groups.values()),
       "source_positive_butterfly_plant_names":sum(bool(v&8) for v in groups.values()),
       "source_WCVP_resolved_iNat_taxon_IDs":len(by_taxon),
       "source_WCVP_iNat_taxa_with_any_positive_moth_visit":
           sum(bool(v&1) for v in by_taxon.values()),
       "source_WCVP_iNat_taxa_with_moth_pollinates":
           sum(bool(v&2) for v in by_taxon.values()),
       "source_WCVP_iNat_taxa_with_Walton":
           sum(bool(v&4) for v in by_taxon.values()),
       "source_WCVP_IDs_without_accepted_binomial":len(botanical_unmatched),
       "zero_mask_not_moth_absence":True,
    }
    return by_taxon,catalog

def source_registry(registry,group):
    require(sha(registry)==REGISTRY_SHA,
            "93.2m original photo target registry byte SHA changed")
    r=json.loads(Path(registry).read_text())
    require(r.get("schema")==REGISTRY_SCHEMA and
            r.get("exact_original_photo_IDs_in_all_disjoint_tiers")==PHOTO_TOTAL and
            r.get("original_named_taxa_from_full_source")==SOURCE_TAXA and
            r.get("source_only_NO_photo_pixels_or_hue_measured") is True and
            len(r["published_tier_archive_SHA256"])==N_TIERS*N_GROUPS,
            "Not the complete 93,209,351 signed original photo-only source")
    return r

def classify_original_photo_rows(rows,taxon_map,stats,seen_obs,tier):
    """One signed archive part at a time, no 93m-photo or UUID RAM materialization."""
    counted=0
    for row in rows:
        counted+=1
        require(row["tier"]==str(tier) and row["photo_id"].isdecimal() and
                row["observation_uuid"] and row["inat_taxon_id"],
                "Original photo/observation/tier source identity invalid")
        t=row["inat_taxon_id"]
        mask=taxon_map.get(t)
        stats["all_source_photos"]+=1
        if mask is None:
            stats["original_photos_unknown_WCVP_accepted_name"]+=1
            continue
        stats["original_photos_accepted_WCVP_annotated"]+=1
        if mask==0:
            # GloBI absence is *not* an observed biological guild-negative.
            stats["accepted_taxon_photos_without_positive_record"]+=1
            continue
        valid=row["valid_observed_date"]=="true"
        geo=row["valid_geopoint"]=="true"
        loader=row["frozen_loader_supported"]=="true"
        eligible=valid and geo and loader
        year=row["observed_on"][:4]
        uuid=row["observation_uuid"]
        prior=seen_obs.get(uuid)
        if prior is None:
            seen_obs[uuid]=mask
            for k,bit in zip(FLAGS,MASKS):
                if mask&bit:stats["independent_source_observation_UUIDs_"+k]+=1
        else:
            require(prior==mask,
                    "An actual original observation UUID has contradictory botanical taxa")
        for k,bit in zip(FLAGS,MASKS):
            if mask&bit:
                stats["original_photo_IDs_"+k]+=1
                if loader:stats["source_image_format_supported_"+k]+=1
                if valid:stats["source_valid_observed_date_"+k]+=1
                if geo:stats["source_valid_photo_geopoint_"+k]+=1
                if eligible:
                    stats["source_RGB_eligible_original_photo_IDs_"+k]+=1
                    if year in SUPPORTED_YEAR:
                        stats["source_RGB_eligible_ALAN_5years_"+k]+=1
                if year in SUPPORTED_YEAR:
                    stats["source_observation_year_2019_2021_2023_2024_2025_"+k]+=1
        if mask&1 and mask&8:
            stats["moth_and_butterfly_both_positive_original_photo_IDs"]+=1
        if mask&1:
            stats["moth_source_original_year_"+year]+=1
    return counted

def scan_group(archives,registry,receipt,wcvp,names,globi,group,out):
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite signed moth source census")
    r=source_registry(registry,group)
    g=json.loads(Path(receipt).read_text())
    require(g.get("schema")=="mois_onepass_full93m_disjoint_photo_target_tiers_v1" and
            g.get("source_snapshot_SHA256")==r["full_source_manifest_sha256"] and
            len(g.get("parts",[]))==PARTS_PER_GROUP and
            set(int(x["part"]) for x in g["parts"])==
              set(range(group*PARTS_PER_GROUP,(group+1)*PARTS_PER_GROUP)) and
            g.get("images_downloaded")==0 and g.get("colour_pixels_measured")==0,
            "Original 93m signed group receipts do not match the global ledger")
    taxon_map,source_catalog=evidence_map(wcvp,names,globi)
    totals=Counter()
    tiers=[]
    for tier in range(N_TIERS):
        filename=f"tier-{tier:02d}-group-{group:02d}.tar"
        source=Path(archives)/filename
        require(source.is_file() and sha(source)==r["published_tier_archive_SHA256"][filename],
                "Signed original photo tier archive SHA changed: "+filename)
        stats=Counter()
        seen_obs={}
        with tarfile.open(source,mode="r") as tar:
            names_present={x.name:x for x in tar.getmembers() if x.isfile()}
            expected={f"tier-{tier:02d}/part-{part:03d}.csv.gz"
                    for part in range(group*PARTS_PER_GROUP,(group+1)*PARTS_PER_GROUP)}
            require(set(names_present)==expected,
                    "One of eight original photo partitions missing or duplicated in "+filename)
            for path in sorted(expected):
                with gzip.GzipFile(fileobj=tar.extractfile(names_present[path])) as gzin:
                    with io.TextIOWrapper(gzin,encoding="utf-8",newline="") as f:
                        reader=csv.DictReader(f)
                        require(tuple(reader.fieldnames or ())==COLUMNS,
                                "Original full 93m source CSV schema changed")
                        classify_original_photo_rows(reader,taxon_map,stats,seen_obs,tier)
        require(stats["all_source_photos"]==g["counts_by_disjoint_tier"][str(tier)],
                "Original tier/group photo count not conserved")
        d={"tier":tier,"group":group,**dict(stats)}
        tiers.append(d)
        totals.update(stats)
        print("ALAN_FULL93M_SOURCE_MOTH_GLOBI_CENSUS_ORIGINAL_TIER",
              json.dumps({"group":group,"tier":tier,
                  "source_photos":stats["all_source_photos"],
                  "moth_visiting_species_photos":stats["original_photo_IDs_moth_visit_positive"],
                  "moth_pollinates_photos":stats["original_photo_IDs_moth_database_pollinates_positive"],
                  "direct_Walton2020_nocturnal_pollen_photos":
                      stats["original_photo_IDs_walton2020_nocturnal_pollen_transport_positive"]},
                  sort_keys=True),flush=True)
    require(totals["all_source_photos"]==g["original_photo_IDs_in_this_group"],
            "All 14 original source photo tiers are not summed exactly once")
    out.mkdir(parents=True)
    by_tier_file=out/f"ALAN_full93m_group{group:02d}_moth_guild_positive_tiers.csv"
    pd.DataFrame(tiers).fillna(0).to_csv(by_tier_file,index=False)
    result={
        "schema":SCHEMA,
        "source_group":group,
        "actual_source_original_photo_IDs_scanned":totals["all_source_photos"],
        "total_global_source_original_photo_IDs":PHOTO_TOTAL,
        "original_source_WCVP_GloBI_catalog":source_catalog,
        "per_group_original_photo_moth_evidence_counts":dict(totals),
        "this_group_by_14_immutably_disjoint_original_photo_tiers":tiers,
        "SHA256":{"registry":REGISTRY_SHA,"accepted":ACCEPTED_SHA,
                  "wcvp":NAMES_SHA,"globi":GLOBI_SHA,
                  "signed_group_receipt":sha(receipt),
                  "by_tier_output":sha(by_tier_file)},
        "moth_visited_is_positive_record_NOT_effective_or_exclusive_pollination":True,
        "GloBI_unrecorded_moth_does_not_establish_day_pollination_or_no_moths":True,
        "unknown_WCVP_taxa_not_excluded_as_biological_moth_absence":True,
        "actual_original_rgb_pixels_downloaded_or_inferred":0,
        "original_flower_colour_ROI9_unmeasured":True,
        "no_individual_photo_IDs_or_original_observation_UUIDs_exported":True,
    }
    (out/f"ALAN_full93m_group{group:02d}_moth_positive_receipt.json").write_text(
        json.dumps(result,sort_keys=True,indent=2)+"\n")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("archives","registry","receipt","wcvp","names","globi","out"):
        p.add_argument("--"+flag,type=Path,required=True)
    p.add_argument("--group",type=int,required=True)
    a=p.parse_args()
    require(0<=a.group<N_GROUPS,"Expected 16 independent group source archives")
    scan_group(a.archives,a.registry,a.receipt,a.wcvp,a.names,a.globi,a.group,a.out)

if __name__=="__main__":
    main()
