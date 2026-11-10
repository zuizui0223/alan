#!/usr/bin/env python3
"""Audit authentic GBIF gecko taxon and 2026-08 frozen Parquet coverage.

No reptile occurrence presence is interpreted as local predation, consumption
of moths, or species native/invasive status. GBIF human photos reflect
detection/observer effort. The GBIF taxon *Gekkonidae* is NOT all geckos;
other families of Gekkota must be checked, and day-gecko Phelsuma is a
different behavioural guild.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.parse,urllib.request
from datetime import datetime, timezone
from pathlib import Path

API="https://api.gbif.org/v1"
YEARS=(2013,2025)
FAMILIES=("Gekkonidae","Phyllodactylidae","Sphaerodactylidae",
          "Eublepharidae","Diplodactylidae","Carphodactylidae","Pygopodidae")
GENERA=("Hemidactylus","Gehyra","Gekko","Tarentola","Phelsuma",
        "Lepidodactylus","Cyrtodactylus")
SPECIES=("Hemidactylus frenatus","Gekko japonicus")
CONTROLS=("Squamata","Scincidae")
RANK={**{x:"FAMILY" for x in FAMILIES},
      **{x:"GENUS" for x in GENERA},
      **{x:"SPECIES" for x in SPECIES},
      "Squamata":"ORDER","Scincidae":"FAMILY"}
COUNTRIES=("JP","TW","PH","ID","AU","NZ","MG","MU","SC","US","CU","DO","FJ","NC")
SCHEMA="ALAN_global_GBIF_gecko_predator_taxonomy_and_record_source_audit_v1"
SNAPSHOT_DOI="10.15468/dl.8yrbe7"
PREVIOUS_MANIFEST_SCHEMA="ALAN_2026_08_GBIF_actual_taxon_partition_public_S3_manifest_v1"

def require(value,reason):
    if not bool(value):raise ValueError("ALAN_REAL_GECKO_GBIF: "+str(reason))

def query(endpoint,params):
    url=API+endpoint+"?"+urllib.parse.urlencode(params,doseq=True)
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-gecko-moth-macroecology-verifiable-source/1.0",
        "Accept":"application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req,timeout=55) as r:
                require(r.status==200,"Live original GBIF API status changed")
                raw=r.read()
            return json.loads(raw),hashlib.sha256(raw).hexdigest(),url
        except Exception:
            if attempt==3:raise
            time.sleep(1+attempt)
    raise AssertionError("GBIF request did not complete")

def resolve(name,rank):
    data,digest,url=query("/species/match",{"name":name,"strict":"true"})
    if data.get("usageKey") and data.get("rank","").upper()==rank and data.get("matchType")!="NONE":
        return {"usageKey":int(data["usageKey"]),"rank":rank,
                "acceptedName":data.get("scientificName"),
                "sourceGBIF":url,"source_sha256":digest,"method":"strict_taxon_match"}
    data,digest,url=query("/species/search",{
        "q":name,"rank":rank,"status":"ACCEPTED","limit":100})
    accepted=[x for x in data.get("results",[])
              if str(x.get("canonicalName","")).lower()==name.lower()
              and str(x.get("rank","")).upper()==rank
              and str(x.get("taxonomicStatus","")).upper()=="ACCEPTED"]
    keys=sorted({int(x["key"]) for x in accepted if x.get("key")})
    require(len(keys)==1,
            f"No unique author-resolved taxon for {name}/{rank}; candidates {keys[:15]}")
    return {"usageKey":keys[0],"rank":rank,"acceptedName":name,
            "sourceGBIF":url,"source_sha256":digest,"method":"exact_accepted_search"}

def occurrence(taxon_key,year_facet=True,country=None):
    args=[("taxonKey",str(taxon_key)),("basisOfRecord","HUMAN_OBSERVATION"),
          ("occurrenceStatus","PRESENT"),("hasCoordinate","true"),
          ("year","2013,2025"),("limit","0")]
    if year_facet:args.extend([("facet","year"),("facetLimit","100")])
    if country:args.append(("country",country))
    r,digest,url=query("/occurrence/search",args)
    require(isinstance(r.get("count"),int) and r["count"]>=0,
            "No real indexed human-observation count returned")
    out={"source_record_count":r["count"],
         "source_query_url":url,"source_response_SHA256":digest}
    if year_facet:
        fv=next((x for x in r.get("facets",[]) if
                 str(x.get("field","")).lower()=="year"),None)
        require(fv is not None,"Live GBIF year facet not returned")
        byyear={str(x["name"]):int(x["count"]) for x in fv.get("counts",[])}
        require(sum(byyear.values())==r["count"],
                "Year-facet source sum not reconcilable with GBIF count")
        out["per_source_year"]=byyear
    return out

def manifest_reptiles(manifest):
    x=json.loads(Path(manifest).read_text())
    require(x["schema"]==PREVIOUS_MANIFEST_SCHEMA and
            x["GBIF_source_DOI"]==SNAPSHOT_DOI and
            x["lepidoptera_original_source_partitions"]==6,
            "Not the actual version-frozen 2026-08 GBIF source manifest")
    bins=x["genuine_lepidoptera_parquet_partition_manifests"] # source has only Lepidoptera
    # Original S3 listing audit stores only Lepidoptera bins in its output,
    # so a new reptile HEAD/list operation is required, NOT a guessed key.
    return {"original_Partial_manifest_cannot_represent_all_Reptilia":True,
            "original_cached_lepidoptera_partitions":len(bins),
            "requires_new_real_S3_reptile_partition_list":True}

def run(manifest,out):
    out=Path(out);require(not out.exists(),"Never overwrite an original gecko opportunity audit")
    taxa={};unknown={}
    for name,rank in RANK.items():
        try:
            taxa[name]={"source_taxonomy":resolve(name,rank)}
            taxa[name]["actual_GBIF_HUMAN_OBSERVATION_2013_2025"]=occurrence(
                taxa[name]["source_taxonomy"]["usageKey"])
            print("ALAN_REAL_GECKO_TAXON_SOURCE",json.dumps({
                "taxon":name,
                "key":taxa[name]["source_taxonomy"]["usageKey"],
                "GBIF_georeferenced_human_records":
                    taxa[name]["actual_GBIF_HUMAN_OBSERVATION_2013_2025"]["source_record_count"]},
                sort_keys=True),flush=True)
        except Exception as exc:
            unknown[name]=str(exc)[:500]
            print("ALAN_GECKO_AMBIGUOUS_SOURCE",json.dumps({
                "taxon":name,"source_error":unknown[name]},sort_keys=True),flush=True)
    require("Gekkonidae" in taxa and "Hemidactylus frenatus" in taxa and
            "Squamata" in taxa,
            "Mandatory live source gecko/human reptile counts unavailable")
    nations={}
    for country in COUNTRIES:
        ent={}
        for taxon in ("Gekkonidae","Hemidactylus frenatus","Scincidae"):
            if taxon in taxa:
                ent[taxon]=occurrence(
                    taxa[taxon]["source_taxonomy"]["usageKey"],
                    country=country,year_facet=False)["source_record_count"]
        nations[country]=ent
        print("ALAN_TRUE_GECKO_COUNTRY_HUMAN_COUNT",
              json.dumps({country:ent},sort_keys=True),flush=True)
    original_manifest=manifest_reptiles(manifest)
    result={
      "schema":SCHEMA,
      "asof_UTC":datetime.now(timezone.utc).isoformat(),
      "GBIF_2026_08_preview_DOI":SNAPSHOT_DOI,
      "GBIF_LIVE_taxon_counts_not_the_2026_08_frozen_snapshot":True,
      "GBIF_observation_year_window":[2013,2025],
      "actual_GBIF_reptile_and_gecko_taxa":taxa,
      "unresolved_taxonomic_source_keys":unknown,
      "real_country_level_counts_NOT_island_counts":nations,
      "source_Parquet_manifest_gap":original_manifest,
      "same_gbif_gecko_and_moth_photos_can_share_observer_effort":True,
      "absence_of_gecko_photograph_is_not_real_gecko_absence":True,
      "some_day_geckos_can_shift_to_artificial_night_foraging":True,
      "Hemidactylus_frenatus_not_automatically_invasive_everywhere":True,
      "taxon_family_guilds_not_established_nocturnal_or_prey_diet":True,
      "no_gecko_eating_moth_observation_or_pollen_transfer_observed":True,
      "no_ALAN_gecko_predation_effect_fitted":True}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("ALAN_GBIF_TRUE_GECKO_WORLDWIDE_SOURCE_PREFLIGHT",
          json.dumps({"taxa_resolved":len(taxa),"taxa_unresolved":len(unknown),
             "Gekkonidae_human_records":taxa["Gekkonidae"][
                "actual_GBIF_HUMAN_OBSERVATION_2013_2025"]["source_record_count"],
             "Hemidactylus_frenatus_human_records":taxa["Hemidactylus frenatus"][
                "actual_GBIF_HUMAN_OBSERVATION_2013_2025"]["source_record_count"]},
            sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.manifest,a.out)
