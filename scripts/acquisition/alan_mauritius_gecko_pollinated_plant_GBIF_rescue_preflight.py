#!/usr/bin/env python3
"""Find alternate *documented* GBIF occurrence opportunities for two endemic gecko-pollinated Mauritius plants.

Official accepted Kew names: Roussea simplex and Ruizia blackburniana
(syn. Trochetia blackburniana). Original 93m iNaturalist photo GSHHG
cohort had zero accepted-name photo sites for both. This script queries
only GBIF public aggregate counts and taxon/protocol metadata, **never
downloads or publishes coordinates of threatened plants**.

GBIF specimen geolocation does NOT prove extant wild plants or exact
light-exposure year; non-photographic surveys are not comparable to
photographic iNaturalist sampling without a source-method covariate.
"""
from __future__ import annotations
import argparse,json,time,hashlib,urllib.parse,urllib.request
from pathlib import Path
from datetime import datetime,timezone

SCHEMA="ALAN_Mauritius_reptile_pollination_original_GBIF_occurrence_rescue_source_v1"
API="https://api.gbif.org/v1"
TAXA=("Roussea simplex","Ruizia blackburniana","Trochetia blackburniana")

def request(api,params):
    url=API+api+"?"+urllib.parse.urlencode(params,doseq=True)
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-island-gecko-endemic-plant-source-audit/1.0",
        "Accept":"application/json"})
    for i in range(4):
        try:
            with urllib.request.urlopen(req,timeout=50) as r:
                if r.status!=200:raise ValueError("Original GBIF HTTP not 200")
                b=r.read()
            return json.loads(b),url,hashlib.sha256(b).hexdigest()
        except Exception:
            if i==3:raise
            time.sleep(i+1)
    raise RuntimeError("unavailable")

def run(out):
    out=Path(out)
    if out.exists():raise FileExistsError(out)
    results={}
    for plant in TAXA:
        m,url,sha=request("/species/match",{"name":plant,"strict":"true"})
        taxon={"source_query":url,
               "source_response_SHA256":sha,
               "source_matchType":m.get("matchType"),
               "source_rank":m.get("rank"),
               "source_usageKey":m.get("usageKey"),
               "source_scientificName":m.get("scientificName"),
               "source_status":m.get("status"),
               "source_acceptedUsageKey":m.get("acceptedUsageKey"),
               "source_confidence":m.get("confidence")}
        if m.get("rank")!="SPECIES" or not isinstance(m.get("usageKey"),int):
            results[plant]={"taxonomy":taxon,"status":"GBIF_taxon_unresolved"}
            continue
        count={}
        for label,params in (
            ("all_Mauritius_1992_2025",{"year":"1992,2025"}),
            ("georeferenced_Mauritius_1992_2025",{
                "year":"1992,2025","hasCoordinate":"true"}),
            ("georeferenced_Mauritius_2013_2025",{
                "year":"2013,2025","hasCoordinate":"true"}),
            ("human_observed_georeferenced_Mauritius_2013_2025",{
                "year":"2013,2025","hasCoordinate":"true",
                "basisOfRecord":"HUMAN_OBSERVATION"}),
            ("specimen_georeferenced_Mauritius_1992_2025",{
                "year":"1992,2025","hasCoordinate":"true",
                "basisOfRecord":"PRESERVED_SPECIMEN"}),
        ):
            dat,loc,dig=request("/occurrence/search",{
                "taxonKey":int(m["usageKey"]),"country":"MU",
                "occurrenceStatus":"PRESENT","limit":"0",**params})
            if not isinstance(dat.get("count"),int):
                raise ValueError("Actual GBIF occurrence source count missing")
            count[label]={"n_real_indexed_records":dat["count"],
                          "source_GBIF_API":loc,"source_SHA256":dig}
        results[plant]={"taxonomy":taxon,"source_record_counts":count,
            "source_precise_location_exported":False,
            "requires_Kew_synonym_taxon_identity_crosswalk":True}
        print("ALAN_MAURITIUS_PLANT_SOURCE_GBIF_ALTERNATIVE",
              json.dumps({plant:{"taxon_key":m.get("usageKey"),
                  "acceptedKey":m.get("acceptedUsageKey"),
                  "counts":{k:v["n_real_indexed_records"]
                            for k,v in count.items()}}},sort_keys=True),flush=True)
    report={"schema":SCHEMA,"asof_UTC":datetime.now(timezone.utc).isoformat(),
        "source_original_Kew_accepted_names":["Roussea simplex","Ruizia blackburniana"],
        "synonym_name_only":["Trochetia blackburniana"],
        "original_ALAN_93m_WCVP_accepted_photo_sites_for_both":0,
        "source_live_GBIF_index_not_locked_historic_snapshot":True,
        "genuine_Mauritius_census":results,
        "no_rare_plant_locations_released":True,
        "specimen_count_not_wild_live_population_or_pollen_transfer":True,
        "no_new_effective_pollination_guild_assumed_from_GBIF_taxonomy":True,
        "no_ALAN_effect_or_risk_estimated_from_metadata":True}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();run(a.out)
