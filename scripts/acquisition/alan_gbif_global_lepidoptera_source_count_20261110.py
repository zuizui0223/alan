#!/usr/bin/env python3
"""Public global moth-versus-butterfly record opportunity census (NOT moth abundance).

Use only real GBIF API record counts, actual GBIF taxon keys, source statuses,
2013-2025 geo-referenced HUMAN_OBSERVATION, source-year breakdown/countries.
A moth is NOT a monophyletic GBIF taxon; a count for Lepidoptera minus
Papilionoidea is a documented *taxonomic remainder proxy*, not an actual moth
night-activity or ecological population estimate. GBIF human observations
are effort/detection-biased; at this stage we deliberately do NOT fit light.
No GBIF login or bulk occurrence download; bounded public count only.
"""
from __future__ import annotations
import argparse, hashlib, json, time, urllib.parse, urllib.request
from pathlib import Path
from datetime import datetime, timezone

BASE="https://api.gbif.org/v1"
RANK={"Lepidoptera":"ORDER","Papilionoidea":"SUPERFAMILY","Insecta":"CLASS",
      "Noctuidae":"FAMILY","Erebidae":"FAMILY","Geometridae":"FAMILY",
      "Sphingidae":"FAMILY","Crambidae":"FAMILY","Pyralidae":"FAMILY",
      "Tortricidae":"FAMILY","Saturniidae":"FAMILY","Notodontidae":"FAMILY",
      "Lasiocampidae":"FAMILY"}
YEARS=tuple(range(2013,2026))
COUNTRIES=("JP","TW","PH","ID","AU","NZ","PG","MG","MU","SC",
           "ES","PT","GR","IT","US","CU","DO","BS","FJ","NC")
SCHEMA="ALAN_GBIF_Lepidoptera_butterfly_moth_source_opportunity_global_v1"

def require(ok, why):
    if not ok:raise ValueError("ALAN_GBIF_MOTH_MACRO: "+str(why))

def get_json(path,query):
    url=BASE+path+"?"+urllib.parse.urlencode(query,doseq=True)
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-Global-Moth-Macro-Source-Audit/1.0 (academic ecological research)",
        "Accept":"application/json"})
    for i in range(4):
        try:
            with urllib.request.urlopen(req,timeout=60) as r:
                require(r.status==200,"Unexpected public GBIF status")
                raw=r.read()
                return json.loads(raw),hashlib.sha256(raw).hexdigest(),url
        except Exception:
            if i==3:raise
            time.sleep(1+2*i)
    raise AssertionError("Unexpected GBIF API empty")

def resolve():
    out={}
    for name,rank in RANK.items():
        data,sha,url=get_json("/species/match",{"name":name,"strict":"true"})
        require(data.get("rank","").upper()==rank and
                data.get("matchType") not in ("NONE",None) and
                data.get("usageKey") is not None,
                f"GBIF taxon key not verified for {name}: {str(data)[:220]}")
        out[name]={"key":int(data["usageKey"]),"rank":rank,
                   "source_accepted_name":data.get("scientificName"),
                   "source_api_url":url,"api_json_SHA256":sha}
    return out

def occ_count(taxon_key,facet=None,country=None,year_start=2013,year_end=2025,
              basis="HUMAN_OBSERVATION"):
    q=[("taxonKey",str(taxon_key)),("hasCoordinate","true"),
       ("year",f"{year_start},{year_end}"),("occurrenceStatus","PRESENT")]
    if basis:q.append(("basisOfRecord",basis))
    if country:q.append(("country",country))
    if facet:q.extend([("facet",facet),("facetLimit","300"),("limit","0")])
    result,checksum,url=get_json("/occurrence/search",q)
    require(isinstance(result.get("count"),int) and result["count"]>=0,
            "Not an actual GBIF occurrence count")
    src={"n_actual_GBIF_occurrence_records":int(result["count"]),
         "source_request":url,"source_response_SHA256":checksum}
    if facet:
        values=result.get("facets",[])
        source=next((v for v in values if v.get("field")==facet),None)
        require(source is not None,"Requested GBIF facet was not returned")
        src["facet"]=facet
        src["facet_values"]={str(x.get("name")):int(x.get("count",0))
                             for x in source.get("counts",[])}
    return src

def check_subtraction(big,small,why):
    require(big>=small,f"Taxonomic source counting inconsistent ({why})")
    return big-small

def run(dest):
    dest=Path(dest)
    require(not dest.exists(),"Never overwrite source count receipt")
    taxonomy=resolve()
    keys={k:v["key"] for k,v in taxonomy.items()}
    samples={}
    for field,name in (("lepidoptera_all","Lepidoptera"),
                       ("butterfly_papilionoidea","Papilionoidea"),
                       ("insects_all","Insecta")):
        samples[field]=occ_count(keys[name],facet="year")
    samples["lepidoptera_non_papilionoidea_record_remainder"]={
       "n_difference_not_independently_taxon_resolved":
          check_subtraction(samples["lepidoptera_all"]["n_actual_GBIF_occurrence_records"],
            samples["butterfly_papilionoidea"]["n_actual_GBIF_occurrence_records"],
            "Lepidoptera-Papilionoidea"),
       "interpretation":"taxonomic remainder proxy, contains nonmoth and unknown lower-rank records"}
    for name in list(RANK)[3:]:
        samples["family_"+name]=occ_count(keys[name])
    countries={}
    for country in COUNTRIES:
        a=occ_count(keys["Lepidoptera"],country=country)
        b=occ_count(keys["Papilionoidea"],country=country)
        countries[country]={
           "lepidoptera_georef_human_records":a["n_actual_GBIF_occurrence_records"],
           "butterfly_georef_human_records":b["n_actual_GBIF_occurrence_records"],
           "other_lepidoptera_remainder_records":check_subtraction(
               a["n_actual_GBIF_occurrence_records"],
               b["n_actual_GBIF_occurrence_records"],country),
           "source_requests":[a["source_request"],b["source_request"]]}
        print("ALAN_GBIF_TRUE_COUNTRY_LEPIDOPTERA_SOURCE",
              json.dumps({country:countries[country]},sort_keys=True),flush=True)
    byyear={}
    ly=samples["lepidoptera_all"]["facet_values"]
    by=samples["butterfly_papilionoidea"]["facet_values"]
    for y in YEARS:
        a=int(ly.get(str(y),0));b=int(by.get(str(y),0))
        byyear[str(y)]={"lepidoptera_records":a,"butterfly_records":b,
            "other_lepidoptera_remainder_records":check_subtraction(a,b,str(y))}
    rep={
        "schema":SCHEMA,"source":"GBIF LIVE public occurrence indexed API, not version-frozen snapshot",
        "query_time_UTC":datetime.now(timezone.utc).isoformat(),
        "taxon_provenance":taxonomy,"time_window":[YEARS[0],YEARS[-1]],
        "source_human_observations_only":True,"source_hasCoordinates":True,
        "source_occurrenceStatus_present":True,
        "global_source_counts":samples,"yearly_source_counts":byyear,
        "preselected_country_source_counts":countries,
        "not_exact_island_coverage":"country includes continental areas; no polygon assignment yet",
        "moth_paraphyletic_not_distinguishable_by_taxon_remainder":True,
        "known_moth_family_counts_are_overlapping_with_Lepidoptera_total_not_new_additive_observations":True,
        "repeated_photographs_not_independent_population_individuals":True,
        "no_moth_biomass_abundance_pollination_or_plant_fitness_assessed":True,
        "no_insect_count_denominator_adjusted_for_observer_effort":True,
        "no_causal_ALAN_effect_fitted":True,
        "gbif_large_scale_collection_requires_version_pinned_bulk_snapshot_or_authenticated_download":
            True
    }
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n")
    print("ALAN_GBIF_TRUE_WORLDWIDE_MOTH_SOURCE_COUNT",json.dumps({
        "global_lepidoptera":samples["lepidoptera_all"]["n_actual_GBIF_occurrence_records"],
        "global_butterfly":samples["butterfly_papilionoidea"]["n_actual_GBIF_occurrence_records"],
        "global_other_lepidoptera":samples["lepidoptera_non_papilionoidea_record_remainder"]
            ["n_difference_not_independently_taxon_resolved"],
        "source_years":len(YEARS),"source_countries":len(countries)},
        sort_keys=True),flush=True)
    return rep

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.out)
