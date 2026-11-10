#!/usr/bin/env python3
"""Bounded worldwide moth vs butterfly record opportunities, not moth abundance.

2013-2025 original georeferenced GBIF HUMAN_OBSERVATION. Explicit butterfly
families = six disjoint taxonomic sources; other Lepidoptera = remainder
NOT an exhaustive strictly nocturnal moth census or local pollination.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.parse,urllib.request
from datetime import datetime,timezone
from pathlib import Path

BASE="https://api.gbif.org/v1"
BUTTERFLY=("Papilionidae","Pieridae","Nymphalidae",
           "Lycaenidae","Riodinidae","Hesperiidae")
MOTH=("Noctuidae","Erebidae","Geometridae","Sphingidae",
      "Crambidae","Pyralidae","Tortricidae","Saturniidae",
      "Notodontidae","Lasiocampidae")
NAMES={"Lepidoptera":"ORDER","Insecta":"CLASS"}
NAMES.update({x:"FAMILY" for x in BUTTERFLY+MOTH})
YEARS=tuple(range(2013,2026))
COUNTRIES=("JP","TW","PH","ID","AU","NZ","PG","MG","MU","SC",
           "ES","PT","GR","IT","US","CU","DO","BS","FJ","NC")
SCHEMA="ALAN_GBIF_Lepidoptera_six_butterfly_families_worldwide_source_v2"

def require(x,m):
    if not x:raise ValueError("ALAN_GLOBAL_GBIF: "+str(m))

def fetch(endpoint,params):
    url=BASE+endpoint+"?"+urllib.parse.urlencode(params,doseq=True)
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-Global-Macroecology-Source-Audit/1.0",
        "Accept":"application/json"})
    for j in range(4):
        try:
            with urllib.request.urlopen(req,timeout=60) as response:
                require(response.status==200,"HTTP source status non200")
                raw=response.read()
            return json.loads(raw),hashlib.sha256(raw).hexdigest(),url
        except Exception:
            if j==3:raise
            time.sleep(j+1)
    raise AssertionError("Missing GBIF source")

def resolve(name,rank):
    d,digest,url=fetch("/species/match",{"name":name,"strict":"true"})
    if d.get("usageKey") and d.get("rank","").upper()==rank and d.get("matchType")!="NONE":
        return {"key":int(d["usageKey"]),"rank":rank,"source_url":url,
                "source_sha256":digest,"resolution":"strict_exact_name"}
    d,digest,url=fetch("/species/search",{"q":name,"rank":rank,
                                          "status":"ACCEPTED","limit":30})
    good=[i for i in d.get("results",[])
          if str(i.get("canonicalName","")).casefold()==name.casefold() and
          str(i.get("rank","")).upper()==rank and
          str(i.get("taxonomicStatus","")).upper() in ("ACCEPTED","DOUBTFUL")]
    keys={int(i["key"]) for i in good if i.get("key")}
    require(len(keys)==1,f"{name} has ambiguous GBIF accepted taxonomy keys {keys}")
    return {"key":next(iter(keys)),"rank":rank,"source_url":url,
            "source_sha256":digest,"resolution":"exact_rank_search"}

def count(key,facet=False,country=None):
    p=[("taxonKey",str(key)),("hasCoordinate","true"),
       ("year","2013,2025"),("occurrenceStatus","PRESENT"),
       ("basisOfRecord","HUMAN_OBSERVATION"),("limit","0")]
    if country:p.append(("country",country))
    if facet:p.extend([("facet","year"),("facetLimit","100")])
    d,digest,url=fetch("/occurrence/search",p)
    require(isinstance(d.get("count"),int) and d["count"]>=0,
            "Source GBIF occurrence count missing or malformed")
    result={"n_actual_records":int(d["count"]),"source":url,"sha256":digest}
    if facet:
        f=next((x for x in d.get("facets",[])
                if x.get("field","").lower()=="year"),None)
        require(f is not None,"GBIF source facet year missing")
        result["by_year"]={str(x["name"]):int(x["count"])
                           for x in f.get("counts",[])}
        require(sum(result["by_year"].values())==result["n_actual_records"],
                "GBIF year facet cannot reconcile source total")
    return result

def group(keys,names,facet=False,country=None):
    n=0;each={};yr={str(y):0 for y in YEARS}
    for x in names:
        r=count(keys[x]["key"],facet=facet,country=country)
        each[x]=r
        n+=r["n_actual_records"]
        if facet:
            for y in YEARS: yr[str(y)]+=int(r["by_year"].get(str(y),0))
    z={"n_actual_records":n,"source_taxon_counts":each}
    if facet:z["by_year"]=yr
    return z

def run(out):
    out=Path(out);require(not out.exists(),"No output overwrite")
    keys={}
    optional_moth_taxonomic_uncertainties={}
    for name,rank in NAMES.items():
        try:
            keys[name]=resolve(name,rank)
            print("ALAN_GBIF_TRUE_SOURCE_TAXON",json.dumps(
                {name:keys[name]},sort_keys=True),flush=True)
        except ValueError as err:
            if name not in MOTH:
                raise
            optional_moth_taxonomic_uncertainties[name]=str(err)
            print("ALAN_GBIF_OPTIONAL_MOTH_TAXON_AMBIGUOUS",
                  json.dumps({name:str(err)[:250]},sort_keys=True),flush=True)
    a=group(keys,("Lepidoptera",),True)
    b=group(keys,BUTTERFLY,True)
    c=group(keys,("Insecta",),True)
    actually_resolved_moth_families=tuple(x for x in MOTH if x in keys)
    require(len(actually_resolved_moth_families)>=5,
            "Too few independently verified true moth family names")
    m=group(keys,actually_resolved_moth_families,True)
    require(0<=b["n_actual_records"]<=a["n_actual_records"]<=c["n_actual_records"],
            "GBIF taxon subset counts cannot exceed order/class")
    countries={}
    for country in COUNTRIES:
        l=group(keys,("Lepidoptera",),country=country)
        f=group(keys,BUTTERFLY,country=country)
        require(f["n_actual_records"]<=l["n_actual_records"],
                "GBIF country butterfly count exceeds Lepidoptera")
        countries[country]={"lepidoptera_records":l["n_actual_records"],
             "six_butterfly_family_records":f["n_actual_records"],
             "other_lepidoptera_remainder":l["n_actual_records"]-f["n_actual_records"],
             "source_lepidoptera":l["source_taxon_counts"],
             "source_butterfly":f["source_taxon_counts"]}
        print("ALAN_GBIF_TRUE_COUNTRY_COUNTS",json.dumps({country:{
            k:v for k,v in countries[country].items() if not k.startswith("source_")
            }},sort_keys=True),flush=True)
    years={}
    for y in YEARS:
        k=str(y);big=a["by_year"][k];small=b["by_year"][k]
        require(big>=small,f"Global year {y} GBIF source families exceed order")
        years[k]={"lepidoptera":big,"six_butterfly_families":small,
                  "other_lepidoptera_remainder":big-small}
    out.parent.mkdir(parents=True,exist_ok=True)
    report={"schema":SCHEMA,
       "source_GBIF_live_index_not_version_frozen":True,
       "source_time_UTC":datetime.now(timezone.utc).isoformat(),
       "source_taxonomy_exact_rank":keys,
       "source_2013_to_2025_human_observation_only":True,
       "source_has_true_coordinates":True,
       "global_lepidoptera_records":a["n_actual_records"],
       "global_six_butterfly_family_records":b["n_actual_records"],
       "global_other_lepidoptera_remainder_records":a["n_actual_records"]-b["n_actual_records"],
       "global_partial_verified_moth_family_records":m["n_actual_records"],
       "actual_resolved_moth_families":list(actually_resolved_moth_families),
       "optional_moth_taxonomic_uncertainties":optional_moth_taxonomic_uncertainties,
       "global_insect_human_observation_records":c["n_actual_records"],
       "original_GBIF_global_taxon_sources":{"lepidoptera":a,"butterfly":b,
             "moth10_families":m,"insect":c},
       "global_13_actual_years":years,
       "20_country_counts_NOT_exact_islands":countries,
       "no_exact_island_assignment_yet":True,
       "other_lepidoptera_remainder_not_verified_moth_or_nocturnal":True,
       "global_records_not_capture_counts_or_observation_effort_corrected":True,
       "no_species_absence_or_effective_pollination_inferred":True,
       "no_global_light_biology_effect_estimated":True,
       "bulk_analysis_should_use_versioned_S3_parquet_or_authenticated_GBIF_download":True}
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("ALAN_GBIF_TRUE_WORLDWIDE_LEPIDOPTERA_SOURCE",
        json.dumps({k:v for k,v in report.items() if k.startswith("global_") and
                    isinstance(v,int)},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    args=p.parse_args();run(args.out)
