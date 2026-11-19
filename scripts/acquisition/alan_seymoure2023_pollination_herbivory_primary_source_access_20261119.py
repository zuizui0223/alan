#!/usr/bin/env python3
"""Publication-level public fulltext-access census for ALAN pollination/herbivory.

Study selection is fixed by source 2023 review Figshare original XLSX,
not by observed study direction. Endpoints are NOT exposure estimates.
Metadata URLs and DOI links never prove a lighting experiment or plant fitness.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from scripts.analysis.alan_seymoure2023_true_workbook_coding_audit_20261119 import (
    SOURCE_SHA256,SOURCE_URL,SOURCE_FILE_ID,SOURCE_NAME,get,source_rows
)
from scripts.analysis.alan_seymoure2023_original_direction_study_units_20261119 import (
    _source_scientific_rows,study_key,doi_norm
)

PRIORITY=("pollination","herbivory")
USER_AGENT="ALAN-biotic-source-discovery/2026-11-19 (scientific metadata)"
MAX_REPLY=1_000_000


def get_json(url):
    req=urllib.request.Request(url,headers={"User-Agent":USER_AGENT,
                                            "Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=12) as r:
        b=r.read(MAX_REPLY+1)
        if len(b)>MAX_REPLY:
            raise ValueError("Metadata too long")
        return json.loads(b)


def discover(meta_rows, query=get_json, strict=True):
    scientific,_=_source_scientific_rows(meta_rows,strict=strict)
    priority=[r for r in scientific if r["interaction_type"] in PRIORITY]
    grouped=defaultdict(list)
    for r in priority:
        grouped[study_key(r)].append(r)
    if strict and (len(priority)!=18 or len(grouped)!=17):
        raise ValueError(f"Original plant interaction review study-grain drift {len(priority)}/{len(grouped)}")
    status=Counter()
    category=Counter()
    resultrows=[]
    for paper, rows in sorted(grouped.items()):
        types=sorted(set(r["interaction_type"] for r in rows))
        category["+".join(types)]+=1
        rawdoi=next((doi_norm(r.get("DOI") or r.get("DOI Link")) for r in rows
                     if doi_norm(r.get("DOI") or r.get("DOI Link"))),"")
        if not rawdoi:
            resultrows.append({"has_DOI":False,"crossref":"NO_DOI","europepmc":"NO_DOI",
                               "original_method_any":any(r.get("method") for r in rows),
                               "interaction_types":types})
            status["no_original_DOI"]+=1
            continue
        entry={"has_DOI":True,"original_method_any":any(r.get("method") for r in rows),
               "interaction_types":types,
               "crossref":"NOT_RETRIEVED","europepmc":"NOT_RETRIEVED",
               "crossref_has_resource_links":False,
               "europepmc_fulltext_openaccess_reported":False,
               "europepmc_pmcid_present":False}
        # Public metadata only; do not parse primary study design from abstract
        # keywords or treat crossref links as open-access proofs.
        try:
            meta=query("https://api.crossref.org/works/"+urllib.parse.quote(rawdoi,safe=""))
            item=meta.get("message",{}) if isinstance(meta,dict) else {}
            if str(item.get("DOI","")).casefold().strip()!=rawdoi:
                entry["crossref"]="SOURCE_DOI_MISMATCH"
            else:
                entry["crossref"]="DOI_VERIFIED"
                entry["crossref_has_resource_links"]=bool(item.get("link"))
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError,KeyError,json.JSONDecodeError) as ex:
            entry["crossref"]="SOURCE_STOP_"+type(ex).__name__
        try:
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/search?"+urllib.parse.urlencode({
                "query":f"DOI:{rawdoi}","format":"json","pageSize":5
            })
            meta=query(url)
            hits=meta.get("resultList",{}).get("result",[])
            exact=[x for x in hits if doi_norm(x.get("doi",""))==rawdoi]
            if exact:
                entry["europepmc"]="EXACT_DOI_RECORD"
                entry["europepmc_fulltext_openaccess_reported"]=any(x.get("isOpenAccess")=="Y" for x in exact)
                entry["europepmc_pmcid_present"]=any(bool(x.get("pmcid")) for x in exact)
            else:
                entry["europepmc"]="DOI_NOT_INDEXED"
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError,KeyError,json.JSONDecodeError) as ex:
            entry["europepmc"]="SOURCE_STOP_"+type(ex).__name__
        status["crossref_"+entry["crossref"]]+=1
        status["europepmc_"+entry["europepmc"]]+=1
        if entry["europepmc_fulltext_openaccess_reported"]:
            status["europepmc_DOI_exact_OA_reported"]+=1
        if entry["europepmc_pmcid_present"]:
            status["europepmc_DOI_exact_PMCID"]+=1
        resultrows.append(entry)
    total=len(resultrows)
    output={
       "schema":"alan_20261119_original_17_pollination_herbivory_paper_public_fulltext_coverage_v1",
       "fixed_original_source_SHA256":SOURCE_SHA256,
       "selection":"All original Seymoure 2023 pollination+herbivory investigations, independent of coded effect directions",
       "original_priority_interaction_rows":len(priority),
       "unique_original_publication_keys":total,
       "original_interaction_group_paper_distribution":dict(sorted(category.items())),
       "original_method_documented_priority_publications":sum(e["original_method_any"] for e in resultrows),
       "original_source_DOI_present_publications":sum(e["has_DOI"] for e in resultrows),
       "metadata_access_status":dict(sorted(status.items())),
       "original_2023_review_total_papers":79,
       "original_original_lighting_treatment_clusters_verified":0,
       "meta_eligible_original_effects_with_cluster_SE_verified":0,
       "scientific_stop":{
         "metadata_DOI_retrieval_is_not_original_experimental_source_reanalysis":True,
         "crossref_has_file_link_is_not_fulltext_OA_proof":True,
         "EuropePMC_DOI_not_indexed_is_not_evidence_fulltext_unavailable":True,
         "no_auto_classification_of_randomized_light_experiment":True,
         "never_derive_plant_reproductive_effect_from_simplified_review_code":True,
         "source_paper_DOI_and_original_rows_not_exported":True,
         "requires_verified_treatment_IDs_and_plant_or_insect_endpoint_variance":True,
       }
    }
    return output


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    meta=json.loads(get(SOURCE_URL))
    xs=[x for x in meta.get("files",[]) if x.get("id")==SOURCE_FILE_ID and x.get("name")==SOURCE_NAME]
    if len(xs)!=1:
        raise ValueError("Original Figshare author source drift")
    b=get(xs[0]["download_url"])
    if hashlib.sha256(b).hexdigest()!=SOURCE_SHA256:
        raise ValueError("Actual original source byte hash mismatch")
    h,rows=source_rows(b)
    if len(h)!=20:
        raise ValueError("Original source table schema drift")
    out=discover(rows,strict=True)
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_PLANT_INTERACTION_17_PAPER_SOURCE_ACCESS_COUNTS",out["metadata_access_status"])
    print("ALAN_PLANT_INTERACTION_17_PAPER_DOI",out["original_source_DOI_present_publications"])
    print("ALAN_PLANT_INTERACTION_17_PAPER_ORIGINAL_METHOD",out["original_method_documented_priority_publications"])


if __name__=="__main__":
    main()
