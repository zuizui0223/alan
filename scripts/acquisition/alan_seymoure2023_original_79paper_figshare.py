#!/usr/bin/env python3
"""Source-authenticated Seymoure et al. 2023 ALAN interaction evidence panel.

The publication reviews 79 articles and 81 investigations, not 81 randomized
independent experiments or comparable effect sizes. Recover published Figshare
supplement securely; export only aggregate workbook schema and original hashes
until a correct unit/variable audit is complete.
"""
import argparse
import hashlib
import io
import json
import urllib.error
import urllib.request
from pathlib import Path
from openpyxl import load_workbook

COLLECTION = 6837630
API = "https://api.figshare.com/v2"
PAPER = "10.1098/rstb.2022.0356"
COLLECTION_DOI = "10.6084/m9.figshare.c.6837630"
LIMIT = 8_000_000


def get(url):
    req = urllib.request.Request(url,headers={
        "User-Agent":"ALAN-79-study-publication-source-audit/1.0",
        "Accept":"application/json,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"
    })
    with urllib.request.urlopen(req, timeout=45) as f:
        data=f.read(LIMIT+1)
        if len(data)>LIMIT:
            raise ValueError("Source file exceeds safe preflight size")
        return data


def sheet_schema(b):
    if not b.startswith(b"PK\x03\x04"):
        raise ValueError("No original Excel file bytes")
    wb=load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    out=[]
    for ws in wb.worksheets:
        header=None
        header_index=None
        nonempty=0
        for ix,row in enumerate(ws.iter_rows(values_only=True),start=1):
            values=[str(x).strip() if x is not None else "" for x in row]
            if not any(values):continue
            nonempty+=1
            if header is None and sum(bool(x) for x in values)>=3:
                header_index=ix
                header=values
        out.append({
            "original_sheet_name":ws.title,
            "first_likely_header_row":header_index,
            "first_likely_headers_only":header or [],
            "source_nonempty_rows":nonempty,
            "no_original_study_rows_exported":True,
        })
    wb.close()
    return out


def audit():
    out={
        "schema":"alan_seymoure2023_79paper_81interaction_figshare_original_audit_v1",
        "source_article_DOI":PAPER,"source_supplement_collection_DOI":COLLECTION_DOI,
        "review_published_coverage":{
            "original_papers":79,"original_distinct_interactions":81,
            "original_predator_prey_investigations":56,"herbivory":11,"pollination":6,
            "detritivory":5,"frugivory_seed_dispersal":1,
            "parasite_host":1,"parasitoid_host":1,
            "paper_author_claims_not_independent_meta_effect_sizes":True,
            "selection_from_movement_ecology_database_not_complete_systematic_ALAN_review":True,
        },
        "source_status":"SOURCE_STOP",
        "figshare_articles":[],
        "retrieval_status":[],
        "raw_rows_or_site_IDs_committed":False,
        "author_original_row_effects_fitted":False,
        "cross_study_meta_effect_fitted":False,
    }
    base=f"{API}/collections/{COLLECTION}/articles"
    try:
        items=json.loads(get(base))
        if not isinstance(items,list):
            raise ValueError("Figshare collection response was not list")
        out["collection_number_of_articles"]=len(items)
        out["retrieval_status"].append({"type":"collection_api","status":"SUCCESS"})
    except (ValueError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as ex:
        out["retrieval_status"].append({"type":"collection_api","status":"SOURCE_STOP",
            "error_class":type(ex).__name__,
            **({"http_status":ex.code} if isinstance(ex,urllib.error.HTTPError) else {})})
        return out
    for entry in items[:12]:
        if not isinstance(entry,dict) or not isinstance(entry.get("id"),int):
            continue
        article_id=entry["id"]
        try:
            doc=json.loads(get(f"{API}/articles/{article_id}"))
            files=doc.get("files",[])
            a={"original_article_id":article_id,
               "title":str(doc.get("title",""))[:280],
               "source_files":[{
                   "original_file_id":f.get("id"),"name":f.get("name"),
                   "size":f.get("size"),"original_MD5":f.get("supplied_md5") or f.get("computed_md5")
               } for f in files]}
            out["figshare_articles"].append(a)
            out["retrieval_status"].append({"type":"article_metadata","article_id":article_id,"status":"SUCCESS"})
            for f in files:
                if not str(f.get("name","")).lower().endswith((".xlsx",".xls")):
                    continue
                size=f.get("size")
                if not isinstance(size,int) or size<100 or size>LIMIT:continue
                url=f.get("download_url")
                if not isinstance(url,str) or not url.startswith("https://"):continue
                try:
                    b=get(url)
                    file_replay={
                        "original_article_id":article_id,
                        "original_file_id":f.get("id"),"source_file_name":f["name"],
                        "source_excel_SHA256":hashlib.sha256(b).hexdigest(),
                        "source_byte_count":len(b),
                        "original_workbook_schema":sheet_schema(b),
                        "raw_study_rows_exported":False}
                    out.setdefault("source_verified_excels",[]).append(file_replay)
                    out["source_status"]="RAW_SUPPLEMENT_XLSX_SCHEMA_VERIFIED"
                except (ValueError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as ex:
                    out["retrieval_status"].append({
                        "type":"original_excel_download","original_file_id":f.get("id"),
                        "status":"SOURCE_STOP","error_class":type(ex).__name__,
                        **({"http_status":ex.code} if isinstance(ex,urllib.error.HTTPError) else {})})
        except (ValueError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as ex:
            out["retrieval_status"].append({
                "type":"article_metadata","article_id":article_id,
                "status":"SOURCE_STOP","error_class":type(ex).__name__})
    if out["source_status"]=="SOURCE_STOP" and out["figshare_articles"]:
        out["source_status"]="METADATA_ONLY"
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    p=a.parse_args()
    x=audit()
    p.out.parent.mkdir(parents=True,exist_ok=True)
    p.out.write_text(json.dumps(x,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_SEYMOURE79_SOURCE_STATUS",x["source_status"])
    print("ALAN_SEYMOURE79_METADATA",[(z["title"],[f["name"] for f in z["source_files"]])
                                      for z in x["figshare_articles"]])
    print("ALAN_SEYMOURE79_SCHEMA",x.get("source_verified_excels",[]))
    print("ALAN_SEYMOURE79_SOURCE_STOPS",x["retrieval_status"])


if __name__=="__main__":
    main()
