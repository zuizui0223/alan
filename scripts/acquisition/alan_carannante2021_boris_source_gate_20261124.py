#!/usr/bin/env python3
"""Public source discovery gate for Carannante et al. 2021 aquatic adults.

DOI 10.1002/2688-8319.12053.
The article's source data availability points to legacy University of Bern
BORIS record https://boris.unibe.ch/151803/. The University moved its
repository software; search only public BORIS/BORIS Portal routes.

Study A: 4 genuine lit/dark paired trap placements, each on 6 nights, with
2 pairs at 3 m and 2 at 40 m from the Gürbe river edge. The 48 catches are
repeated subsamples; ~100k insect individuals are not independent ALAN units.
Study B: light traps at multiple distances without dark controls, not a
direct ALAN-vs-dark experiment.

A source landing record, thesis, published fitted statistics, or a DOI
citation is NOT author Study A individual trap-night data. Do not adopt a
new programme into a quantitative meta unless real public data bytes,
sample keys, provenance and original darkness comparisons are verified.
"""
import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlparse
from urllib.request import Request,urlopen

DOI="10.1002/2688-8319.12053"
EXPECTED_HOSTS={"boris.unibe.ch","boris-portal.unibe.ch","www.boris.unibe.ch"}
SOURCES=(
    "https://boris.unibe.ch/151803/",
    "https://boris-portal.unibe.ch/server/api/discover/search/objects?query=Carannante",
    "https://boris-portal.unibe.ch/server/api/discover/search/objects?query=151803",
    "https://boris-portal.unibe.ch/search?query=Carannante",
)
MAX_RESPONSE=450_000

class CandidateRefs(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls=[]
    def handle_starttag(self,tag,attrs):
        if tag not in ("a","link"):return
        d=dict(attrs)
        href=d.get("href") or ""
        if href.lower().split("?")[0].endswith((".csv",".tsv",".xls",".xlsx",".zip",".rds",".txt")):
            self.urls.append(href)

def extract_name_hints(obj,limit=15):
    names=[]
    def recurse(o):
        if len(names)>=limit:return
        if isinstance(o,dict):
            for k,v in o.items():
                if k in {"name","title","dc.title","label"} and isinstance(v,str):
                    if "carannante" in v.lower() or "aquatic insects" in v.lower():
                        names.append(v[:150])
                elif isinstance(v,(dict,list)):recurse(v)
        elif isinstance(o,list):
            for v in o:recurse(v)
    recurse(obj)
    return names

def probe(url):
    u=urlparse(url)
    if u.scheme!="https" or u.hostname not in EXPECTED_HOSTS:
        raise ValueError("STOP: never access unrelated or insecure external web origins")
    result={"url":url,"status":"SOURCE_ACCESS_STOP"}
    try:
        req=Request(url,headers={"Accept":"application/json,text/html;q=0.9,*/*;q=0.1",
                                 "User-Agent":"ALAN-Public-University-Data-Access-Audit/2026-11-24"})
        with urlopen(req,timeout=12) as conn:
            endpoint=conn.geturl()
            target=urlparse(endpoint)
            if target.scheme!="https" or target.hostname not in EXPECTED_HOSTS:
                result.update({"status":"PUBLIC_REDIRECT_NOT_TRUSTED_STOP",
                               "redirect_host":target.hostname})
                return result
            mime=conn.headers.get("Content-Type","")
            body=conn.read(MAX_RESPONSE+1)
            if len(body)>MAX_RESPONSE:
                raise ValueError("STOP: over 450kb source metadata")
            result.update({"status":"PUBLIC_METADATA_OR_LANDING_PAGE_ACCESSIBLE",
                           "http_status":conn.status,"content_type":mime,
                           "final_url":endpoint,"bytes_examined":len(body)})
            if "json" in mime:
                metadata=json.loads(body.decode("utf-8"))
                result["dataset_title_hints"]=extract_name_hints(metadata)
                result["public_metadata_not_original_raw_data"]=True
            elif "html" in mime:
                txt=body.decode("utf-8",errors="replace")
                links=CandidateRefs()
                links.feed(txt)
                result["literal_data_link_suffix_candidates"]=links.urls[:15]
                lower=txt.lower()
                if any(flag in lower for flag in (
                    "temporary access restriction",
                    "public access to boris-portal has been temporarily suspended",
                    "access is currently granted exclusively",
                    "university of bern network",
                    "currently restricted to",
                )):
                    result["status"]="PUBLIC_ACCESS_RESTRICTED_HTML_SOFT_200"
                elif "/server/api/" in url:
                    result["status"]="API_RETURNED_HTML_SHELL_NOT_AUTHOR_DATA"
                else:
                    result["status"]="PUBLIC_LANDING_HTML_ONLY_NOT_AUTHOR_DATA"
                result["public_metadata_not_original_raw_data"]=True
            else:
                result["status"]="NON_METADATA_RESPONSE_NOT_ADOPTED"
    except (HTTPError,URLError,TimeoutError,OSError,ValueError) as e:
        result["error"]=f"{type(e).__name__}: {str(e)[:165]}"
    return result

def report():
    attempts=[probe(x) for x in SOURCES]
    return {
        "article_doi":DOI,
        "paper_author_data_record":"https://boris.unibe.ch/151803/",
        "actual_new_trap_data_downloaded":False,
        "independent_lit_dark_pair_n":4,
        "repeated_sampling_nights":6,
        "study_A_paired_trap_observations":48,
        "study_B_no_unlit_control":True,
        "source_access_probes":attempts,
        "status":"PUBLIC_SOURCE_DISCOVERY_ONLY_NOT_TRAP_DATA",
        "university_legacy_record_redirects_to_item":"https://boris-portal.unibe.ch/items/f559ab1e-e693-4fbc-b7fe-683134f6dbe1",
        "source_gate_warning":"HTTP 200 HTML from the API search path is not an authenticated source JSON or actual data; BORIS can return a soft-200 public-access restriction. No access restriction will be bypassed.",
        "new_meta_k_increment":0,
        "next_required":"University dataset actual CSV/XLSX with pair, date, distance, light, order and catch count; hash and full 4-pair design gate before lnRR or RRR.",
        "critical_sampling_warning":"Four spatial pair installations, not 48 independent observations or 99,852 individual assignments.",
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=Path("results/alan_carannante2021_boris_source_gate_20261124.json"))
    args=ap.parse_args()
    x=report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(x,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("CARANNANTE_BORIS_SOURCE_GATE",json.dumps({
        "status":x["status"],"new_k":x["new_meta_k_increment"],
        "attempts":x["source_access_probes"]
    },ensure_ascii=False))

if __name__=="__main__":
    main()
