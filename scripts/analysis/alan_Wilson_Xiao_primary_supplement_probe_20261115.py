#!/usr/bin/env python3
"""Wilson 2021 yucca and Xiao 2025 ALAN flower original-source acquisition gate.

Only public publisher/official institutional and article DOI metadata endpoints.
Inspect file signatures, archive member lists, byte hashes and URL availability.
Do NOT report a literature DOI or nominal supplementary ZIP as raw treatment
cells unless archive actually opens and contains source fruit/plant records.

The Xiao paper has a 40-plot ALAN x animal-exclusion field design, measures
FRUITING INFLORESCENCES and phenology; it is NOT a proven 40-replicate
day-v-night insect pollination exclusion. We never relabel mammal/herbivore/
other animal exclusion as a diel pollinator-access gate without source methods.
"""

from __future__ import annotations
from html.parser import HTMLParser
import hashlib
import io
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone

PAPERS={
    "wilson2021":{
        "doi":"10.1093/icb/icab010",
        "endpoints":{
            "oup_article":"https://academic.oup.com/icb/article/61/3/1122/6173993",
            "oup_article_abstract":"https://academic.oup.com/icb/article-abstract/61/3/1122/6173993",
            "oup_supp_zip_page":"https://academic.oup.com/icb/article/61/3/1122/6173993#supplementary-data",
            "oup_supp_zip_guess_1":"https://academic.oup.com/icb/article/61/3/1122/6173993/suppl_file/icab010_supplementary_data.zip",
            "oup_supp_zip_guess_2":"https://academic.oup.com/icb/article/61/3/1122/6173993/suppl_file/icab010_supplementary_data.zip?download=1",
            "oup_supp_zip_guess_3":"https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/icb/61/3/10.1093_icb_icab010/1/icab010_supplementary_data.zip",
            "calpoly_original_thesis":"https://digitalcommons.calpoly.edu/theses/2224/",
            "crossref_original":"https://api.crossref.org/works/10.1093/icb/icab010",
        },
        "correct_unit":"5 field sites, 15 randomly lit/dark transects (8/7), 58 fruit-recruitment plants",
        "plant_outcome":"mature fruit count per plant, not site-independent binary fruit set",
    },
    "xiao2025":{
        "doi":"10.1016/j.biocon.2025.111472",
        "endpoints":{
            "elsevier_official":"https://www.sciencedirect.com/science/article/pii/S0006320725005099",
            "elsevier_supp_pdf":"https://ars.els-cdn.com/content/image/1-s2.0-S0006320725005099-mmc1.pdf",
            "elsevier_supp_docx":"https://ars.els-cdn.com/content/image/1-s2.0-S0006320725005099-mmc1.docx",
            "elsevier_supp_xlsx":"https://ars.els-cdn.com/content/image/1-s2.0-S0006320725005099-mmc1.xlsx",
            "crossref_original":"https://api.crossref.org/works/10.1016/j.biocon.2025.111472",
            "openalex_original":"https://api.openalex.org/works/https://doi.org/10.1016/j.biocon.2025.111472",
        },
        "correct_unit":"40 original field 4-square-meter plots; ALAN vs dark crossed with animal-exclusion vs open. Plot allocation/animal-exclusion definition/individual reproduction cells to inspect.",
        "plant_outcome":"fruit-bearing inflorescence count per plant; not initially tagged flower binary fruit set",
    }
}
CAP=15_000_000
UA="ALAN-Publisher-Supplement-Provenance/2026-11-15 (+https://github.com/zuizui0223/alan)"

class RelevantLinks(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag.lower() not in ("a","link","meta"):return
        val=d.get("href") or d.get("content") or ""
        if val and re.search(r"icab010|supplem|mmc[0-9]|supplement|download|figshare|data[._\-]",val,re.I):
            self.links.append(val[:600])


def request(url):
    info={"host":urllib.parse.urlsplit(url).hostname,"path":urllib.parse.urlsplit(url).path,
          "status":None,"content_type":None,"bytes":None,"sha256":None,"error":None}
    try:
        req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"*/*"})
        with urllib.request.urlopen(req,timeout=35) as res:
            info["status"]=res.status
            info["final_host"]=urllib.parse.urlsplit(res.geturl()).hostname
            info["content_type"]=res.headers.get("Content-Type","")
            data=res.read(CAP+1)
        if len(data)>CAP:
            info["error"]="OVER_15MB_CAP";return info,b""
        info["bytes"]=len(data)
        info["sha256"]=hashlib.sha256(data).hexdigest()
        return info,data
    except urllib.error.HTTPError as e:
        info["status"]=e.code
        info["error"]="HTTP_"+str(e.code)
    except (urllib.error.URLError,ValueError,TimeoutError) as e:
        info["error"]=type(e).__name__+":"+str(e)[:140]
    return info,b""


def inspect(data,info):
    if data[:2]==b"PK" and zipfile.is_zipfile(io.BytesIO(data)):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            members=[]
            for member in z.infolist()[:200]:
                members.append({"name":member.filename[:200],"bytes":member.file_size})
            info["zip_valid"]=True
            info["zip_members"]=members
            info["plausible_experimental_csv"]=any(
                m["name"].lower().endswith((".csv",".tsv",".txt",".xlsx")) for m in members)
            info["genuine_reproduction_cell_rows_extracted"]=False
    elif data[:4]==b"%PDF":
        info["pdf_valid"]=True
        info["genuine_reproduction_cell_rows_extracted"]=False
    elif data.lstrip().startswith((b"{",b"[")):
        try:
            j=json.loads(data)
            if isinstance(j,dict):
                info["top_json_keys"]=list(j)[:20]
                if "message" in j and isinstance(j["message"],dict):
                    info["publisher_metadata_links"]=j["message"].get("link",[])[:8]
                if "best_oa_location" in j:
                    oa=j.get("best_oa_location") or {}
                    info["best_oa_location"]={
                        "landing_page_url":oa.get("landing_page_url"),
                        "pdf_url":oa.get("pdf_url")}
        except ValueError:info["json_parse_error"]=True
    elif "text/html" in (info["content_type"] or "") or b"<html" in data[:200].lower():
        parser=RelevantLinks()
        parser.feed(data.decode("utf-8","replace"))
        info["relevant_public_html_links"]=sorted(set(parser.links))[:50]


def run():
    out={"schema":"ALAN_Wilson2021_Xiao2025_public_reproductive_supplement_probe_v1",
         "created_utc":datetime.now(timezone.utc).isoformat(),
         "status":"SOURCE_ACCESS_CHECK_NOT_EFFECT_EXTRACTION",
         "source_archives_raw_downloads_not_redistributed":True,
         "papers":{},
         "new_site_cluster_valid_plant_logRR_and_SE":0,
         "new_valid_same_population_dark_V":0,
         "new_main_A_beta_V":None}
    for name,paper in PAPERS.items():
        checks={}
        for key,url in paper["endpoints"].items():
            info,data=request(url)
            if data:inspect(data,info)
            checks[key]=info
        out["papers"][name]={
            "doi":paper["doi"],
            "preexisting_original_experimental_unit":paper["correct_unit"],
            "reproductive_measure":paper["plant_outcome"],
            "checks":checks,
            "genuinely_recovered_randomized_treatment_cell_table":False,
            "reproduction_lit_dark_logRR_with_site_SE":None,
            "nighttime_insect_dependence_measured":False
        }
    target=pathlib.Path("results/alan_Wilson_Xiao_official_supplement_source_probe_20261115.json")
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "http":{paper:{k:v["status"] for k,v in d["checks"].items()} for paper,d in out["papers"].items()},
        "valid_zips":{paper:[k for k,v in d["checks"].items() if v.get("zip_valid")] for paper,d in out["papers"].items()},
        "effect_sizes_acquired":out["new_site_cluster_valid_plant_logRR_and_SE"],
        "beta_V":out["new_main_A_beta_V"],
    },ensure_ascii=False))


if __name__=="__main__":
    run()
