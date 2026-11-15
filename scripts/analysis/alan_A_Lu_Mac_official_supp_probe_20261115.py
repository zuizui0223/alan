#!/usr/bin/env python3
"""Discover a fresh, publicly authorized Lu 2026 supplementary DOCX link.

Distinguish publisher-landing 200 from actual DOCX 200. Do not reuse an
expired silverchair signed CDN URL and do not attempt credential bypass.
Only public URLs are queried and only availability, hashes, and content
structure are stored. If an official DOCX is genuinely downloadable,
record whether Table S4b text exists before any effect is extracted.
Macgregor Wiley Appendix S1 is also checked via publisher's exact
supplement download endpoint, never inferred as an original raw CSV.
"""
from __future__ import annotations
from html.parser import HTMLParser
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree
import zipfile

SOURCES={
 "lu_jpe_publisher":"https://academic.oup.com/jpe/advance-article/doi/10.1093/jpe/rtag160/8731923?login=false",
 "lu_jpe_publisher_alt":"https://academic.oup.com/jpe/advance-article/doi/10.1093/jpe/rtag160/8731923?searchresult=1",
 "lu_jpe_chinese_journal":"https://www.jpe.ac.cn/EN/10.1093/jpe/rtag160",
 "lu_silverchair_static_backfile":"https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/jpe/PAP/10.1093_jpe_rtag160/1/rtag160_supplementary_data.docx",
 "mac_wiley_exact_supp":"https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecs2.2550&file=ecs22550-sup-0001-AppendixS1.pdf",
 "mac_wiley_supp_info":"https://esajournals.onlinelibrary.wiley.com/doi/supinfo/10.1002/ecs2.2550",
}
EXTRA_SUFFIX=("rtag160_supplementary_data.docx","ecs22550-sup-0001-AppendixS1.pdf")
MAX_BYTES=9_000_000
USER_AGENT="ALAN-PublicSupplementaryDiscovery/2026-11-15 (https://github.com/zuizui0223/alan)"

class Links(HTMLParser):
    def __init__(self):
        super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower() not in ("a","link"):
            return
        props=dict(attrs)
        url=props.get("href") or ""
        if url and re.search(r"supplem|s4b|appendix|rtag160|ecs22550|download",url,re.I):
            self.links.append(url)

def get(url:str,limit=MAX_BYTES):
    out={"requested":url[:250],"http_code":None,"content_type":None,
         "final_host":None,"bytes":None,"sha256":None,"error":None}
    data=b""
    try:
        req=urllib.request.Request(url,headers={"User-Agent":USER_AGENT,
            "Accept":"application/vnd.openxmlformats-officedocument.wordprocessingml.document,application/pdf,text/html,*/*"})
        with urllib.request.urlopen(req,timeout=35) as r:
            out["http_code"]=r.status
            out["content_type"]=r.headers.get("Content-Type","")
            out["final_host"]=urllib.parse.urlsplit(r.geturl()).netloc
            data=r.read(limit+1)
        if len(data)>limit:
            out["error"]="max_bytes_exceeded"
            return out,b""
        out["bytes"]=len(data)
        out["sha256"]=hashlib.sha256(data).hexdigest()
    except urllib.error.HTTPError as err:
        out["http_code"]=err.code
        out["error"]="HTTP_"+str(err.code)
    except Exception as err:
        out["error"]=type(err).__name__+": "+str(err)[:170]
    return out,data

def docx_summary(data:bytes)->dict:
    base={"valid_word_docx":False,"table_count":None,
          "has_table_S4b_literal":False,"has_reproductive_seed_text":False,
          "numeric_table_cell_count":None,"raw_plant_rows_recovered":False}
    if not zipfile.is_zipfile(io.BytesIO(data)):
        return base
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        if "word/document.xml" not in z.namelist():return base
        root=ElementTree.fromstring(z.read("word/document.xml"))
        nsp="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
        table_nodes=root.findall(".//"+nsp+"tbl")
        cells=root.findall(".//"+nsp+"tc")
        text=" ".join(x.text or "" for x in root.findall(".//"+nsp+"t"))
        base.update(valid_word_docx=True,table_count=len(table_nodes),
                    has_table_S4b_literal=bool(re.search(r"(Table|Tab[.]?)\s*S\s*4\s*b",text,re.I)),
                    has_reproductive_seed_text=bool(re.search(r"seed|ovule|pollen",text,re.I)),
                    numeric_table_cell_count=sum(bool(re.search(r"\d"," ".join(x.text or "" for x in c.findall(".//"+nsp+"t")))) for c in cells))
    return base

def run():
    out={"schema":"ALAN_A_Lu_Mac_official_2026_fresh_supplement_probe_v1",
         "datetime_utc":datetime.now(timezone.utc).isoformat(),
         "no_bypass_or_authentication":True,
         "probe_type":"AVAILABILITY_AND_STRUCTURE_ONLY",
         "sources":{},
         "discovered_authorized_urls":[],
         "lu_S4b_numeric_treatment_cells_extracted":False,
         "mac_source_seed_CSVs_extracted":False,
         "real_meta_beta_V":None}
    candidates=[]
    for name,url in SOURCES.items():
        entry,data=get(url)
        if data and data.startswith(b"PK"):
            entry["docx_structure"]=docx_summary(data)
        if data and data.startswith(b"%PDF"):
            entry["pdf_bytes_valid"]=True
            # Explicitly no PDF page extraction in this workflow.
        if data and b"<" in data[:100] and "text/html" in (entry["content_type"] or "").lower():
            parser=Links()
            try:
                parser.feed(data.decode("utf-8","replace"))
                resolved=[urllib.parse.urljoin(url,v) for v in parser.links]
                entry["supplement_related_link_count"]=len(resolved)
                candidates.extend(u for u in resolved if
                   any(e.lower() in urllib.parse.urlsplit(u).path.lower() for e in EXTRA_SUFFIX))
            except Exception as exc:
                entry["html_parse_error"]=type(exc).__name__
        out["sources"][name]=entry
    candidates=list(dict.fromkeys(candidates))
    out["discovered_authorized_urls"]=[
      {"host":urllib.parse.urlsplit(u).netloc,"path":urllib.parse.urlsplit(u).path,
       "signed_query_present":bool(urllib.parse.urlsplit(u).query)} for u in candidates[:15]
    ]
    for n,url in enumerate(candidates[:10]):
        # Follow only URLs discovered in the original publisher/official site.
        entry,data=get(url)
        if data.startswith(b"PK"):
            entry["docx_structure"]=docx_summary(data)
        if data.startswith(b"%PDF"):
            entry["pdf_bytes_valid"]=True
        # Avoid storing expiry/signature tokens, private URLs or copyrighted pages.
        entry["requested"]=(urllib.parse.urlsplit(url).scheme+
              "://"+urllib.parse.urlsplit(url).netloc+urllib.parse.urlsplit(url).path)
        out["sources"]["discovered_download_"+str(n)]=entry
    out["downloaded_valid_lu_docx"]=any(
       bool(z.get("docx_structure",{}).get("valid_word_docx"))
       for z in out["sources"].values())
    out["downloaded_valid_mac_supp_pdf"]=any(
       z.get("pdf_bytes_valid") and "mac_" in k for k,z in out["sources"].items())
    dest=Path("results/alan_A_Lu_Mac_fresh_official_supp_gate_20261115.json")
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_http_statuses":{name:z["http_code"] for name,z in out["sources"].items()},
        "discovered_authorized_links":len(out["discovered_authorized_urls"]),
        "lu_docx":out["downloaded_valid_lu_docx"],
        "mac_pdf":out["downloaded_valid_mac_supp_pdf"],
        "genuine_beta_V":out["real_meta_beta_V"]
    },ensure_ascii=False))

if __name__=="__main__":
    run()
