#!/usr/bin/env python3
"""Authentic Xiao et al. 2025 official Elsevier Supplement DOCX: numeric table audit.

Published direct manipulated ALAN×animal-exclusion 40-plot field study,
DOI 10.1016/j.biocon.2025.111472; archived officially as Elsevier
S0006320725005099-mmc1.docx. No source DOCX bytes or paragraphs published;
retain compact paraphrased factor labels, selected model table cell numbers,
and original file SHA, not full copyrighted appendix text.

CRITICAL: animal exclusion != pollinator exclusion or day/night exclusion.
Fruiting inflorescences per plant != newly matured fruit per original flower.
"""

from __future__ import annotations
import hashlib
import io
import json
import re
import urllib.request
import zipfile
from pathlib import Path
from xml.etree import ElementTree

URL="https://ars.els-cdn.com/content/image/1-s2.0-S0006320725005099-mmc1.docx"
SHA256="127950cecd9ac78c0927b89df06b307e92a7a9445b0c054aa40a2e9d272e8a9f"
NS="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
MAX_BYTES=3_000_000
KEYWORD=re.compile(r"\b(?:alan|artificial light|light pollution|night|control|cage|animal|exclu\w*|fruit\w*|seed\w*|flower\w*|inflorescen\w*|interact\w*|growth|biomass|develop\w*|pheno\w*)",re.I)
NUMERIC=re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")

def get_authentic()->bytes:
    req=urllib.request.Request(URL,headers={
       "User-Agent":"ALAN-Official-Publisher-Supplement-Science-Audit/1.0 (https://github.com/zuizui0223/alan)"
    })
    with urllib.request.urlopen(req,timeout=35) as response:
        data=response.read(MAX_BYTES+1)
    if len(data)>MAX_BYTES:raise RuntimeError("Above source size bound")
    digest=hashlib.sha256(data).hexdigest()
    if digest!=SHA256:raise RuntimeError("Original publication supplement hash drift "+digest)
    if not zipfile.is_zipfile(io.BytesIO(data)):raise RuntimeError("Not a valid DOCX zip")
    return data

def text(paragraph)->str:
    return " ".join("".join(n.text or "" for n in r.findall(".//"+NS+"t")) for r in paragraph.findall(".//"+NS+"r")).strip()

def paragraph_text(node)->str:
    return "".join(n.text or "" for n in node.findall(".//"+NS+"t")).strip()

def extract(data:bytes)->dict:
    with zipfile.ZipFile(io.BytesIO(data)) as zipobj:
        root=ElementTree.fromstring(zipobj.read("word/document.xml"))
        embedded_names=[n for n in zipobj.namelist() if n.startswith("word/media/")]
    body=root.find(".//"+NS+"body")
    if body is None:raise ValueError("No DOCX body")
    paragraphs=[]
    tables=[]
    for child in body:
        if child.tag==NS+"p":
            value=paragraph_text(child)
            if value:paragraphs.append(value)
        if child.tag==NS+"tbl":
            ix=len(tables)
            rows=[]
            for tr in child.findall("./"+NS+"tr"):
                values=[paragraph_text(cell).strip() for cell in tr.findall("./"+NS+"tc")]
                if any(values):rows.append(values)
            caption=next((p for p in reversed(paragraphs[-5:]) if re.search(r"table|supplement|effect|model|fruit|seed",p,re.I)),"")
            selected=[]
            headers=[v[:72] for v in rows[0][:10]] if rows else []
            for i,r in enumerate(rows):
                brief=[re.sub(r"\s+"," ",v)[:90] for v in r[:10]]
                if any(KEYWORD.search(v) for v in r[:2]):
                    selected.append({"row_index":i,"first_factor_values":brief[:2],
                                     "numeric_cell_tokens":[NUMERIC.findall(v)[:4] for v in r[2:10]],
                                     "source_numeric_row_not_model_verified":True})
            # Some original appendices store generic model rows, e.g.
            # 'Light treatment' or 'Animal treatment', in first column.
            tables.append({
                "table_index":ix,"row_count":len(rows),
                "max_column_count":max((len(r) for r in rows),default=0),
                "header_cells_bounded":headers,
                "caption_context_bounded":caption[:160],
                "selected_factor_rows_bounded":selected[:35],
                "n_factor_rows_selected":len(selected),
                "is_verified_reproductive_response_from_this_table":False
            })
    selected_paragraph_keywords=[]
    for i,txt in enumerate(paragraphs):
        if KEYWORD.search(txt) and (
            re.search(r"table|Figure S\d|sampling|plot|fruiting|inflorescence|animal exclusion|interaction",txt,re.I)
        ):
            selected_paragraph_keywords.append({
                "paragraph_index":i,
                "contains_fruit":bool(re.search(r"fruit\w*|seed\w*",txt,re.I)),
                "contains_plot":bool(re.search("plot|exclusion",txt,re.I)),
                "prefix_for_human_locator":re.sub(r"\s+"," ",txt)[:90],
            })
    return {
        "schema":"ALAN_Xiao2025_official_Elsevier_supplement_structural_numeric_audit_v1",
        "original_article_doi":"10.1016/j.biocon.2025.111472",
        "official_supplement_url":URL,
        "official_supplement_SHA256":SHA256,
        "official_supplement_bytes":len(data),
        "document_paragraphs":len(paragraphs),
        "document_tables":len(tables),
        "embedded_media_items":len(embedded_names),
        "n_keyword_paragraphs":len(selected_paragraph_keywords),
        "keyword_locators":selected_paragraph_keywords[:25],
        "table_shapes_and_bounded_source_numerics":tables,
        "paper_design":"40 original 4m2 plots, light vs dark × animal exclusion/open. Animal exclusion's protected taxa remain unverified.",
        "source_outcome":"number of fruiting inflorescences; NOT per original flower mature-fruit success",
        "source_is_mature_fruit_yield_experiment":False,
        "source_published_table_factor_beta_SE_reconstructed":False,
        "independent_light_cluster_se_estimated":False,
        "nighttime_insect_dependence_estimated":False,
        "meta_beta_V":None,
        "limitations":"This parses original supplemental tables; table layout, coefficients, estimand and animal exclusion meaning require content reading before numerical promotion."
    }

def main():
    report=extract(get_authentic())
    dest=Path("results/alan_Xiao2025_true_supplement_model_source_audit_20261115.json")
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "verified_SHA256":report["official_supplement_SHA256"],
        "paragraphs":report["document_paragraphs"],
        "tables":report["document_tables"],
        "media":report["embedded_media_items"],
        "table_receipts":report["table_shapes_and_bounded_source_numerics"],
        "independent_y_SE":report["independent_light_cluster_se_estimated"]
    },ensure_ascii=False))

if __name__=="__main__":
    main()
