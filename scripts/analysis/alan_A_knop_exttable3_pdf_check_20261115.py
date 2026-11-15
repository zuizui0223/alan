#!/usr/bin/env python3
"""Knop 2017 original Extended Data Table 3 text-layer independent source check.

Extract ONLY provenance and a handful of interaction coefficients from a
government-hosted bibliography bundle reproducing the Nature paper. The
full third-party PDF or rendered copyrighted pages are never committed.
If the source PDF is inaccessible/textless, mark unavailable rather than
guessing signs from search snippets.
"""
from __future__ import annotations
import hashlib
import json
import pathlib
import re
import urllib.error
import urllib.request

URL="https://planesynormas.mma.gob.cl/archivos/2021/proyectos/ANTECEDENTES_ENTREGADOS_POR_ROC-OIKONOS.pdf"
LIMIT=125_000_000
MAX_HTML_ERRORS=3

def download_pdf()->tuple[dict,bytes]:
    manifest={"url":URL,"http_status":None,"sha256":None,"bytes":None,"error":None}
    try:
        req=urllib.request.Request(URL,headers={"User-Agent":"ALAN-Public-Source-Audit/2026-11-15 (https://github.com/zuizui0223/alan)"})
        with urllib.request.urlopen(req,timeout=100) as response:
            manifest["http_status"]=response.status
            chunks=[]
            total=0
            while True:
                chunk=response.read(1024*1024)
                if not chunk:break
                total+=len(chunk)
                if total>LIMIT:
                    manifest["error"]="SIZE_CAP_125_MB"
                    return manifest,b""
                chunks.append(chunk)
        data=b"".join(chunks)
        if not data.startswith(b"%PDF"):
            manifest["error"]="NOT_PDF_SIGNATURE"
            return manifest,b""
        manifest["bytes"]=len(data)
        manifest["sha256"]=hashlib.sha256(data).hexdigest()
        return manifest,data
    except urllib.error.HTTPError as err:
        manifest["http_status"]=err.code;manifest["error"]="HTTP_"+str(err.code)
    except (urllib.error.URLError,TimeoutError,ValueError) as err:
        manifest["error"]=type(err).__name__+":"+str(err)[:160]
    return manifest,b""


def extract_table_text(pdf_bytes:bytes)->dict:
    import fitz
    doc=fitz.open(stream=pdf_bytes,filetype="pdf")
    found=[]
    for page_index in range(len(doc)):
        text=doc[page_index].get_text(sort=True)
        if "treatLamp" in text or ("baggingOpen" in text and "EST" in text):
            lines=[re.sub(r"\s+"," ",z).strip() for z in text.splitlines()]
            hits=[]
            for i,line in enumerate(lines):
                if re.search(r"treatLamp|baggingOpen|developed fruits|^Intercept$",line,re.I):
                    hits.append((i,line))
            # Avoid redistributing the page, even as text; keep only concise
            # numeric table-row names and a bounded original coefficient subset.
            found.append({"page_index_zero_based":page_index,
                          "first_relevant_lines":[{"line":j,"wording":v[:110]} for j,v in hits[:12]],
                          "full_page_text_preserved":False})
    return {"pdf_pages":len(doc),"pages_matching_terms":len(found),"relevant_page_metadata":found[:3]}


def main():
    manifest,data=download_pdf()
    report={
      "schema":"ALAN_A_Knop2017_ExtendedDataTable3_original_pdf_text_audit_v1",
      "document_kind":"PUBLIC_CHILE_GOVERNMENT_ARCHIVE_REPRINT_OF_NATURE_ORIGINAL",
      "source_article_doi":"10.1038/nature23288",
      "verified_public_copy_download":bool(data),
      "original_study":"Knop et al. Nature 2017, Cirsium oleraceum",
      "source_receipt":manifest,
      "table_extraction":None,
      "known_web_text_outcome":"Earlier web-indexed text reported fruit proportion lighting x bagging beta=-1.40 SE=0.592 p=0.018; fruit-count interaction beta=+0.76 SE=0.309 p=.014, apparently contradicting original figure sign and requiring visual/model factor-level reconciliation.",
      "scale_gate":"Link-scale GLMM coefficients, not fruit log response ratios and not diel-night dependence.",
      "reconstructed_night_dependence_V":None,
      "do_not_trust_raw_sign_without_original_layout":True
    }
    if data:
        try:
            report["table_extraction"]=extract_table_text(data)
        except Exception as err:
            report["extraction_error"]=type(err).__name__+":"+str(err)[:200]
    out=pathlib.Path("results/alan_A_knop2017_exttable3_public_pdf_receipt_20261115.json")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
      "source_available":report["verified_public_copy_download"],
      "download_bytes":manifest["bytes"],"error":manifest["error"],
      "matched_pages":(report["table_extraction"] or {}).get("pages_matching_terms"),
      "extract":report["table_extraction"],
    },ensure_ascii=False))

if __name__=="__main__":
    main()
