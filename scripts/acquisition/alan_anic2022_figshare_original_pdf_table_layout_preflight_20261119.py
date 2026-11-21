#!/usr/bin/env python3
"""Authentic Anic 2022 Figshare 24-page original data PDF table preflight.

This is FIRST-PASS table identification and extraction feasibility, NOT a fit.
The published Figshare attachment is a PDF (~147 kB, 24 pages) containing
numeric-heavy rows. Extract an original file SHA, original page lengths,
treatment/organism header tokens and line SHAPES only. No plant biomass,
identifier, plot coordinates, species×plot measurements or fabricated model
covariances are committed.

PDF->table is fragile: numeric body text can wrap and columns align by
position, not by ASCII space. Only fit a new original-data GLM after verified
source table headers, unique physical trough IDs and sampled independent
light cell counts are confirmed. Never OCR PDFs if text exists.
"""
import argparse
import hashlib
import io
import json
import re
import urllib.request
from pathlib import Path

FIGSHARE_API="https://api.figshare.com/v2/articles/20304039"
FILE_ID=36265164
FILE_SHA256="baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0"
DOI="10.1002/ece3.9157"
ALLOWED_WORDS=(
    "treatment","light","trophic","mesocosm","trough","species",
    "biomass","shoot","root","carbon","nitrogen",
    "lotus","holcus","plant","dry","replicate","data",
)
MAX_DOWNLOAD_BYTES=1_000_000


def get(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-original-Anic2022-source-PDF-reading/2026-11-19",
        "Accept":"application/pdf,application/json,*/*"})
    with urllib.request.urlopen(req,timeout=35) as f:
        b=f.read(MAX_DOWNLOAD_BYTES+1)
        if len(b)>MAX_DOWNLOAD_BYTES:raise ValueError("Original data file oversized")
        return b


def masking(line):
    """Keep only structural words and numeric token SHAPES, not raw values."""
    # Remove real experimental IDs, punctuation contents and numeric cells.
    s=re.sub(r"(?<![A-Za-z])[+-]?\d+(?:[.,]\d+)?", "<NUM>",str(line))
    s=re.sub(r"\b(?:[A-Za-z]*\d+[A-Za-z0-9]*)\b","<CODE>",s)
    return s.strip()[:180]


def study_pdf_preflight(original,strict=True):
    from pypdf import PdfReader
    if hashlib.sha256(original).hexdigest()!=FILE_SHA256:
        raise ValueError("Original Figshare archived PDF SHA mismatch")
    if not original.startswith(b"%PDF-"):
        raise ValueError("Not actually the author original data PDF")
    reader=PdfReader(io.BytesIO(original))
    if strict and len(reader.pages)!=24:
        raise ValueError("Original Figshare 24-page source has drifted")
    report=[]
    total_numeric=0
    for i,page in enumerate(reader.pages):
        text=page.extract_text() or ""
        text_layout=page.extract_text(extraction_mode="layout") or ""
        if not text.strip() and strict:
            raise ValueError("Original data PDF page has no selectable text; do not OCR")
        lines=[x.rstrip() for x in text.splitlines() if x.strip()]
        lay=[x.rstrip() for x in text_layout.splitlines() if x.strip()]
        numeric_heavy=sum(len(re.findall(r"[+-]?\d+(?:\.\d+)?",line))>=4 for line in lines)
        total_numeric+=numeric_heavy
        word_match=[w for w in ALLOWED_WORDS if re.search(r"\b"+w+r"\b",text,re.I)]
        # Include the first eight MASKED lines to identify separated header
        # labels. No original numeric content, trough IDs or plot rows.
        selected_lines=[masking(x) for x in lines[:8]]
        page_report={
          "original_page_zero_index":i,
          "text_characters":len(text),
          "nonempty_standard_text_lines":len(lines),
          "nonempty_layout_lines":len(lay),
          "lines_with_at_least_four_numeric_tokens":numeric_heavy,
          "source_keyword_hits":word_match,
          "first_eight_original_line_shapes_masked":selected_lines,
          "numeric_measurements_exported":False,
          "original_plot_ID_exported":False,
        }
        report.append(page_report)
    return {
        "schema":"alan_original_anic2022_figshare_24page_numeric_PDF_table_preflight_v1",
        "source_primary_DOI":DOI,
        "source_figshare_article_id":20304039,
        "source_original_file_id":FILE_ID,
        "original_file_SHA256":FILE_SHA256,
        "original_file_bytes":len(original),
        "original_24page_data_document_pages":len(report),
        "original_pdf_total_numeric_dense_lines":total_numeric,
        "page_shape_only":report,
        "source_status":"AUTHENTIC_PDF_BYTES_AND_SELECTABLE_TEXT_VERIFIED",
        "plant_species_x_physical_mesocosm_rows_extracted":False,
        "biological_values_transcribed":False,
        "author_original_54_trough_light_trophic_assignment_verified":False,
        "new_original_per_species_light_x_trophic_effect_fitted":False,
        "PDF_source_data_not_independent_light_experiment":True,
        "source_reference_headers_require_layout_validation_before_model":True
    }


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    meta=json.loads(get(FIGSHARE_API))
    files=[f for f in meta.get("files",[]) if f.get("id")==FILE_ID]
    if len(files)!=1:raise ValueError("Exact published Figshare file ID not found")
    raw=get(files[0]["download_url"])
    d=study_pdf_preflight(raw,strict=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n")
    print("ANIC2022_PDF_SOURCE_SHA",d["original_file_SHA256"])
    print("ANIC2022_PDF_NUMERIC_SOURCE_LINES",d["original_pdf_total_numeric_dense_lines"])
    print("ANIC2022_PDF_PAGE_FIRST_LINE_SHAPES",[
      {"page":p["original_page_zero_index"],
       "rows":p["nonempty_standard_text_lines"],
       "shapes":p["first_eight_original_line_shapes_masked"],
       "terms":p["source_keyword_hits"]}
      for p in d["page_shape_only"] if p["original_page_zero_index"] in (0,1,12,13,14,23)
    ])
    print("ANIC2022_PDF_ORIGINAL_REANALYSIS_GATE",d["new_original_per_species_light_x_trophic_effect_fitted"])


if __name__=="__main__":
    main()
