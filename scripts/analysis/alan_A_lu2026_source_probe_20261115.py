#!/usr/bin/env python3
"""Audit PUBLIC Lu 2026 reproduction paper source locations without access bypass.

Queries published Crossref/OpenAlex/EuropePMC metadata and standard publisher
landing + supplementary file URLs. Records HTTP status, file hashes, metadata
and whether an authentic DOCX contains source tables; no copyrighted full text
or treatment rows are redistributed.
"""
from __future__ import annotations

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
from html.parser import HTMLParser
from xml.etree import ElementTree

DOI = "10.1093/jpe/rtag160"
URLS = {
    "jpe_author_official": "https://www.jpe.ac.cn/EN/10.1093/jpe/rtag160",
    "oup_publisher": "https://academic.oup.com/jpe/advance-article/doi/10.1093/jpe/rtag160/8731923",
    "oup_supp_advance": "https://academic.oup.com/jpe/advance-article/doi/suppl/10.1093/jpe/rtag160/suppl_file/rtag160_supplementary_data.docx",
    "oup_supp_article": "https://academic.oup.com/jpe/article/doi/suppl/10.1093/jpe/rtag160/suppl_file/rtag160_supplementary_data.docx",
    "oup_accepted_pdf": "https://academic.oup.com/jpe/advance-article-pdf/doi/10.1093/jpe/rtag160/8731923/rtag160.pdf",
    "crossref_metadata": "https://api.crossref.org/works/10.1093/jpe/rtag160",
    "openalex_metadata": "https://api.openalex.org/works/https://doi.org/10.1093/jpe/rtag160",
    "europepmc_metadata": "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:10.1093/jpe/rtag160&format=json",
}
MAX_BYTES = 18 * 1024 * 1024


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
    def handle_starttag(self, tag, attrs):
        if tag != "a":
            return
        d = dict(attrs)
        url = d.get("href", "")
        if re.search(r"rtag160|supplementary|suppl_file|\\.docx|\\.pdf", url, re.I):
            self.links.append(url)


def public_get(url: str) -> tuple[dict, bytes]:
    record = {"url": url, "status": None, "size": None, "sha256": None,
              "content_type": None, "error": None, "final_host": None}
    try:
        req = urllib.request.Request(url, headers={"User-Agent":
            "ALAN-evidence-check/1.0 (scientific source verification; https://github.com/zuizui0223/alan)",
            "Accept": "*/*"})
        with urllib.request.urlopen(req, timeout=25) as resp:
            final_url = resp.geturl()
            record["final_host"] = urllib.parse.urlsplit(final_url).netloc
            record["status"] = resp.status
            record["content_type"] = resp.headers.get("Content-Type", "")
            data = resp.read(MAX_BYTES + 1)
            if len(data) > MAX_BYTES:
                record["error"] = "FILE_OVER_18_MB_LIMIT"
                return record, b""
            record["size"] = len(data)
            record["sha256"] = hashlib.sha256(data).hexdigest()
            return record, data
    except urllib.error.HTTPError as exc:
        record["status"] = exc.code
        record["error"] = "PUBLIC_ACCESS_HTTP_" + str(exc.code)
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        record["error"] = type(exc).__name__ + ": " + str(exc)[:200]
    return record, b""


def inspect_docx(data: bytes) -> dict:
    out = {"is_zip": zipfile.is_zipfile(io.BytesIO(data)), "is_real_docx": False,
           "document_tables": None, "text_paragraphs": None,
           "numeric_cells": None, "treatment_terms_present": []}
    if not out["is_zip"]:
        return out
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        if "word/document.xml" not in z.namelist():
            return out
        out["is_real_docx"] = True
        ns = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
        root = ElementTree.fromstring(z.read("word/document.xml"))
        tables = root.findall(".//" + ns + "tbl")
        paragraphs = root.findall(".//" + ns + "p")
        texts = [node.text or "" for node in root.findall(".//" + ns + "t")]
        joined = " ".join(texts).lower()
        out["document_tables"] = len(tables)
        out["text_paragraphs"] = len(paragraphs)
        out["numeric_cells"] = sum(bool(re.search(r"\\d", " ".join(n.text or "" for n in cell.findall(".//" + ns + "t"))))
                                   for cell in root.findall(".//" + ns + "tc"))
        for keyword in ("seed set", "pedicularis", "euphrasia", "illuminat", "dark", "pollination"):
            if keyword in joined:
                out["treatment_terms_present"].append(keyword)
    return out


def main() -> None:
    results = {}
    discovered_links = []
    for key, url in URLS.items():
        record, data = public_get(url)
        if data and "json" in key:
            try:
                doc = json.loads(data)
                if key == "crossref_metadata":
                    msg = doc.get("message", {})
                    record["publisher"] = msg.get("publisher")
                    record["published"] = str(msg.get("published", {}))[:200]
                    record["publisher_links"] = [
                        {"url": x.get("URL"), "content_type": x.get("content-type")}
                        for x in msg.get("link", [])[:8]]
                elif key == "openalex_metadata":
                    best = doc.get("best_oa_location") or {}
                    record["oa_location"] = {"landing_page_url": best.get("landing_page_url"),
                                              "pdf_url": best.get("pdf_url")}
                    record["open_access"] = doc.get("open_access")
                elif key == "europepmc_metadata":
                    record["result_count"] = doc.get("hitCount")
                    record["pmcid"] = [
                        x.get("pmcid") for x in (doc.get("resultList", {}).get("result") or [])[:5]]
            except (ValueError, TypeError):
                record["metadata_parse_error"] = True
        if data and "text/html" in (record.get("content_type") or "").lower():
            p = LinkParser()
            try:
                p.feed(data.decode("utf-8", "replace"))
                record["article_relevant_links"] = sorted(set(p.links))[:50]
                discovered_links += record["article_relevant_links"]
            except ValueError:
                record["article_relevant_links"] = []
        if data and (data[:2] == b"PK" or "docx" in key):
            record["docx_inspection"] = inspect_docx(data)
        if data and data.startswith(b"%PDF"):
            record["pdf_available_for_separate_methods_review"] = True
        results[key] = record
    report = {
        "schema": "ALAN_A_Lu2026_public_source_probe_v1",
        "date": datetime.now(timezone.utc).isoformat(),
        "doi": DOI,
        "status": "SOURCE_ACCESS_AUDIT_ONLY_NOT_RAW_REANALYSIS",
        "no_credentials_or_bypass": True,
        "checks": results,
        "relevant_links_in_publisher_html": sorted(set(discovered_links))[:60],
        "reproductive_treatment_cells_recovered": None,
        "recomputed_night_x_light_beta": None,
        "interpretation": "Public HTTP errors are access limitations; source rows are not inferred from abstract percentages. A DOCX must be authenticated and table grain audited before any new model."
    }
    pathlib.Path("results/alan_A_lu2026_public_source_probe_20261115.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: {"status": val["status"], "size": val["size"],
                             "type": val["content_type"],
                             "docx": val.get("docx_inspection")}
                      for key, val in results.items()}, ensure_ascii=False))


if __name__ == "__main__":
    main()
