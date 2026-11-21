#!/usr/bin/env python3
"""Haas-Desmarais et al. 2025/2026 herbivory-pollination 171 studies SUPPLEMENT source gate.

Paper 10.1093/aob/mcaf258 / PMC13095882, original supplement ZIP.
Goal: find an actually machine-readable 1,348-case dataset, including original
study ID, plant binomial, herbivore damage type, pollination/fitness outcome,
Hedges g and sampling SE/variance; before any wetness moderator is estimated.

Do not treat supplementary PDF descriptions as 1,348 accessible rows.
No causal ALAN effect in this meta. Stop if raw rows cannot be source authenticated.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,urllib.error,urllib.request,zipfile
from pathlib import Path

SCHEMA="alan_herbivory_pollinator_171_original_source_gate_20261120_v1"
DOI="10.1093/aob/mcaf258"
URLS=[
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13095882/supplementaryFiles",
  "https://pmc.ncbi.nlm.nih.gov/articles/instance/13095882/bin/mcaf258_supplementary_data.zip",
]
MAX_COMPRESSED=14<<20
MAX_MEMBER=4<<20

def retrieve(url):
    req=urllib.request.Request(url,headers={
      "User-Agent":"ecological-meta-source-audit alan 2026",
      "Accept":"application/zip,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=55) as res:
        blob=res.read(MAX_COMPRESSED+1)
        if len(blob)>MAX_COMPRESSED:
            raise ValueError("SUPPLEMENT_ZIP_OVER_SIZE")
        if not zipfile.is_zipfile(io.BytesIO(blob)):
            raise ValueError("NOT_ZIP_CONTENT")
        return blob,res.geturl()

def inspect_docx(blob):
    """Inspect only archived OOXML table structure/section titles; not raw data."""
    from xml.etree import ElementTree as ET
    W="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        if z.testzip() is not None:raise ValueError("BAD_DOCX_ARCHIVE_CRC")
        root=ET.fromstring(z.read("word/document.xml"))
        body=root.find(f".//{W}body")
        if body is None:raise RuntimeError("NO_WORD_BODY")
        def plain(elem):
            return "".join(t.text or "" for t in elem.iter(f"{W}t")).strip()
        tables=body.findall(f".//{W}tbl")
        table_info=[]
        for t in tables:
            rows=t.findall(f"./{W}tr")
            cells=lambda tr:[plain(v)[:110] for v in tr.findall(f"./{W}tc")]
            table_info.append({
              "row_count":len(rows),
              "first_row_cell_headings":cells(rows[0])[:20] if rows else [],
              "second_row_cells":cells(rows[1])[:8] if len(rows)>1 else []
            })
        paragraphs=[plain(p) for p in body.findall(f"./{W}p")]
        headings=[p[:160] for p in paragraphs if p and (
                  p.startswith("Appendix") or p.startswith("Table") or
                  p.startswith("Fig") or "data" in p.lower())][:30]
        return {"type":"source_docx","paragraphs":len(paragraphs),
            "tables":len(tables),"tables_row_counts":[t["row_count"] for t in table_info],
            "tables_headers_preview":table_info[:12],
            "titles_data_references":headings,
            "images_or_drawings_are_not_effect_size_rows":True,
            "study_level_meta_raw_case_data_not_inferred_from_manuscript_text":True}

def inspect_zip(blob,max_depth=1,depth=0):
    out=[]
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        if z.testzip() is not None:raise ValueError("INVALID_ZIP_CRC")
        for info in z.infolist():
            if info.is_dir():continue
            label=info.filename
            item={"name":label,"bytes":info.file_size,"crc32":f"{info.CRC:08x}"}
            typ=Path(label).suffix.lower()
            item["extension"]=typ
            if info.file_size>MAX_MEMBER:
                item["status"]="EXCESSIVE_MEMBER_SIZE_NO_READ"
            elif typ in (".csv",".tsv",".txt"):
                b=z.read(info)
                item["sha256"]=hashlib.sha256(b).hexdigest()
                try:
                    head=b[:300000].decode("utf-8-sig")
                    for sep in (",","\t",";"):
                        reader=csv.DictReader(io.StringIO(head),delimiter=sep)
                        names=reader.fieldnames or []
                        if len(names)>3:
                            item["status"]="MACHINE_READABLE_DELIMITED"
                            item["columns"]=names[:120]
                            item["delimiter"]=sep
                            item["source_row_count_estimate"]=b.count(b"\n")-1
                            break
                    if "status" not in item:
                        item["status"]="NO_WELL_FORMED_DELIMITED_HEADER"
                except UnicodeError:
                    item["status"]="TEXT_ENCODING_NOT_UTF8"
            elif typ==".docx":
                b=z.read(info)
                item["sha256"]=hashlib.sha256(b).hexdigest()
                item["status"]="ORIGINAL_WORD_SUPPLEMENT_STRUCTURE"
                item["layout"]=inspect_docx(b)
            elif typ==".zip" and depth<max_depth:
                inner=z.read(info)
                if zipfile.is_zipfile(io.BytesIO(inner)):
                    item["status"]="NESTED_ZIP_VERIFIED"
                    item["nested"]=inspect_zip(inner,max_depth,depth+1)
                else:item["status"]="NESTED_ZIP_INVALID"
            else:
                item["status"]="ORIGINAL_BINARY_DOCUMENT_NOT_A_CSV"
            out.append(item)
    return out

def run():
    out={"schema":SCHEMA,"article_doi":DOI,
         "original_supplement_identifier":"mcaf258_supplementary_data.zip",
         "source_reported_original_studies":171,
         "source_reported_aggregate_effect_cells":1348,
         "source_has_ALAN_manipulation":False,
         "paper_meta_effect_not_an_ALAN_wet_species_moderator":True,
         "source_original_raw_effect_cases_loaded":False,
         "source_urls_tried":[]}
    for url in URLS:
        try:
            blob,final=retrieve(url)
            out["zip_SHA256"]=hashlib.sha256(blob).hexdigest()
            out["zip_bytes"]=len(blob)
            out["official_source_download_url"]=url
            out["resolved_url"]=final
            out["source_members"]=inspect_zip(blob)
            out["original_zip_CRC_checked"]=True
            # Structural audit, not a full proof of 1348 independent effects.
            discovered=[]
            def flatten(entries):
                for x in entries:
                    if x["status"]=="MACHINE_READABLE_DELIMITED":
                        discovered.append(x)
                    flatten(x.get("nested",[]))
            flatten(out["source_members"])
            out["machine_readable_candidate_files"]=len(discovered)
            out["status"]="PUBLISHED_ARCHIVE_SOURCE_VERIFIED"
            out["plant_moisture_by_ALAN_moderator_fitted"]=False
            out["plant_moisture_by_herbivore_moderator_fitted"]=False
            return out
        except Exception as e:
            out["source_urls_tried"].append({
                "url":url,"error_type":type(e).__name__,
                "http_code":getattr(e,"code",None),
                "details":str(e)[:200] if not hasattr(e,"code") else None})
    out["status"]="OFFICIAL_SUPPLEMENT_ACCESS_STOP"
    out["plant_moisture_by_ALAN_moderator_fitted"]=False
    out["plant_moisture_by_herbivore_moderator_fitted"]=False
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    j=run()
    p=Path(args.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(j,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("HERBIVORY171_ORIGINAL_META_SUPPLEMENT_SOURCE_GATE "+json.dumps(j,sort_keys=True),flush=True)

if __name__=="__main__":main()
