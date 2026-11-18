#!/usr/bin/env python3
"""Original Lu et al 2026 five paired-plot pollination network source gate.

Paper DOI 10.1002/ece3.74368, PMC13645961, original Supporting Information
'Data S1' (ECE3-16-e74368-s001.zip) with author Data.xlsx and R code.
Try only official publication hosts; if source HTTP-403 or ZIP access blocked,
write an explicit source STOP; do not invent plant or pollinator-level effects.
Any later model must treat 5 spatial lit/control blocks as independent units,
not 2160 visitor events, 29 plants, or 4 repeated surveys.

The separate JPE DOI 10.1093/jpe/rtag160 is by the same first authors and
research institution; do not double-count its source records as independent
without actual plot/date identifiers.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import time
import urllib.request
import urllib.error
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

SCHEMA="alan_lu2026_alpine_pollination_real_DataS1_source_gate_v1"
DOI="10.1002/ece3.74368"
SUPP="ECE3-16-e74368-s001.zip"
MAX=15*1024*1024
URLS=[
 "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles",
 "https://pmc.ncbi.nlm.nih.gov/articles/instance/13645961/bin/ECE3-16-e74368-s001.zip",
 "https://onlinelibrary.wiley.com/doi/suppl/10.1002/ece3.74368/supinfo/ece374368-sup-0002-DataS1.zip"
]

def request(url):
    req=urllib.request.Request(url,headers={
      "User-Agent":"ALAN ecology source verification research contact via repository",
      "Accept":"application/zip,application/octet-stream,application/x-zip-compressed,*/*"})
    with urllib.request.urlopen(req,timeout=42) as response:
        final=response.url
        content=response.read(MAX+1)
        return response.status,response.headers.get("Content-Type",""),final,content

def recursively_inventory_archive(archive,depth=0):
    if depth>3:return {"status":"max_nested_depth"}
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        bad=z.testzip()
        if bad:raise ValueError("published ZIP CRC failed")
        listing=[]
        for a in z.infolist():
            if a.is_dir():continue
            row={"filename":a.filename,"bytes":a.file_size,"CRC32":f"{a.CRC:08x}"}
            if a.file_size>MAX:
                row["status"]="too_large_skip"
            elif a.filename.lower().endswith(".zip"):
                row["nested"]=recursively_inventory_archive(z.read(a),depth+1)
            elif a.filename.lower().endswith(".xlsx"):
                x=z.read(a)
                if zipfile.is_zipfile(io.BytesIO(x)):
                    with zipfile.ZipFile(io.BytesIO(x)) as book:
                        nam="xl/workbook.xml"
                        if nam in book.namelist():
                            ns={"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
                            xml=ET.fromstring(book.read(nam))
                            row["sheet_names"]= [
                              el.attrib["name"] for el in xml.findall(".//s:sheet",ns)]
                        row["xlsx_inner_checksum_sha256"]=hashlib.sha256(x).hexdigest()
                else:row["status"]="file_extension_not_true_xlsx"
            elif a.filename.lower().endswith(".r"):
                text=z.read(a)[:750000].decode("utf-8",errors="replace")
                row["r_source_has_original_data_excel"]=("Data.xlsx" in text)
                row["r_source_n_chars"]=len(text)
            listing.append(row)
        return {"type":"true_original_zip","items":listing}

def execute():
    out={
      "schema":SCHEMA,
      "paper_DOI":DOI,
      "PMC_article":"PMC13645961",
      "original_supplementary_file":SUPP,
      "source_paper_original_independent_treatment_control_pairs":5,
      "source_paper_original_observation_year":2024,
      "raw_reference_DOI_has_five_independent_spatial_blocks":True,
      "other_Lu_2026_seed_fitness_JPE_DOI":"10.1093/jpe/rtag160",
      "same_research_group_independence_not_verified":True,
      "source_gates":[],
      "meta_pollinator_replacement_coefficient_fitted":False,
      "species_reproductive_outcome_fitted":False,
    }
    for url in URLS:
        status={"url":url}
        try:
            http,mime,final,payload=request(url)
            status.update(http_status=http,mime=mime,
                          exact_downloaded_bytes=len(payload),
                          "redirected_url_if_different":final if final!=url else None)
            if len(payload)>MAX:
                status["status"]="remote_source_too_large_STOP"
            elif zipfile.is_zipfile(io.BytesIO(payload)):
                status["sha256"]=hashlib.sha256(payload).hexdigest()
                status["contents"]=recursively_inventory_archive(payload)
                status["status"]="original_published_supplementary_zip_verified"
                out["source_gates"].append(status)
                out["status"]="original_published_supplementary_zip_verified"
                break
            else:
                status["status"]="returned_nonzip_source_STOP"
        except urllib.error.HTTPError as e:
            status.update(status="remote_source_access_STOP",http_status=e.code)
        except (urllib.error.URLError,TimeoutError,ValueError,zipfile.BadZipFile) as e:
            status.update(status="source_not_accessible_or_corrupt_STOP",
                          error_type=type(e).__name__)
        out["source_gates"].append(status)
    else:
        out["status"]="source_access_inconclusive_STOP"
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    d=execute()
    o=Path(a.out);o.parent.mkdir(parents=True,exist_ok=True)
    o.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("ALAN_LU2026_ORIGINAL_POLLINATION_SOURCE_GATE",
          json.dumps(d,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
