#!/usr/bin/env python3
"""Check ONLY public Wiley-supporting-information routes for Perkin 2014 raw XLSX.

Never bypass access controls, do not copy restricted source bytes into repo.
If source is unavailable, emit DATA_ACCESS_STOP and leave meta k unchanged.
"""
import hashlib
import io
import json
import os
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET
import zipfile

FILE="fwb12426-sup-0004-AppendixS2.xlsx"
URLS=[
    "https://onlinelibrary.wiley.com/doi/suppl/10.1111/fwb.12426/supinfo/"+FILE,
    "https://onlinelibrary.wiley.com/doi/suppl/10.1111/fwb.12426/supinfo/"+FILE+"?download=1",
    "https://onlinelibrary.wiley.com/doi/supinfo/10.1111/fwb.12426",
]
OUT=Path(os.environ.get("PERKIN_GATE_OUTPUT","results/alan_perkin2014_supplement_access_20261124.json"))
MAX_BYTES=8_000_000
def main():
    records=[]
    result=None
    for url in URLS:
        entry={"url":url}
        try:
            req=Request(url,headers={
                "User-Agent":"ALANResearchSourceAccessGate/1.0 (+public DOI 10.1111/fwb.12426)",
                "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,application/octet-stream,*/*"
            })
            with urlopen(req, timeout=10) as obj:
                entry.update({"http_status":obj.status,"final_url":obj.geturl(),"content_type":obj.headers.get("Content-Type","")})
                data=obj.read(MAX_BYTES+1)
            if len(data)>MAX_BYTES:
                raise RuntimeError("content greater than 8MB stop")
            entry["byte_count"]=len(data)
            if not zipfile.is_zipfile(io.BytesIO(data)):
                entry["status"]="NON_XLSX_OR_HTML_STOP"
            else:
                with zipfile.ZipFile(io.BytesIO(data)) as z:
                    if "xl/workbook.xml" not in z.namelist():
                        entry["status"]="ZIP_NOT_XLSX_STOP"
                    else:
                        ns={"x":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
                        top=ET.fromstring(z.read("xl/workbook.xml"))
                        sheets=[x.get("name") for x in top.findall(".//x:sheets/x:sheet",ns)]
                        entry.update({"status":"PUBLIC_XLSX_BYTES_VERIFIED",
                            "sha256":hashlib.sha256(data).hexdigest(),
                            "sheet_names":sheets,
                            "archive_file_names_count":len(z.namelist())})
                        result=entry
            records.append(entry)
            if result:
                break
        except (HTTPError,URLError,TimeoutError,OSError,RuntimeError,ValueError) as ex:
            entry["status"]="DATA_ACCESS_STOP"
            entry["reason"]=f"{type(ex).__name__}: {str(ex)[:220]}"
            records.append(entry)
    response={
        "publication_doi":"10.1111/fwb.12426",
        "source_filename":FILE,
        "status":"PUBLIC_XLSX_BYTES_VERIFIED" if result else "DATA_ACCESS_STOP",
        "attempts":records,
        "raw_bytes_redistributed":False,
        "meta_k_added":0,
        "qualifier":"Metadata/sheet discovery only. Real unit design, row values, sample covariance, and independent meta eligibility still require source-scale audit."
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(response,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":response["status"],"attempts":records},ensure_ascii=False))


if __name__=="__main__":
    main()
