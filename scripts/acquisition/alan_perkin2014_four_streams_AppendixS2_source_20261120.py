#!/usr/bin/env python3
"""Perkin et al. 2014 4-stream independent field experiment, genuine Appendix S2.

Original peer-reviewed paper 10.1111/fwb.12426 lists Appendix S2 XLSX
(fwb12426-sup-0004-AppendixS2.xlsx): "Original data from all reaches".
Try exact official Wiley publisher supplement paths. Publisher GET may
fail 403; if so SOURCE_STOP and no new effect/variance claim. Do not
replace source with an unverified mirror or manually invent source cells.

If source downloaded, output original content SHA, worksheet names,
row counts, column-header candidates and number of unique coded stream
IDs if obvious, not raw plot IDs or biological observation rows. True
4 independently paired streams can only be fitted after checking design,
time, experimental arm and residual/error covariance.
"""
import argparse
import hashlib
import io
import json
import urllib.error
import urllib.request
from pathlib import Path

DOI="10.1111/fwb.12426"
APPENDIX_NAME="fwb12426-sup-0004-AppendixS2.xlsx"
WILEY_PUBLISHER_CANDIDATES=[
  "https://onlinelibrary.wiley.com/doi/suppl/10.1111/fwb.12426/supinfo/"+APPENDIX_NAME,
  "https://onlinelibrary.wiley.com/doi/supinfo/10.1111/fwb.12426/"+APPENDIX_NAME,
  "https://onlinelibrary.wiley.com/doi/suppl/10.1111/fwb.12426/supinfo/"+APPENDIX_NAME+"?download=1",
]
MAX_BYTES=1_500_000

def pull(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"ALAN-four-original-paired-streams-source-audit/2026-11-20",
        "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"})
    with urllib.request.urlopen(req,timeout=25) as r:
        data=r.read(MAX_BYTES+1)
        if len(data)>MAX_BYTES:raise ValueError("Source original supplement file too large")
        return data

def workbook_source_shape(b):
    from openpyxl import load_workbook
    if not b.startswith(b"PK\x03\x04"):
        raise ValueError("Source response is not XLSX")
    book=load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    tabs=[]
    for sheet in book.worksheets:
        count=0
        first_headers=[]
        candidate_fields=set()
        for row in sheet.iter_rows(values_only=True):
            cells=[str(v).strip() if v is not None else "" for v in row]
            if not any(cells):continue
            count+=1
            if count<=4:
                first_headers.append({"nonempty_cells":sum(bool(v) for v in cells),
                                      "candidate_header_labels":[c for c in cells if c and
                                          (c.lower().replace(" ","").replace("_","") in
                                          {"river","stream","site","reach","treatment","light","night",
                                           "drift","emergence","density","year","sample","block","time"})]})
            if count>50000:
                raise ValueError("Unexpected original XLSX row volume")
        tabs.append({"original_sheet":sheet.title,"original_nonempty_rows_including_headers":count,
                     "early_header_shapes_only":first_headers,
                     "original_numeric_cells_not_exported":True,
                     "original_stream_or_reach_ID_values_not_exported":True})
    book.close()
    return tabs

def audit(fetch=pull):
    result={
        "schema":"alan_perkin2014_original_four_paired_streams_AppendixS2_XLSX_source_gate_v1",
        "original_publication_DOI":DOI,
        "original_publisher_supporting_filename":APPENDIX_NAME,
        "original_independent_stream_pair_count":4,
        "supporting_information_announced_by_publisher":True,
        "source_status":"SOURCE_STOP",
        "official_publisher_source_attempts":[],
        "source_original_workbook_sha256":None,
        "true_source_independent_stream_means_or_covariance_verified":False,
        "night_drift_original_field_meta_effect_fitted":False,
        "aquatic_emergence_original_field_meta_effect_fitted":False,
        "no_raw_site_or_stream_observations_exported":True,
    }
    for url in WILEY_PUBLISHER_CANDIDATES:
        try:
            b=fetch(url)
            sheets=workbook_source_shape(b)
            if len(sheets)<1:raise ValueError("Original workbook contains no worksheets")
            result["source_status"]="ORIGINAL_WILEY_APPX_S2_BYTES_AND_SHEETS_VERIFIED"
            result["source_original_workbook_sha256"]=hashlib.sha256(b).hexdigest()
            result["source_original_workbook_bytes"]=len(b)
            result["source_original_worksheet_structure_only"]=sheets
            result["official_publisher_source_success_url"]=url
            result["official_publisher_source_attempts"].append({"source":"Wiley publisher","status":"ORIGINAL_BYTES_VERIFIED"})
            break
        except (ValueError,ImportError,urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            result["official_publisher_source_attempts"].append({
                "source_url":url,"status":"SOURCE_STOP",
                "error_type":type(e).__name__,
                **({"http_status":e.code} if isinstance(e,urllib.error.HTTPError) else {})
            })
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=audit()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("PERKIN2014_FOUR_TRUE_STREAMS_SOURCE_STATUS",result["source_status"])
    print("PERKIN2014_EXACT_PUBLISHER_XLSX_ATTEMPTS",result["official_publisher_source_attempts"])
    print("PERKIN2014_SOURCE_WORKBOOK_STRUCTURE",result.get("source_original_worksheet_structure_only",[]))

if __name__=="__main__":
    main()
