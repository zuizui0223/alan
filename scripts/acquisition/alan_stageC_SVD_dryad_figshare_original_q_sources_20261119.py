#!/usr/bin/env python3
"""Stage C: original single visit pollen delivery primary dataset source gate.

Official Dryad 10.5061/dryad.05qfttf1n: Földesi et al. 2020,
28 papers, 127 honey-bee-referenced SVD records: 94 bees, 33 flies.
Original source contains NO nocturnal moth SVD comparison.

Official Figshare 10.6084/m9.figshare.26502400: Ollerton et al. 2025,
Knautia arvensis single-visit butterfly/hoverfly/bumblebee,
NOT nocturnal moths. Same flower but wrong Lepidoptera guild to
treat Lu 2026 nocturnal moths as effective substitute data.

Acquire exact official metadata and original source bytes if available.
Print/store only file names, checksums, table columns and category cardinality
without original biological observations/plot or visitor ID rows.
No imputation of q, no qfly/qmoth parameter or compensated plant fraction.
"""
import argparse
import csv
import hashlib
import io
import json
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

DRYAD_DOI="10.5061/dryad.05qfttf1n"
FIGSHARE_ID=26502400
DOI_ORIGINAL_META="10.1111/1365-2664.13798"
DOI_ORIGINAL_PRIMARY="10.1111/jen.13345"
MAX_BYTES=6_000_000


def get(url):
    req=urllib.request.Request(url,headers={
      "User-Agent":"ALAN-stage-C-pollen-efficiency-meta-source/2026-11-19",
      "Accept":"application/json,application/zip,text/csv,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=24) as f:
        b=f.read(MAX_BYTES+1)
        if len(b)>MAX_BYTES:raise ValueError("Overlarge original pollen source")
        return b


def parse_table_header(name,raw):
    suffix=Path(name).suffix.lower()
    if suffix in (".csv",".tsv"):
        delimiter="\t" if suffix==".tsv" else ","
        text=raw.decode("utf-8-sig",errors="replace")
        rr=list(csv.reader(io.StringIO(text,newline=""),delimiter=delimiter))
        if not rr:raise ValueError("Published original CSV empty")
        h=rr[0]
        return {"format":suffix,"original_column_names":h,
                "original_nonempty_data_rows":sum(any(str(x).strip() for x in row) for row in rr[1:]),
                "original_biological_values_not_exported":True}
    if suffix==".xlsx":
        from openpyxl import load_workbook
        if not raw.startswith(b"PK\x03\x04"):raise ValueError("Not original XLSX file")
        wb=load_workbook(io.BytesIO(raw),read_only=True,data_only=True)
        summary=[]
        for sh in wb.worksheets:
            rows=iter(sh.iter_rows(values_only=True))
            header=list(next(rows,[]))
            count=sum(any(v is not None and str(v).strip() for v in row) for row in rows)
            summary.append({"sheet":sh.title,
                            "first_row_column_names_or_cells":[str(x)[:100] if x is not None else "" for x in header],
                            "original_nonempty_rows_after_first":count,
                            "original_records_not_exported":True})
        wb.close()
        return {"format":suffix,"original_sheets_shape_only":summary}
    if suffix==".zip":
        if not raw.startswith(b"PK\x03\x04"):raise ValueError("Not original zip")
        out=[]
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            members=z.infolist()
            if len(members)>80 or sum(x.file_size for x in members)>MAX_BYTES:
                raise ValueError("Unexpectedly large nested source")
            for x in members:
                if not x.is_dir() and Path(x.filename).suffix.lower() in (".csv",".tsv",".xlsx"):
                    bb=z.read(x)
                    out.append({"member_name":x.filename,"member_SHA256":hashlib.sha256(bb).hexdigest(),
                                "member_schema":parse_table_header(x.filename,bb)})
        return {"format":suffix,"source_members":out}
    return {"format":suffix or "unknown","numeric_source_parser_status":"FILE_EXTENSION_UNSUPPORTED"}


def original_source_preflight():
    d={
      "schema":"alan_stageC_original_moth_fly_q_two_existing_meta_primary_datasets_preflight_v1",
      "original_meta_study_DOI":DOI_ORIGINAL_META,
      "original_meta_dryad_DOI":DRYAD_DOI,
      "original_meta_published_28_studies":28,
      "original_meta_published_127_observations":127,
      "original_meta_published_bee_taxa":94,
      "original_meta_published_fly_taxa":33,
      "original_meta_published_moth_taxa":0,
      "meta_source_original_format_status":"SOURCE_STOP",
      "original_primary_DOI":DOI_ORIGINAL_PRIMARY,
      "original_primary_figshare_ID":FIGSHARE_ID,
      "source_primary_pollinator_groups":["butterflies and diurnal moths","hoverflies","bumblebees"],
      "source_primary_nocturnal_moth_singlevisit_measurement_verified":False,
      "primary_source_original_format_status":"SOURCE_STOP",
      "no_moth_fly_same_plant_efficiency_ratio_estimated":True,
      "no_lu2026_plant_taxon_SVD_efficiency_verified":True,
      "all_original_single_visit_observation_rows_omitted":True,
      "any_completed_stageC_meta_analysis":False,
    }
    doi=urllib.parse.quote("doi:"+DRYAD_DOI,safe="")
    try:
        b=get(f"https://datadryad.org/api/v2/datasets/{doi}/download")
        if not b.startswith(b"PK\x03\x04"):raise ValueError("Original Dryad source not ZIP")
        d["meta_source_original_format_status"]="GENUINE_ARCHIVE_VERIFIED"
        d["meta_original_Dryad_source_ZIP_SHA256"]=hashlib.sha256(b).hexdigest()
        d["meta_source_original_file_summary"]=parse_table_header("source-dryad.zip",b)
    except (ValueError,UnicodeError,urllib.error.URLError,TimeoutError,zipfile.BadZipFile,ImportError) as exc:
        d["meta_source_original_format_status"]="SOURCE_STOP"
        d["meta_source_error_class"]=type(exc).__name__
        if isinstance(exc,urllib.error.HTTPError):d["meta_source_http_status"]=exc.code
    try:
        b=get(f"https://api.figshare.com/v2/articles/{FIGSHARE_ID}")
        meta=json.loads(b)
        if int(meta.get("id",-1))!=FIGSHARE_ID:raise ValueError("Original Figshare record ID mismatch")
        d["primary_source_metadata_SHA256"]=hashlib.sha256(b).hexdigest()
        d["primary_source_title"]=meta.get("title")
        files=meta.get("files",[])
        records=[]
        for entry in files[:40]:
            data={"source_file_ID":entry.get("id"),
                  "source_file_name":entry.get("name"),
                  "declared_size":entry.get("size"),
                  "source_download_status":"SOURCE_STOP"}
            url=entry.get("download_url","")
            size=entry.get("size",None)
            if isinstance(size,int) and 0<size<=MAX_BYTES and isinstance(url,str) and url.startswith("https://"):
                try:
                    file=get(url)
                    data["source_download_status"]="ACTUAL_ORIGINAL_BYTES_VERIFIED"
                    data["original_file_SHA256"]=hashlib.sha256(file).hexdigest()
                    data["file_schema"]=parse_table_header(str(entry.get("name") or ""),file)
                except (ValueError,UnicodeError,urllib.error.URLError,TimeoutError,zipfile.BadZipFile,ImportError) as ex:
                    data["source_error_class"]=type(ex).__name__
                    if isinstance(ex,urllib.error.HTTPError):data["source_http_status"]=ex.code
            records.append(data)
        d["primary_source_original_format_status"]="ORIGINAL_FIGSHARE_METADATA_VERIFIED"
        d["primary_source_file_summary"]=records
        if any(r["source_download_status"]=="ACTUAL_ORIGINAL_BYTES_VERIFIED" for r in records):
            d["primary_source_original_format_status"]="ORIGINAL_FILE_BYTES_VERIFIED"
    except (ValueError,KeyError,urllib.error.URLError,TimeoutError,ImportError) as exc:
        d["primary_source_original_format_status"]="SOURCE_STOP"
        d["primary_source_error_class"]=type(exc).__name__
        if isinstance(exc,urllib.error.HTTPError):d["primary_source_http_status"]=exc.code
    return d


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z=original_source_preflight()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_STAGEC_SVD_PUBLIC_META_SOURCE_STATUS",z["meta_source_original_format_status"])
    print("ALAN_STAGEC_SVD_PUBLIC_META_SOURCE_TABLE_STRUCTURE",z.get("meta_source_original_file_summary"))
    print("ALAN_STAGEC_KNAUTIA_HOVERFLY_BUTTERFLY_ORIGINAL",z["primary_source_original_format_status"],z.get("primary_source_file_summary"))
    print("ALAN_STAGEC_TRUE_SAME_PLANT_MOTH_FLY_RATIO_NOT_MEASURED",z["no_moth_fly_same_plant_efficiency_ratio_estimated"])


if __name__=="__main__":
    main()
