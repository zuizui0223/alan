#!/usr/bin/env python3
"""Acquire public ALAN field experiment *manifest* without republishing raw files.

Primary original Lu et al. 2026 Ecology and Evolution DOI 10.1002/ece3.74368,
PMCID PMC13645961, publisher supporting Data S1 ZIP, original Data.xlsx.
Secondary Swiss Grognuz et al. 2025 Mendeley DOI 10.17632/xjmn5jhxyn.2.
This source discovery reports real workbook names, sheets, header fields, sizes,
checksums, file lookup statuses. It is NOT a pollination or phenology meta-model.
Never mark missing public API access as ecological absence, and do not put raw
copyrighted original study rows into repository artifacts.
"""
from __future__ import annotations
import argparse, hashlib, json, io, urllib.request, urllib.error, zipfile
from pathlib import Path
from datetime import datetime, timezone

PUBLISHED_SUPP = [
 "https://pmc.ncbi.nlm.nih.gov/articles/PMC13645961/bin/ece374368-sup-0002-DataS1.zip",
 "https://pmc.ncbi.nlm.nih.gov/articles/instance/13645961/bin/ece374368-sup-0002-DataS1.zip",
 "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles"
]
SWISS_ENDPOINTS = [
 "https://api.data.mendeley.com/datasets/xjmn5jhxyn?version=2",
 "https://api.data.mendeley.com/datasets/publics/xjmn5jhxyn/files?version=2",
 "https://api.data.mendeley.com/datasets/xjmn5jhxyn/files?version=2"
]
HEADERS={"User-Agent":"ALAN-open-primary-repro-source/1.0 (scientific metadata QA)",
         "Accept":"application/json, application/zip, */*"}

def request(url, max_bytes, timeout=28):
    r=urllib.request.Request(url,headers=HEADERS)
    try:
        with urllib.request.urlopen(r,timeout=timeout) as f:
            data=f.read(max_bytes+1)
            if len(data)>max_bytes:
                return None,{"status":"byte_limit_exceeded","bound":max_bytes}
            return data,{"status":"downloaded","bytes":len(data),
                         "sha256":hashlib.sha256(data).hexdigest(),
                         "content_type":str(f.headers.get("Content-Type",""))[:100]}
    except urllib.error.HTTPError as e:
        return None,{"status":"http_"+str(e.code)}
    except (urllib.error.URLError, TimeoutError,OSError,ValueError) as e:
        return None,{"status":"network_unavailable","error_type":type(e).__name__}

def xlsx_manifest(b, size_limit=120_000_000):
    try:
        from openpyxl import load_workbook
        w=load_workbook(io.BytesIO(b),read_only=True,data_only=True)
        info=[]
        for s in w.worksheets:
            header=[]
            # Schema only, never public source plant×site records or data rows.
            for x in s.iter_rows(min_row=1,max_row=2,values_only=True):
                header=[str(v)[:125] if v is not None else None for v in x[:70]]
                if sum(v is not None for v in x)>=2:
                    break
            info.append({"sheet":s.title,"rows_including_header":s.max_row,
                         "columns":s.max_column,"first_readable_header_row":header,
                         "sample_row_values_included":False})
        w.close()
        return {"status":"valid_xlsx","sheets":info}
    except Exception as e:
        return {"status":"xlsx_read_error","error_type":type(e).__name__}

def manifest_zip_bytes(blob,label):
    if not zipfile.is_zipfile(io.BytesIO(blob)):
        return {"status":"not_zip","source":label,"bytes":len(blob)}
    f=zipfile.ZipFile(io.BytesIO(blob))
    entries=[]
    for item in f.infolist():
        if item.is_dir():continue
        z={"filename":item.filename,"uncompressed_size":item.file_size,
           "compressed_size":item.compress_size,
           "CRC32":format(item.CRC,"08x")}
        if item.file_size>120_000_000:
            z["skipped"]="oversized_uncompressed_file"
        elif item.filename.lower().endswith((".xlsx",".xlsm")):
            raw=f.read(item)
            z["sha256"]=hashlib.sha256(raw).hexdigest()
            z["workbook"]=xlsx_manifest(raw)
        entries.append(z)
    f.close()
    return {"status":"valid_zip","source":label,"original_entries":entries}

def main(out):
    out=Path(out)
    if out.exists():raise FileExistsError("Refuse overwriting original source audit")
    lu={"published_article_DOI":"10.1002/ece3.74368",
        "PMCID":"PMC13645961","original_supplement_filename":
           "ece374368-sup-0002-DataS1.zip",
        "candidate_downloads":[]}
    for url in PUBLISHED_SUPP:
        data,receipt=request(url,28_000_000)
        lu["candidate_downloads"].append({"source_url":url,**receipt})
        if data is not None and zipfile.is_zipfile(io.BytesIO(data)):
            lu["original_source"]=manifest_zip_bytes(data,url)
            lu["actual_successful_source_url"]=url
            break
    if "original_source" not in lu:
        lu["original_source"]={"status":"NO_VERIFIED_PUBLIC_ZIP_ACQUIRED",
                               "no_synthetic_network_results":True}
    swiss={"published_article_DOI":"10.1016/j.biocon.2025.111650",
           "published_raw_dataset_DOI":"10.17632/xjmn5jhxyn.2",
           "original_source_version":2,
           "API_candidate_receipts":[]}
    for url in SWISS_ENDPOINTS:
        data,receipt=request(url,3_000_000)
        a={"source_url":url,**receipt}
        if data is not None:
            try:
                j=json.loads(data.decode("utf-8"))
                if isinstance(j,dict):
                    a["parsed_dataset_title"]=str(j.get("name",""))[:190]
                    a["parsed_version"]=j.get("version")
                    files=j.get("files")
                else: files=j
                if isinstance(files,list):
                    a["public_file_manifest"]=[{
                      "filename":str(z.get("filename",""))[:150],
                      "id":str(z.get("id",""))[:100],
                      "size":z.get("size"),
                      "sha256":(z.get("content_details") or {}).get("sha256_hash")
                    }for z in files[:150] if isinstance(z,dict)]
                    a["n_source_files"]=len(files)
                a["json_parse"]="ok"
            except (ValueError,UnicodeError) as e:
                a["json_parse_error_type"]=type(e).__name__
        swiss["API_candidate_receipts"].append(a)
    result={"schema":"ALAN_primary_manipulated_light_daynight_syndrome_rawsource_manifest_v1",
         "status":"ACTUAL_PUBLIC_SOURCE_PROBE_ONLY_NO_MODEL_ESTIMATED",
         "performed_utc":datetime.now(timezone.utc).isoformat(),
         "sources":{"Lu2026":lu,"Swiss_ALAN_phenology":swiss},
         "source_gates":{"real_local_altered_Lu_visitation_NOT_real_seed_set":True,
             "do_not_derive_true_syndrome_from_visiting_moth_name_only":True,
             "no_automatic_same_site_year_phenology_overlap":True,
             "Mendeley_API_access_may_require_auth":True,
             "never_infer_syndrome_sensitivity_from_differential_observation_effort":True,
             "no_user_provided_private_source_values_exported":True,
             "not_causal_effect_of_predator_or_island_shift":True}
         }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ORIGINAL_LU_SWISS_PRIMARY_PUBLIC_DATA_MANIFEST",
          json.dumps({"Lu2026":{"source":lu["original_source"],
                                 "receipts":lu["candidate_downloads"]},
                     "Swiss":swiss,"status":result["status"]},
                     ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    main(a.out)
