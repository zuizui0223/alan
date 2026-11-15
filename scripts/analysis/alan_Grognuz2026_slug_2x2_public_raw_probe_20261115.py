#!/usr/bin/env python3
"""Access and inventory the exact public raw-data files for a direct 2x2
nightlight x nocturnal slug manipulation with measured daytime visits and seeds.

Grognuz & Knop (2026) publicly indexed Mendeley dataset:
10.17632/drhzznh4ps.1, published 2026-04-10 (CC BY 4.0).
Prior audit PR75 had author-reported paths only because Mendeley API returned
401. This is source discovery, not fabricated refit or missing=0 treatment.

When raw files can be authenticated, download and inspect their structure,
not results; preserve raw data only on ephemeral Actions runners. Files are
never stored or copied into repository/artifacts. Public URLs only; no auth.
"""
from __future__ import annotations
import argparse,concurrent.futures,hashlib,io,json,re,time,urllib.error,urllib.request,urllib.parse,zipfile
from pathlib import Path
from html.parser import HTMLParser

DATASET_ID="drhzznh4ps"
DOI="10.17632/drhzznh4ps.1"
ENDPOINTS={
    "v1_HTML":"https://data.mendeley.com/datasets/drhzznh4ps/1",
    "public_dataset_v1":"https://api.data.mendeley.com/datasets/drhzznh4ps?version=1",
    "public_dataset_snapshot":"https://api.data.mendeley.com/datasets/drhzznh4ps/snapshot?version=1",
    "public_dataset_files":"https://api.data.mendeley.com/datasets/publics/drhzznh4ps/files?version=1",
    "public_dataset_files_alt":"https://api.data.mendeley.com/datasets/drhzznh4ps/files?version=1",
    "public_dataset_zip":"https://data.mendeley.com/public-files/datasets/drhzznh4ps/zip/file_downloaded",
    "api_zip":"https://api.data.mendeley.com/datasets/drhzznh4ps/zip/file_downloaded",
    "article_preprint_repository":"https://data.mendeley.com/datasets/drhzznh4ps",
}
UA="ALAN-independent-public-2x2-reproduction-source-audit/1.0"
MAX_BYTES=38_000_000

def request(url,limit=MAX_BYTES):
    if not isinstance(url,str) or not url.startswith(("https://","http://")):
        return None,{"status":"URL_ERROR","http_code":None,"error_kind":"ValueError"}
    req=urllib.request.Request(url,headers={
        "User-Agent":UA,
        "Accept":"application/vnd.mendeley-public-dataset.1+json,application/json,text/html,application/zip,*/*",
    })
    try:
        with urllib.request.urlopen(req,timeout=38) as r:
            raw=r.read(limit+1)
            if len(raw)>limit:
                return None,{"status":"byte_limit_exceeded","size_bound":limit}
            return raw,{"status":"HTTP_"+str(r.status),"http_code":r.status,
                "downloaded_bytes":len(raw),"final_host":urllib.parse.urlparse(r.url).hostname,
                "content_type":r.headers.get("Content-Type"),
                "sha256":hashlib.sha256(raw).hexdigest()}
    except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,OSError,ValueError) as e:
        return None,{"status":"HTTP_ERROR" if hasattr(e,"code") else "URL_ERROR",
                     "http_code":getattr(e,"code",None),"error_kind":type(e).__name__}

def retrieve_api():
    output=[]
    files={}
    candidate_id_urls=[]
    def fetch_one(item):
        k,url=item
        return k,url,*request(url,38_000_000 if "zip" in k.lower() else 5_000_000)
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        items=list(ex.map(fetch_one,ENDPOINTS.items()))
    for k,url,blob,receipt in items:
        meta={"source":k,"requested_public_url":url,**receipt}
        if blob:
            if zipfile.is_zipfile(io.BytesIO(blob)):
                with zipfile.ZipFile(io.BytesIO(blob)) as z:
                    entries=[{"filename":p.filename,"bytes":p.file_size,
                             "is_dir":p.is_dir()}
                             for p in z.infolist()[:300]]
                meta["public_zip_original_file_manifest"]=entries
                meta["raw_ZIP_acquired_and_ephemeral_only"]=True
            else:
                try:
                    v=json.loads(blob.decode("utf-8"))
                except (ValueError,UnicodeError):v=None
                if isinstance(v,(dict,list)):
                    meta["returned_JSON"]=True
                    meta["top_keys"]=list(v)[:30] if isinstance(v,dict) else []
                    if isinstance(v,dict):
                        meta["dataset_title"]=v.get("name") or v.get("title")
                        meta["dataset_version"]=v.get("version")
                        meta["status_field"]=v.get("status")
                        for subfield in ("files","content","items","results"):
                            if isinstance(v.get(subfield),list):
                                v=v[subfield];break
                    if isinstance(v,list):
                        meta["n_list_entries"]=len(v)
                        for obj in v[:100]:
                            if not isinstance(obj,dict):continue
                            nm=str(obj.get("filename",obj.get("name","")))
                            fileid=obj.get("id")
                            sha=(obj.get("content_details") or {}).get("sha256_hash")
                            size=obj.get("size") or (obj.get("content_details") or {}).get("size")
                            files[(str(fileid),nm)]={"filename":nm,"id":fileid,
                                "bytes_reported":size,"sha256_reported":sha,
                                "source_api":k}
                elif blob[:30].lstrip().startswith(b"<!") or b"<html" in blob[:300].lower():
                    text=blob.decode("utf-8",errors="ignore")
                    # The dataset JSON may be embedded as NextData, not an
                    # official API payload. Do NOT guess file ids from filenames.
                    ids=sorted(set(re.findall(
                        r'(?:file[_-]?id|/files/)[/"=: ]+([0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12})',
                        text,re.IGNORECASE)))
                    meta["html_public_file_guid_candidates"]=ids[:30]
                    meta["has_files_section"]="files" in text.lower()
        output.append(meta)
    return output,list(files.values())

def main(out):
    path=Path(out)
    if path.exists():raise FileExistsError("source audit would overwrite")
    receipts,metadata=retrieve_api()
    eligible=[m for m in metadata if m["id"] and m["filename"]]
    resolved=[]
    for f in eligible[:30]:
        ident=str(f["id"])
        if not re.fullmatch(r"[0-9a-fA-F-]{36}",ident):continue
        uri=f"https://data.mendeley.com/public-files/datasets/{DATASET_ID}/files/{ident}/file_downloaded"
        raw,receipt=request(uri,MAX_BYTES)
        item={**f,"direct_public_download":uri,**receipt}
        if raw:
            item["raw_file_acquired_and_ephemeral_only"]=True
            if f["sha256_reported"]:
                item["publisher_sha256_verified"]=(receipt.get("sha256")==f["sha256_reported"])
            if zipfile.is_zipfile(io.BytesIO(raw)):
                with zipfile.ZipFile(io.BytesIO(raw)) as z:
                    item["zip_manifest_only"]=[{"filename":x.filename,"bytes":x.file_size} for x in z.infolist()[:50]]
            if raw.startswith((b"PK",b"<?xml",b"{",b"[")):
                item["file_magic"]="zip_or_XML_or_JSON"
            elif b"," in raw[:300]:
                item["file_magic"]="csv_candidate"
            elif b"\t" in raw[:300]:
                item["file_magic"]="tsv_candidate"
        resolved.append(item)
    result={"schema":"ALAN_GrognuzKnop_2026_2x2_slug_light_seed_public_raw_discovery_v1",
      "dataset_DOI":DOI,"manuscript_DOI":"10.13140/RG.2.2.36822.10563",
      "published_dataset_version":1,"source_receipts":receipts,
      "api_listed_files":metadata,"requested_actual_files":resolved,
      "confirmed_public_authentic_raw_files":sum(
        bool(x.get("raw_file_acquired_and_ephemeral_only") and
        x.get("publisher_sha256_verified",True)) for x in resolved),
      "raw_plantlevel_seeds_and_pollen_data_refitted_here":False,
      "preprint_peerreview_status":"NOT_CONFIRMED",
      "missing_data_not_zero":True,
      "null_source_api_access_is_NOT_ecological_result":True,
      "no_synthetic_seed_effect_or_fruit_count":True}
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GROGNUZ_2026_PUBLIC_MANUSCRIPT_DATA_ACCESS",json.dumps({
        "endpoints":{x["source"]:(x["status"],x.get("http_code"),x.get("downloaded_bytes")) for x in receipts},
        "api_listed_files":len(metadata),"raw_downloads":result["confirmed_public_authentic_raw_files"],
        "filenames":[x["filename"] for x in metadata]},ensure_ascii=False))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    main(p.parse_args().out)
