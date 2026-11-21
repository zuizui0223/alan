#!/usr/bin/env python3
"""Source-only gate for 2026 freshwater detritivore ALAN experiments (3 taxa).

Published Zenodo DOI 10.5281/zenodo.19022243 referenced in
Frontiers Ecology and Evolution 2026 DOI 10.3389/fevo.2026.1830790.

Read and SHA official public record metadata, safe-size original files.
No effect size/SE is calculated until per-treatment raw records, units,
replication hierarchy and control identity can be reconstructed.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, urllib.request, urllib.error, zipfile
from pathlib import Path
SCHEMA="alan_zenodo19022243_freshwater_feeding_original_source_gate_v1"
RECORD=19022243
META=f"https://zenodo.org/api/records/{RECORD}"
MAXB=6<<20

def get(url):
    req=urllib.request.Request(url,headers={
      "User-Agent":"ALAN ecological mechanism source audit 2026",
      "Accept":"application/json,text/csv,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=60) as r:
        b=r.read(MAXB+1)
        if len(b)>MAXB:raise RuntimeError("SIZE_6MB_STOP")
        return b,r.geturl()

def inspect(name,buf):
    extension=Path(name).suffix.lower()
    if extension in (".csv",".tsv",".txt"):
        for enc in ("utf-8-sig","latin1"):
            try:s=buf.decode(enc)
            except UnicodeError:continue
            for sep in (",","\t",";"):
                data=list(csv.DictReader(io.StringIO(s),delimiter=sep))
                if not data or len(data[0])<3:continue
                keys=list(data[0])
                return {"type":"source_delimited","rows":len(data),"columns":keys[:120],
                        "separator":sep,"encoding":enc,
                        "candidate_factor_levels":{
                         k:sorted({r[k].strip() for r in data})[:30]
                         for k in keys
                         if ("light" in k.lower() or "species" in k.lower()
                             or "treatment" in k.lower() or "lux" in k.lower())
                            and len({r[k].strip() for r in data})<25
                        }}
        return {"type":"text_unparsed"}
    if extension==".xlsx":
        with zipfile.ZipFile(io.BytesIO(buf)) as z:
            return {"type":"xlsx","excel_sheet_blob_count":sum(
                q.startswith("xl/worksheets/") and q.endswith(".xml")
                for q in z.namelist()),"bytes":len(buf)}
    if extension==".zip":
        with zipfile.ZipFile(io.BytesIO(buf)) as z:
            return {"type":"zip","archive_members":[{
                "name":q.filename,"size":q.file_size} for q in z.infolist()[:50]]}
    return {"type":"original_other_binary","extension":extension,"bytes":len(buf)}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    d={"schema":SCHEMA,"source_doi":"10.5281/zenodo.19022243",
       "experiment_article_doi":"10.3389/fevo.2026.1830790",
       "species":"Potamopyrgus antipodarum;Sericostoma vittatum;Chironomus riparius",
       "original_source_url":META,"no_new_effect_size_fitted":True,
       "wet_plant_moisture_moderator_unavailable":True}
    try:
        data,final=get(META)
        meta=json.loads(data)
        d.update(status="METADATA_VERIFIED",meta_sha256=hashlib.sha256(data).hexdigest(),
                 verified_url=final,metadata_title=meta.get("metadata",{}).get("title"))
        fs=meta.get("files",[])
        d["files"]=[]
        for item in fs:
            fname=item.get("key",item.get("filename",""))
            url=(item.get("links") or {}).get("self") or (item.get("links") or {}).get("content")
            entry={"filename":fname,"original_size":item.get("size"),
                   "original_checksum":item.get("checksum"),
                   "source_original_link":url}
            if url and item.get("size",0)<=MAXB:
                try:
                    blob,resolved=get(url)
                    entry.update(status="SOURCE_BYTES_VERIFIED",sha256=hashlib.sha256(blob).hexdigest(),
                         schema=inspect(fname,blob),resolved_url=resolved)
                except Exception as e:
                    entry.update(status="SOURCE_BYTE_DOWNLOAD_STOP",
                                 error_type=type(e).__name__,
                                 http_code=getattr(e,"code",None))
            else:entry["status"]="NO_OR_LARGE_FILE_SOURCE_STOP"
            d["files"].append(entry)
        if len(fs)==0:d["status"]="NO_SOURCE_FILES_STOP"
    except Exception as e:
        d.update(status="ZENODO_METADATA_ACCESS_STOP",
                 reason_type=type(e).__name__,
                 http_code=getattr(e,"code",None))
    x=Path(a.out);x.parent.mkdir(parents=True,exist_ok=True)
    x.write_text(json.dumps(d,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("FRESHWATER2026_ZENODO_SOURCE_GATE "+json.dumps(d,sort_keys=True),flush=True)
if __name__=="__main__":main()
