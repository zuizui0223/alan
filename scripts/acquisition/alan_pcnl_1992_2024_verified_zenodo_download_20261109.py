#!/usr/bin/env python3
"""Download independently published PCNL v3 annual rasters with fixed Zenodo MD5.

The output is RAW annual 1992-2024 harmonized DMSP-OLS/VIIRS raster data,
not photo colour or a biological endpoint. Store raster originals only in
ephemeral CI; preserve checksums, byte sizes and source metadata in receipts.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import urllib.error
import urllib.request
import time

RECORD="https://zenodo.org/records/17013414"
SCHEMA="ALAN_PCNLPrior1992_2024_Zenodo_record17013414_MD5_v1"

def require(ok, msg):
    if not ok:
        raise ValueError("ALAN_PCNL_REAL_RASTER: "+str(msg))

def md5(path):
    h=hashlib.md5()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):
            h.update(block)
    return h.hexdigest()

def manifest_read(path):
    a=json.loads(Path(path).read_text())
    require(a.get("schema")==SCHEMA and
            a.get("source_url")==RECORD and
            a.get("source_doi")=="10.5281/zenodo.17013414" and
            set(a.get("files",{}))=={f"PCNL{yr}.tif" for yr in range(1992,2025)},
            "Noncanonical PCNL v3 1992-2024 manifest")
    require(all(len(x)==32 and all(c in "0123456789abcdef" for c in x)
                for x in a["files"].values()),"Bad published MD5")
    return a

def fetch_one(name,expected,target):
    target=Path(target)
    if target.is_file():
        require(md5(target)==expected,"Preexisting source file MD5 incorrect")
        return {"name":name,"size_bytes":target.stat().st_size,
                "md5":expected,"status":"verified_existing"}
    url=f"{RECORD}/files/{name}?download=1"
    temporary=target.with_suffix(".tmp")
    for trial in range(3):
        h=hashlib.md5()
        n=0
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"ALAN-Macroecology-Source-Audit/1.0"})
            with urllib.request.urlopen(req,timeout=120) as response:
                require(response.status==200,"Zenodo served incomplete HTTP")
                with temporary.open("wb") as f:
                    while True:
                        part=response.read(4<<20)
                        if not part:break
                        f.write(part);h.update(part);n+=len(part)
            require(n>1_000_000,"Not an original PCNL raster download")
            require(h.hexdigest()==expected,
                    f"Published MD5 mismatch {name}, got {h.hexdigest()}")
            temporary.rename(target)
            return {"name":name,"size_bytes":n,
                    "md5":expected,"status":"verified_download"}
        except (OSError,urllib.error.URLError,ValueError) as exc:
            temporary.unlink(missing_ok=True)
            if trial==2:raise
            time.sleep(1+trial)
    raise AssertionError("Unreachable downloader guard")

def run(manifest,directory,receipt):
    m=manifest_read(manifest)
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    require(not Path(receipt).exists(),"Refuse replace existing PCNL receipt")
    values=[]
    for name,expected in m["files"].items():
        item=fetch_one(name,expected,directory/name)
        values.append(item)
        print("ALAN_REAL_PCNL_YEAR_MD5_VERIFIED",json.dumps(item,sort_keys=True),flush=True)
    require(len(values)==33,"Original annual 1992-2024 PCNL not complete")
    output={"schema":"ALAN_real_PCNL_1992_2024_complete_33_raster_MD5_receipt_v1",
            "record":RECORD,"doi":"10.5281/zenodo.17013414",
            "year_rasters":values,"original_raster_year_count":33,
            "raster_byte_sum":sum(x["size_bytes"] for x in values),
            "unaligned_raw_DMSP_VIIRS_not_combined":True,
            "raster_origin_unchanged_and_no_photo_ROI":True}
    Path(receipt).write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
    return output

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    a=p.parse_args()
    run(a.manifest,a.directory,a.receipt)

if __name__=="__main__":
    main()
