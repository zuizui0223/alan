#!/usr/bin/env python3
"""Official Dryad Bucher et al 2023 species biomass source availability and grain gate.

DOI 10.5061/dryad.hhmgqnknt:
 ALAN_specis_biomass.csv (published original filename), file 2619566
 ALAN_biomass.csv 2619567
 ALAN_traits.csv 2619568
 README.md 2619572

Read public original data only. Do not infer plant responses, wetness effects or
replicate count until actual species × Ecotron unit × light × harvest raw rows
and exact species-name mapping have been source verified.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,urllib.request,urllib.error
from pathlib import Path
SCHEMA="alan_bucher2023_dryad_original_species_source_gate_20261120_v1"
URL="https://datadryad.org"
FILES={"ALAN_specis_biomass.csv":2619566,
       "ALAN_biomass.csv":2619567,
       "ALAN_traits.csv":2619568,
       "README.md":2619572}
PREFIX="https://datadryad.org/api/v2"
MAX=1024*1024

def get(u):
    request=urllib.request.Request(u,headers={
       "User-Agent":"Ecological research reproducibility alan/2026 OpenAI",
       "Accept":"text/csv,text/plain,application/octet-stream,application/json,*/*"})
    with urllib.request.urlopen(request,timeout=25) as response:
        b=response.read(MAX+1)
        if len(b)>MAX:raise RuntimeError("SOURCE_LIMIT_1M_STOP")
        return b,response.geturl(),response.headers.get("Content-Type","")

def inspect(name,raw):
    if name.endswith(".md"):
        return {"type":"markdown","bytes":len(raw),
                "first_line":raw.decode("utf-8-sig",errors="replace").splitlines()[0][:130]}
    s=raw.decode("utf-8-sig")
    if s.lstrip().lower().startswith(("<!doctype","<html")):
        raise RuntimeError("SOURCE_HTML_NOT_CSV")
    for sep in (",",";","\t"):
        rows=list(csv.DictReader(io.StringIO(s),delimiter=sep))
        if len(rows)>0 and len(rows[0])>2:
            keys=list(rows[0])
            r={"type":"original_delimited","rows":len(rows),"columns":keys,"delimiter":sep}
            for col in keys:
                if col.lower().strip() in ("unit","ecounit","light","treatment","harvest","species"):
                    values={str(x.get(col,"")).strip() for x in rows}
                    r[col]={"levels":len(values),"unique_values":sorted(values)[:120]}
            return r
    raise RuntimeError("NO_TRUSTWORTHY_TABULAR_DATA")

def execute():
    report={"schema":SCHEMA,
      "doi":"10.5061/dryad.hhmgqnknt",
      "article":"10.1098/rstb.2022.0358",
      "dryad_url":URL+"/dataset/doi:10.5061/dryad.hhmgqnknt",
      "published_file_identifiers":FILES,
      "original_ecotron_units_not_independent_plants_or_harvest_rows":True,
      "original_raw_rows_committed_to_git":False,
      "wet_trait_coefficient_fitted":False,
      "species_light_interaction_fitted":False,
      "files":{}}
    for name,id in FILES.items():
        options=[
           f"{URL}/stash/downloads/file_stream/{id}",
           f"{URL}/stash/downloads/file_stream/{id}?download=1",
           f"{URL}/downloads/file_stream/{id}",
           f"{PREFIX}/files/{id}/download",
           f"{PREFIX}/files/{id}/download?download=1",
        ]
        trials=[]
        found=None
        for u in options:
            try:
                b,final,ctype=get(u)
                d=inspect(name,b)
                found={"status":"ORIGINAL_SOURCE_SCHEME_VERIFIED","sha256":hashlib.sha256(b).hexdigest(),
                       "size":len(b),"content_type":ctype,"schema":d,
                       "original_source_url":u,"resolved_url":final}
                break
            except Exception as e:
                trials.append({"endpoint":u,"error_type":type(e).__name__,
                    "status_code":getattr(e,"code",None)})
        report["files"][name]=found or {
           "status":"SOURCE_ACCESS_STOP","endpoint_trials":trials}
    success=all(v["status"]=="ORIGINAL_SOURCE_SCHEME_VERIFIED" for v in report["files"].values())
    report["status"]="SOURCE_AUDIT_COMPLETE" if success else "SOURCE_ACCESS_INCOMPLETE_STOP"
    if success:
        s=report["files"]["ALAN_specis_biomass.csv"]["schema"]
        report["source_measured_species_labels"]=s.get("species",{}).get("unique_values",[])
        report["source_measured_unique_ecounits"]=s.get("unit",{}).get("levels")
        report["source_light_regimes"]=s.get("light",{}).get("unique_values",[])
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    o=execute()
    dest=Path(args.out);dest.parent.mkdir(exist_ok=True,parents=True)
    dest.write_text(json.dumps(o,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("BUCHER_DRYAD_ORIGINAL_ACCESS_AND_ECOUNIT_SOURCE_GATE "+json.dumps(o,sort_keys=True),flush=True)
if __name__=="__main__":
    main()
