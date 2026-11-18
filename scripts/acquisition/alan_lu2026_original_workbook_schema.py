#!/usr/bin/env python3
"""Read genuine Lu 2026 original Data S1: XLSX XML source schema, no imputation.

Source Europe PMC / PMC13645961 supplemental package. Only author-original
Date.xlsx (actual original typo), sheet metadata and readme are inspected.
The original five lit/control spatial blocks are units, NOT thousands of
recorded plant-pollinator visitation events. Source R code stays as published.
No pooled cross-article ALAN effect is estimated by this preflight.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
import urllib.request
import zipfile
from pathlib import PurePosixPath,Path
from xml.etree import ElementTree as ET

SUPPLEMENT_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13645961/supplementaryFiles"
ORIGINAL_WORKBOOK_SHA="ba7e15e66986f64943a38fd294e49d7128cf7719d318578f6b805d29cde8f98f"
ORIGINAL_SUPP_CRC="6e247a04"
ORIGINAL_INNER_CRC="54082a49"
SCHEMA="alan_lu2026_original_fiveblock_plant_pollinator_sheet_schema_v1"
MAIN="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
DOCREL="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE="http://schemas.openxmlformats.org/package/2006/relationships"

def source_open():
    with urllib.request.urlopen(urllib.request.Request(
        SUPPLEMENT_URL,headers={"User-Agent":"ALAN source-verification ecology"}),timeout=90) as r:
        data=r.read(15<<20)
    assert zipfile.is_zipfile(io.BytesIO(data)),"EuropePMC original source not an archive"
    with zipfile.ZipFile(io.BytesIO(data)) as f:
        assert f.testzip() is None
        target=f.getinfo("ECE3-16-e74368-s001.zip")
        assert f"{target.CRC:08x}"==ORIGINAL_SUPP_CRC
        nested=f.read(target)
    with zipfile.ZipFile(io.BytesIO(nested)) as z:
        assert z.testzip() is None
        x=z.getinfo("Date.xlsx")
        assert f"{x.CRC:08x}"==ORIGINAL_INNER_CRC
        workbook=z.read(x)
        rscript=z.read("R analysis-ECE.R").decode("utf-8-sig",errors="replace")
    assert hashlib.sha256(workbook).hexdigest()==ORIGINAL_WORKBOOK_SHA
    return workbook,rscript,{"original_pmc_supplement_url":SUPPLEMENT_URL,
        "source_workbook_sha256":ORIGINAL_WORKBOOK_SHA,
        "source_raw_bytes_excel":len(workbook),
        "source_R_code_bytes":len(rscript.encode("utf-8"))}

def sheet_data(book):
    with zipfile.ZipFile(io.BytesIO(book)) as z:
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        paths={p.attrib["Id"]:p.attrib["Target"] for p in rel.findall(f"{{{PACKAGE}}}Relationship")}
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{{{MAIN}}}si"):
                strings.append("".join(t.text or "" for t in si.findall(f".//{{{MAIN}}}t")))
        sheets={}
        for sh in wb.findall(f".//{{{MAIN}}}sheet"):
            name=sh.attrib["name"]
            key=sh.attrib[f"{{{DOCREL}}}id"]
            target=paths[key].lstrip("/")
            if not target.startswith("xl/"):target="xl/"+target
            root=ET.fromstring(z.read(target))
            rows=[]
            for rr in root.findall(f".//{{{MAIN}}}sheetData/{{{MAIN}}}row"):
                vals={}
                for c in rr.findall(f"{{{MAIN}}}c"):
                    address=c.attrib["r"]
                    col=re.match(r"[A-Z]+",address)
                    if col is None:continue
                    i=0
                    for chr in col.group():
                        i=i*26+ord(chr)-64
                    i-=1
                    t=c.attrib.get("t","")
                    v=c.find(f"{{{MAIN}}}v")
                    raw=v.text or "" if v is not None else ""
                    if t=="s" and raw:raw=strings[int(raw)]
                    elif t=="inlineStr":
                        inline=c.find(f"{{{MAIN}}}is")
                        raw="".join(q.text or "" for q in inline.findall(f".//{{{MAIN}}}t")) if inline is not None else ""
                    vals[i]=raw
                if vals:
                    n=max(vals)+1
                    rows.append([vals.get(j,"") for j in range(n)])
            sheets[name]=rows
    return sheets

def report(workbook,rscript,src):
    sheets=sheet_data(workbook)
    assert "ReadMe" in sheets and "visitation" in sheets and "species visitation" in sheets
    stat={}
    for k,rows in sheets.items():
        nonempty=len(rows)
        # ReadMe is author-public variable dictionary. No individual IDs
        # or raw visitor rows are logged in the source publication summary.
        descriptions=[row[:8] for row in rows[:40]] if k=="ReadMe" else []
        stat[k]={"source_rows_nonempty":nonempty,
                 "first_five_row_field_counts":[len(x) for x in rows[:5]],
                 "header_row":[str(x)[:100] for x in rows[0][:35]] if rows else [],
                 "second_row_preview":[str(x)[:70] for x in rows[1][:25]] if len(rows)>1 else [],
                 "readme_dictionary_rows":descriptions}
    selected=[q.strip()[:180] for q in rscript.splitlines()
              if any(s in q.lower() for s in ("read_excel","sheet =","sheet=","read.xlsx","visitation"))]
    return {"schema":SCHEMA,**src,"original_independent_lit_dark_site_pairs":5,
            "raw_original_sheets":stat,"R_original_sheet_reference_lines":selected[:40],
            "status":"original_workbook_exact_SHA_verified_and_sheet_schema_read",
            "source_plant_visitor_effects_fitted":False,
            "source_visits_cannot_be_treated_as_2160_independent_light_treatments":True,
            "potential_overlap_JPE_rtag160_study_not_resolved":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    workbook,rscript,src=source_open()
    ans=report(workbook,rscript,src)
    file=Path(a.out);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(ans,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("LU2026_SOURCE_EXACT_WORKBOOK_SCHEMA",json.dumps(
      {k:{"n_rows":v["source_rows_nonempty"],"columns":v["header_row"],
          "preview":v["second_row_preview"]}
       for k,v in ans["raw_original_sheets"].items()},ensure_ascii=False,sort_keys=True))
    print("LU2026_README_AUTHOR_EXPLANATIONS",
          json.dumps(ans["raw_original_sheets"]["ReadMe"]["readme_dictionary_rows"],
                     ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
