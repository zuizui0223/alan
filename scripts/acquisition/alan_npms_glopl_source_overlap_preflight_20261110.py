#!/usr/bin/env python3
"""Audit the *real* NPMS 2015-2025 vegetation package and GloPL experiment data.

This is a data-support inspection, not an ALAN→pollen-limitation or
pollen-limitation→plant-cover model. NPMS repeated plots measure vegetation,
NOT measured pollen receipt. GloPL supplementation log response ratios
measure reproductive pollen limitation at their original experiment sites,
NOT a contemporaneous NPMS plot or moth-specific dependency.

The geographical and taxonomic join MUST fail closed when plot-coordinate
or species-name columns are not unambiguously identified. No zero-filled
absence or artificial per-site pollen-limitation values.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd

NPMS_DOI = "10.5285/64caeed0-160f-48f1-b182-cd24a1d7dca8"
GLOPL_DOI = "10.5061/dryad.dt437"
GLOPL_SHA = "264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de"
SCHEMA = "ALAN_NPMS2015_2025_GloPL_source_support_preflight_v1"

def require(condition, message):
    if not condition:
        raise ValueError("ALAN_NPMS_GLOPL: " + message)

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()

def read_csv_member(z, member, *, nrows=None, chunksize=None):
    errors = []
    for enc in ("utf-8-sig", "latin-1"):
        try:
            with z.open(member) as src:
                # Explicit storage schema: do not guess biological meanings.
                return (pd.read_csv(src, encoding=enc, dtype=str,
                                    nrows=nrows, chunksize=chunksize,
                                    low_memory=False), enc)
        except (UnicodeError, pd.errors.ParserError, ValueError) as exc:
            errors.append(f"{enc}:{type(exc).__name__}")
    raise ValueError("Unable to parse original NPMS CSV: "+member+" "+",".join(errors))

def npms_inventory(filename, *, count_rows=True):
    f = Path(filename)
    require(f.exists() and zipfile.is_zipfile(f), "NPMS data are not an actual ZIP archive")
    data = []
    with zipfile.ZipFile(f) as z:
        require(z.testzip() is None, "NPMS ZIP CRC failed")
        info = [x for x in z.infolist() if not x.is_dir()]
        require(len(info) > 0 and sum(x.file_size for x in info) > 10,
                "NPMS source ZIP empty")
        for member in info:
            item = {"name": member.filename, "uncompressed_bytes":member.file_size,
                    "crc32": format(member.CRC, "08x")}
            if member.filename.lower().endswith((".csv", ".tsv")):
                try:
                    sample, encoding = read_csv_member(z,member.filename,nrows=3)
                    item.update(parsed=True, encoding=encoding,
                                columns=list(sample.columns),
                                sample_n_columns=len(sample.columns))
                    if count_rows:
                        block, _ = read_csv_member(z,member.filename,chunksize=75000)
                        item["parsed_rows"] = sum(len(x) for x in block)
                except Exception as exc:
                    item.update(parsed=False, reason=str(exc)[:200])
            data.append(item)
    return {"NPMS_DOI":NPMS_DOI, "NPMS_archive_SHA256":sha256(f),
            "compressed_source_bytes":f.stat().st_size,
            "members":data,"csv_file_count":sum(bool(x.get("parsed")) for x in data),
            "unparsed_CSVs":[x["name"] for x in data if
                              x["name"].lower().endswith(".csv") and not x.get("parsed")]}

def glopl_inventory(csv):
    require(sha256(csv) == GLOPL_SHA,
            "Not the previously pinned 2018 Dryad-identical GloPL original")
    x = pd.read_csv(csv, dtype=str, encoding="latin-1").fillna("")
    require(len(x) == 2969,
            "GloPL original published experimental-case row count changed")
    required = {"Latitude","Longitude","Species_accepted_names","Year","PL_Effect_Size"}
    require(required.issubset(x.columns), "Original accepted GloPL fields are missing")
    coordinates = pd.DataFrame({
        "lat":pd.to_numeric(x.Latitude,errors="coerce"),
        "lon":pd.to_numeric(x.Longitude,errors="coerce")})
    valid = coordinates.lat.between(-90,90)&coordinates.lon.between(-180,180)
    require(valid.all(), "Source GloPL coordinates unexpectedly incomplete")
    species = x.Species_accepted_names.astype(str).str.split().str.join(" ")
    require(species.ne("").all(), "Source GloPL accepted plant names missing")
    # GloPL Year may be PUBLICATION YEAR; do not use as experimental date.
    pl = pd.to_numeric(x.PL_Effect_Size, errors="coerce")
    return {
        "GloPL_DOI": GLOPL_DOI,
        "source_SHA256": GLOPL_SHA,
        "n_experiment_records":len(x),
        "n_species_accepted_names":int(species.nunique()),
        "n_distinct_geographic_points":int(coordinates.drop_duplicates().shape[0]),
        "n_finite_pollen_limitation_log_response_ratio":int(np.isfinite(pl).sum()),
        "columns":list(x.columns),
        "Year_field_NOT_assumed_experiment_year":True,
        "source_is_not_local_moth_pollination_identity":True
    }

def run(npms, glopl, output):
    out = Path(output)
    require(not out.exists(), "Source audit refuses overwrite")
    n = npms_inventory(npms)
    g = glopl_inventory(glopl)
    report = {"schema":SCHEMA,"npms":n,"glopl":g,
        "NPMS_plot_is_NOT_GloPL_experiment_location":True,
        "species_level_PL_transfer_requires_external_validity_check":True,
        "NPMS_vegetation_changes_are_not_fruit_seed_output":True,
        "GloPL_species_lacking_record_have_unknown_PL_not_zero":True,
        "plant_visitor_GloBI_positive_not_effective_moth_dependence":True,
        "validated_NPMS_plot_overlap_count":None,
        "validated_NPMS_GloPL_species_overlap_count":None,
        "ALAN_PL_interaction_regression_fitted":False}
    out.mkdir(parents=True)
    (out/"npms_glopl_actual_source_schema_audit.json").write_text(
        json.dumps(report,sort_keys=True,indent=2)+"\n")
    (out/"npms_package_csv_structure.csv").write_text(
        "file,rows,bytes,columns\n"+"\n".join(
           '"'+f["name"].replace('"','""')+'",'+str(f.get("parsed_rows",""))+
           ","+str(f["uncompressed_bytes"])+',"'+" | ".join(
                f.get("columns",[])).replace('"','""')+'"' for f in n["members"]
            if f["name"].lower().endswith(".csv"))+"\n")
    print("ALAN_NPMS_GLOPL_REAL_SOURCE_SCHEMA_GATE",
          json.dumps({"npms_sha":n["NPMS_archive_SHA256"],
            "npms_csvs":[{"name":f["name"],"rows":f.get("parsed_rows"),
                          "columns":f.get("columns")}
                         for f in n["members"] if f["name"].lower().endswith(".csv")],
            "glopl":{k:g[k] for k in (
                "n_experiment_records","n_species_accepted_names",
                "n_distinct_geographic_points")}},sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("npms","glopl","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.npms,a.glopl,a.output)

if __name__=="__main__":
    main()
