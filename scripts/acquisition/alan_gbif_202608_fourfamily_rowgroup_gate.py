#!/usr/bin/env python3
"""Verify genuine individual GBIF occurrence row schema in four Lepidoptera families.

Test only the first original row group of one file per class, not an
ecological or globally representative sample. Never expose raw contributor
coordinates or identifying information in GitHub logs. Source metadata and
data categories are enough to plan exact island point-in-polygon extraction.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pyarrow.parquet as pq
from pyarrow.fs import S3FileSystem
import pandas as pd

FAMILIES=("Noctuidae","Geometridae","Nymphalidae","Pieridae")
DOI="10.15468/dl.8yrbe7"
SCHEMA="ALAN_GBIF_202608_real_fourfamily_row_group_technical_gate_v1"

def require(ok,reason):
    if not ok:raise ValueError("ALAN_GBIF_REAL_ROWS: "+str(reason))

def find_col(columns,n):
    q=next((x for x in columns if x.casefold()==n.casefold()),None)
    require(q is not None,"GBIF source missing "+n)
    return q

def run(manifest,out):
    out=Path(out);require(not out.exists(),"Do not overwrite technical source gate")
    d=json.loads(Path(manifest).read_text())
    require(d["schema"]=="ALAN_GBIF_2026_08_original_lepidoptera_Parquet_frozen_footers_v1"
            and d["published_GBIF_snapshot_DOI"]==DOI,
            "Not original source footer manifest")
    fs=S3FileSystem(anonymous=True,region="us-east-1")
    diag={}
    for fam in FAMILIES:
        choices=[(name,info) for name,info in d["original_partition_footers"].items()
                 if name.endswith("_"+fam)]
        require(len(choices)==1,"Original taxon-specific source partition missing: "+fam)
        name,info=choices[0]
        file=info["original_parquet_files"][0]
        reader=fs.open_input_file("gbif-public-data/"+file["source_key"])
        pqf=pq.ParquetFile(reader)
        cols=pqf.schema_arrow.names
        fields=("gbifid","year","decimallatitude","decimallongitude",
                "basisofrecord","occurrencestatus","family","specieskey",
                "datasetkey")
        actual=[find_col(cols,x) for x in fields]
        tab=pqf.read_row_group(0,columns=actual)
        # Only a small sample needed for technical schema. Do not print rows.
        frame=tab.slice(0,min(10000,tab.num_rows)).to_pandas()
        counts={str(k):int(v) for k,v in frame[find_col(cols,"basisofrecord")].value_counts().head(12).items()}
        yy=pd.to_numeric(frame[find_col(cols,"year")],errors="coerce")
        lat=pd.to_numeric(frame[find_col(cols,"decimallatitude")],errors="coerce")
        lon=pd.to_numeric(frame[find_col(cols,"decimallongitude")],errors="coerce")
        diag[fam]={
            "source_partition":name,
            "original_parquet_key":file["source_key"],
            "original_file_compressed_bytes":file["source_bytes"],
            "source_total_family_original_rows":info["total_original_rows"],
            "source_first_group_row_count":tab.num_rows,
            "tested_first_rows":len(frame),
            "original_required_field_mapping":dict(zip(fields,actual)),
            "source_basis_of_record_sample_counts":counts,
            "sample_valid_2013_2025_year_rows":int(yy.between(2013,2025).sum()),
            "sample_nonmissing_real_GBIF_coordinates":int((lat.between(-90,90)&lon.between(-180,180)).sum()),
            "sample_distinct_source_datasets":int(frame[find_col(cols,"datasetkey")].nunique()),
            "no_original_coordinates_or_observer_ids_redistributed":True,
        }
        print("ALAN_GBIF_FOURFAMILY_TRUE_ROW_GROUP",json.dumps({
            fam:{k:v for k,v in diag[fam].items() if k!="original_required_field_mapping"}
        },sort_keys=True),flush=True)
        reader.close()
    report={"schema":SCHEMA,"GBIF_source_DOI":DOI,
            "actual_family_source_technical_gates":diag,
            "no_global_island_moth_occurrence_statistics_claimed":True,
            "no_unobserved_moth_zero_or_effective_pollination":True}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();run(a.manifest,a.out)
