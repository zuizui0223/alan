#!/usr/bin/env python3
"""Original-source 2018 Penryn 54 mesocosms: genuine 13-CSV experimental grain.

Directly downloads the official Dryad dataset ZIP, pins SHA256 and 13
individual source CSV checksums; performs an ephemeral in-memory join across
54 original trough IDs, 3 lighting and 3 trophic conditions. Publishes
only source-grain checks and treatment-cell aggregates (no raw plot IDs or
original author rows).

Original source DID NOT include 'slugs2015.csv' in the 2018 author's
mollusc treatment-effect fit: 2015 follows 2014 decline & 2015
reintroduction. Do not include 2015 slugs in any inferential analysis.
Longitudinal counts are measured repeatedly within plot, not more
independent lighting experiments.

Source: 10.5061/dryad.8217n76, version id 15138, official dataset ZIP
SHA256 5dc5bc93e48758364d068c113edb84e05b81f99b22827c4c6ee0c14f42ba3130.
"""
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

SOURCE_DOI="10.5061/dryad.8217n76"
ARCHIVE_SHA256="5dc5bc93e48758364d068c113edb84e05b81f99b22827c4c6ee0c14f42ba3130"
SOURCE_VERSION_ID=15138
ARCHIVE_URL="https://datadryad.org/api/v2/datasets/"+urllib.parse.quote("doi:"+SOURCE_DOI,safe="")+"/download"
MANIFEST={
"CNratio.csv":"62ea01bc0983b64c98fe243f50586a96f7ebbccc897bd2b34d13296bcfb04a1e",
"point_quadrat_110713.csv":"347be45e1bfd7bffb522a99a687cc4f1a27ac83a7787d5b7451e2d5832622494",
"point_quadrat_240714.csv":"d71ff4fbb0dd95e4e297495fd5b5fb4f0c7e4c75401f54442ac34d2b749b41dd",
"point_quadrat_220715.csv":"5af785c16415ffdd75a371153ad34defa47f9d6f2512031b072ff2d4c7bf05f1",
"lotus_p_2013.csv":"7941b16f4e8468fdaaabb95887f7749264a10b7429976518ca200109b94aa05e",
"lotus_p_2014.csv":"a07d1b553febb6456e8d603b80096822690502b3ab7c55c2415f66c1cadda13c",
"lotus_p_2015.csv":"e26c3169856bbe682b9eb3a9eaa6f67fb250490bb4664ef1f6000787c688876a",
"PeaAphid2013.csv":"0a59187cae45850f1285ff4c1d29cc75de77055b00974b09e72e20f8d174c6a0",
"PeaAphid2014.csv":"901b5663f5dc4707f0d79eaab45c566abb94c990e6fc649235882899c0f8abe3",
"PeaAphid2015.csv":"e2933e2a36ee234d96e56029920e0efd445f2b6a71ee298bfb76a33e45573eb4",
"slugs2013.csv":"9604e932c701daa8c832ced68ea5d49f91309fb5ca2dc8874cd5b33da4a07f3e",
"slugs2014.csv":"8c4fab2456d9725117e0438f371adfa921d82f460e6dd4e33f45acfc47ba3f31",
"slugs2015.csv":"3333846bc09d475dfd84a6b2268e7c21701ede3f70e9e05d73b9787a918e2da6",
}
OUT_SCHEMA="alan_2018_original_13CSV_54trough_repeated_physical_factor_and_outcome_gate_v1"
SOURCE_ROWS={"CNratio.csv":51,"slugs2013.csv":36,"slugs2014.csv":36,"slugs2015.csv":36}
for group in ("lotus_p_2013.csv","lotus_p_2014.csv","lotus_p_2015.csv",
              "PeaAphid2013.csv","PeaAphid2014.csv","PeaAphid2015.csv",
              "point_quadrat_110713.csv","point_quadrat_240714.csv","point_quadrat_220715.csv"):
    SOURCE_ROWS[group]=54

ID_ALIASES={"","plot","mesocosm"}
LIGHT_ALIASES={"l treatment","l.treatment","l_treatment"}
TROPH_ALIASES={"t treatment","t.treatment","t_treatment"}
SOURCE_SLUG_FIT_YEARS=(2013,2014)


def get(url):
    """Use the exact official Dryad archive downloader already proven in CI.

    The first independent workflow transiently returned HTTP401 because it
    used different request headers. Do not fall back to untrusted mirrors or
    assume an HTTP401 means original data are scientifically missing.
    """
    from scripts.acquisition.alan_penryn54_crosspaper_figshare_dryad_primary_audit_20261119 import pull
    return pull(url,limit=2_000_000)


def parse_zip_archive(data,strict=True):
    if strict and hashlib.sha256(data).hexdigest()!=ARCHIVE_SHA256:
        raise ValueError("Author-original Dryad archive SHA256 drift")
    if not data.startswith(b"PK\x03\x04"):
        raise ValueError("Original source was not ZIP")
    tables={}
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        members={Path(p).name:p for p in z.namelist() if p.lower().endswith(".csv")}
        if strict and set(members)!=set(MANIFEST):
            raise ValueError("Expected 13 source CSVs and no extras")
        for name,p in members.items():
            raw=z.read(p)
            if strict and hashlib.sha256(raw).hexdigest()!=MANIFEST[name]:
                raise ValueError(f"Original source CSV SHA drift {name}")
            content=raw.decode("utf-8-sig",errors="replace")
            reader=csv.DictReader(io.StringIO(content))
            if not reader.fieldnames:
                raise ValueError(f"Source header empty {name}")
            if len(set(reader.fieldnames))!=len(reader.fieldnames):
                raise ValueError("Original field names collide")
            rows=list(reader)
            if strict and len(rows)!=SOURCE_ROWS[name]:
                raise ValueError(f"Source experimental row count changed {name}={len(rows)}")
            tables[name]={"fields":reader.fieldnames,"rows":rows}
    return tables


def id_norm(v):
    v=(v or "").strip().upper()
    if re.fullmatch(r"\d+(?:\.0+)?",v):
        v=str(int(float(v)))
    if not v:raise ValueError("Original mesocosm ID is missing")
    return v


def source_factor_field(fields,alias_set,name):
    opts=[f for f in fields if (f or "").strip().lower() in alias_set]
    if len(opts)!=1:raise ValueError(f"Cannot identify original {name} field")
    return opts[0]


def table_plot_factor_grain(tables,strict=True):
    out={}
    all_allocations=defaultdict(set)
    categorical=defaultdict(set)
    data_per_file={}
    source_legacy_name_warnings=[]
    for name,data in sorted(tables.items()):
        fields=data["fields"]
        lightkey=source_factor_field(fields,LIGHT_ALIASES,"light")
        trophkey=source_factor_field(fields,TROPH_ALIASES,"trophic")
        id_opts=[f for f in fields if (f or "").strip().lower() in ID_ALIASES]
        if name=="slugs2015.csv":
            # Column "Date" is a legacy CSV artefact and NOT an id/year;
            # plot is explicitly stored second.
            id_opts=[f for f in id_opts if (f or "").lower()=="plot"]
        if len(id_opts)!=1:raise ValueError(f"Missing physical plot ID in {name}")
        idkey=id_opts[0]
        seen=set()
        table_alloc={}
        for row in data["rows"]:
            pid=id_norm(row.get(idkey))
            if pid in seen:raise ValueError(f"Duplicated physical mesocosm row in {name}")
            seen.add(pid)
            l=(row.get(lightkey) or "").strip()
            t=(row.get(trophkey) or "").strip()
            if not l or not t:raise ValueError("Original light×trophic labels missing")
            table_alloc[pid]=(l,t)
            all_allocations[pid].add((l,t))
            categorical["light"].add(l)
            categorical["trophic"].add(t)
        out[name]={
            "original_physical_unit_count":len(seen),
            "original_header_row_count":len(data["rows"]),
            "source_id_column_name":idkey or "<ORIGINAL_UNNAMED_FIRST_COLUMN>",
            "original_light_column_name":lightkey,
            "original_trophic_column_name":trophkey,
            "n_source_light_codes":len({v[0] for v in table_alloc.values()}),
            "n_source_trophic_codes":len({v[1] for v in table_alloc.values()}),
            "original_unique_sampled_factorial_cells":len(set(table_alloc.values())),
            "no_original_plot_ids_exported":True,
        }
        data_per_file[name]=table_alloc
    global_consistent=sum(len(x)==1 for x in all_allocations.values())
    conflicts=sum(len(x)>1 for x in all_allocations.values())
    observed_cohort={pid:next(iter(v)) for pid,v in all_allocations.items() if len(v)==1}
    cells=Counter(observed_cohort.values())
    # If the same original physical mesocosm key is reused differently
    # between files, refuse to make cross-stage statements.
    verified=(len(observed_cohort)==54 and conflicts==0 and len(cells)==9 and
              set(cells.values())=={6} and
              len(categorical["light"])==3 and len(categorical["trophic"])==3)
    if strict and not verified:
        # Do not fail source fetching, but do not claim joint experiment.
        source_status="FACTORIAL_SOURCE_STOP_OR_SOURCE_CODE_INCONSISTENCY"
    else:source_status="ORIGINAL_54_SOURCE_FACTORIAL_PLOT_JOIN_VERIFIED"
    return {
       "status":source_status,
       "original_unique_physical_IDs_across_13_source_CSVs":len(all_allocations),
       "source_same_plot_different_treatment_code_conflict_count":conflicts,
       "cross_file_exact_original_plot_factor_allocations_consistent":verified,
       "original_source_factor_codes":{
         "light":sorted(categorical["light"]),"trophic":sorted(categorical["trophic"])
       },
       "all_original_physical_light_trophic_cell_counts_only":[
         {"light_original_code":l,"trophic_original_code":t,
          "unique_original_troughs":n}
         for (l,t),n in sorted(cells.items())
       ],
       "original_13_source_CSV_physical_grain":out,
       "source_plot_identity_recovered_in_memory":True,
       "source_original_plot_IDs_committed_or_exported":False,
       "source_original_biological_rows_exported":False,
       "eligible_original_slug_years_for_authors_2018_model":[2013,2014],
       "post_reintroduction_slugs_2015_excluded_from_authors_original_fit":True,
       "no_source_2016_biomass_join_yet":True,
       "no_new_causal_medation_or_light_x_trophic_glmm_fitted":True,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    data=get(ARCHIVE_URL)
    tables=parse_zip_archive(data,strict=True)
    report=table_plot_factor_grain(tables,strict=True)
    report.update({
        "schema":OUT_SCHEMA,
        "original_dataset_DOI":SOURCE_DOI,
        "original_dataset_version_id":SOURCE_VERSION_ID,
        "original_zip_sha256":ARCHIVE_SHA256,
        "source_all_13_source_member_hashes_verified":True,
        "source_all_13_source_row_counts_verified":True,
        "physical_light_experiment_programmes":1,
        "original_54_3x3_troughs_if_consistent":54,
    })
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_PENRYN54_ORIGINAL_SOURCE_FACTOR_CODES",report["original_source_factor_codes"])
    print("ALAN_PENRYN54_ORIGINAL_SOURCE_FULL_54_PLOT_ASSIGNMENT_STATUS",report["status"])
    print("ALAN_PENRYN54_ORIGINAL_SOURCE_PHYSICAL_CELL_COUNTS",report["all_original_physical_light_trophic_cell_counts_only"])
    print("ALAN_PENRYN54_ORIGINAL_SOURCE_13_CSV_FILE_GRAIN",report["original_13_source_CSV_physical_grain"])
    print("ALAN_PENRYN54_SOURCE_FIT_STATUS",report["no_new_causal_medation_or_light_x_trophic_glmm_fitted"])


if __name__=="__main__":
    main()
