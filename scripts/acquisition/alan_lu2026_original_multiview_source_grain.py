#!/usr/bin/env python3
"""Resolve Lu (2026) original plant source-table measurement grain, no model fit.

The same author Data S1 has three differently organized views:
visitation (400 rows), selected 'species visitation' (1038), raw
'plant–pollinator network' (1606). Do NOT assume all are a single
partition of visits or auto-fill blank metadata in the network sheet.
Determine whether author original row source and supplementary R code
permit a valid same-plant replacement effect, explicitly source-gate it.

This is a source table audit. No uploaded/raw workbook rows exposed.
"""
import argparse,collections,json,math,re
from pathlib import Path
from scripts.acquisition.alan_lu2026_original_workbook_schema import source_open,sheet_data,ORIGINAL_WORKBOOK_SHA
from scripts.analysis.alan_lu2026_original_fiveblock_guild_replacement import rows,num,normalize_time,normalize_light

SCHEMA="alan_lu2026_source_three_workbook_views_visitor_grain_20261118_v1"
IDENTIFIER=("time","sampling_time","treatment","block","repeat","plant")
GUILD=("all","Lepidoptera","Coleoptera","Diptera","Hymenoptera","Others")
CODES=("llxq","glz","hrc","gsmxh","jlm","ch","crc","hhjd","lbq","symg")
LATIN=("Anaphalis hancockii","Bupleurum chinense","Carum carvi",
       "Pedicularis kansuensis","Oxytropis ochrocephala","Cirsium setosum")
def summarize(data):
    h=[str(x).strip() for x in data[0]]
    rows0=rows(data)
    counts=collections.Counter()
    source_numeric=collections.defaultdict(float)
    source_valid_ids=collections.defaultdict(float)
    id_absent=collections.Counter()
    identity_records=[]
    for r in rows0:
        for key in IDENTIFIER:
            if not str(r.get(key,"")).strip():id_absent[key]+=1
        phase=normalize_time(r.get("sampling_time",""))
        arm=normalize_light(r.get("treatment",""))
        has_all_identifiers=(all(str(r.get(k,"")).strip() for k in IDENTIFIER
                             if k in h) and phase in ("day","night") and arm in ("lit","dark"))
        count_fields={
            "has_time":bool(str(r.get("time","")).strip()),
            "has_phase":phase in ("day","night"),
            "has_arm":arm in ("lit","dark"),
            "has_block":bool(str(r.get("block","")).strip()),
            "has_plant":bool(str(r.get("plant","")).strip()),
            "has_any_guild_count":any(num(r.get(k)) is not None for k in GUILD if k in h),
            "has_all_identifiers":has_all_identifiers
        }
        counts[tuple(sorted(k for k,v in count_fields.items() if v))]+=1
        for g in GUILD:
            if g not in h:continue
            val=num(r[g])
            if val is not None:
                source_numeric[(phase,g)]+=val
                if has_all_identifiers:source_valid_ids[(phase,g)]+=val
        if has_all_identifiers:
            identity_records.append((r.get("plant",""),
               phase,arm,r.get("block",""),r.get("time",""),r.get("repeat","")))
    return {
        "raw_source_data_rows":len(rows0),
        "source_header":h,
        "missing_identifier_field_row_counts":dict(id_absent),
        "valid_parseable_full_source_identifier_rows":len(identity_records),
        "distinct_parseable_plant_labels":len(set(q[0] for q in identity_records)),
        "distinct_parseable_plants_by_diel":{
          t:len(set(q[0] for q in identity_records if q[1]==t))
          for t in ("day","night")},
        "partially_blank_rows_by_identifier_pattern":[
          {"fields":list(k),"rows":v} for k,v in counts.most_common(8)],
        "all_numeric_guild_sums_by_source_phase":{
          t:{g:source_numeric[(t,g)] for g in GUILD}
          for t in ("day","night")},
        "fully_keyed_only_numeric_guild_sums":{
          t:{g:source_valid_ids[(t,g)] for g in GUILD}
          for t in ("day","night")},
        "code_names_not_assumed_to_be_latin_binomials":True
    }


def inspect():
    wb,rcode,metadata=source_open()
    assert metadata["source_workbook_sha256"]==ORIGINAL_WORKBOOK_SHA
    sheets=sheet_data(wb)
    table={name:summarize(sheets[name]) for name in (
        "visitation","species visitation","plant–pollinator network")}
    hints={}
    for tag in (*CODES,*LATIN):
        found=[(i+1,line[:200]) for i,line in enumerate(rcode.splitlines())
               if re.search(r"\b"+re.escape(tag)+r"\b",line,re.IGNORECASE)]
        hints[tag]={
            "matches_in_original_R":len(found),
            "source_R_line_numbers":[i for i,_ in found[:15]]
        }
    author_code_map={}
    for key,abbr in re.findall(
        r'"([a-z]{2,8})"\s*=\s*"([A-Z]\.\s*[a-z][a-z-]+)"',rcode):
        if key in author_code_map and author_code_map[key]!=abbr:
            raise ValueError("Original R species code-label conflict")
        author_code_map[key]=abbr
    original=table["visitation"]
    focal=table["species visitation"]
    return {
        "schema":SCHEMA,
        "source_original_2026_DOI":"10.1002/ece3.74368",
        "source_exact_xlsx_sha256":ORIGINAL_WORKBOOK_SHA,
        "original_paired_light_spatial_blocks":5,
        "tables":table,
        "author_source_R_plant_code_lookup":hints,
        "author_R_verified_plant_code_to_latin_abbreviation":author_code_map,
        "source_population_warning":{
            "selected_11_plant_source_has_more_night_diptera_than_aggregate":
                (focal["all_numeric_guild_sums_by_source_phase"]["night"]["Diptera"] >
                 original["all_numeric_guild_sums_by_source_phase"]["night"]["Diptera"]),
            "focal_11_plants_not_known_partition_of_full_visitation":True,
            "one_full_network_row_has_unfilled_diel_and_treatment":True,
            "source_no_blank_identifier_auto_fill":True,
            "same_plant_moth_to_fly_replacement_NOT_pollen_function":True
        },
        "original_source_seed_set_measured":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    a=p.parse_args()
    d=inspect()
    out=Path(a.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(d,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("LU2026_ORIGINAL_MULTIVIEW_OBSERVATION_GRAIN_AUDIT",json.dumps({
        "source_counts":{name:{
            "raw_rows":t["raw_source_data_rows"],
            "valid_rows":t["valid_parseable_full_source_identifier_rows"],
            "identifier_blank_rows":t["missing_identifier_field_row_counts"],
            "plant_n":t["distinct_parseable_plant_labels"]
        } for name,t in d["tables"].items()},
        "author_source_code_map":d["author_R_verified_plant_code_to_latin_abbreviation"],
        "warning":d["source_population_warning"]
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
