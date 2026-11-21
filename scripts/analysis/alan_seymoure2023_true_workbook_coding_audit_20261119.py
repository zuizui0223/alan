#!/usr/bin/env python3
"""Source-authenticated 2023 ALAN interaction review: study-grain and coding audit.

Never interpret an original spreadsheet row as an independent illumination
experiment. Retain 79 publications versus 81 reported investigations separately
from the original nonempty workbook row count. No p-values / meta effects until
the direction coding, paper identity, sampling and comparators are audited.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import io
import json
import re
import urllib.request
from pathlib import Path

SOURCE_ARTICLE = 24147963
SOURCE_FILE_ID = 42366918
SOURCE_NAME = "rstb20220356_si_002.xlsx"
SOURCE_SHA256 = "a4c105c50f7b53cb2c58b287d065b79432b6ec7c308c35b1c4b873e78a4a359e"
SOURCE_PAPER = "10.1098/rstb.2022.0356"
REVIEW_PUBLISHED_PAPERS = 79
REVIEW_PUBLISHED_INTERACTIONS = 81
SOURCE_URL = f"https://api.figshare.com/v2/articles/{SOURCE_ARTICLE}"
REQUIRED = [
    "interaction_type", "method", "Consumer_taxon", "Resource_taxon",
    "interaction strength (IS) reduction", "IS enhancement", "strength sum",
    "temporal mismatch induced", "spatial mismatch induced", "effect sum",
    "Authors", "Article Title", "Publication Year", "DOI",
]


def get(url):
    request = urllib.request.Request(url, headers={
        "User-Agent": "ALAN-biotic-Figshare-original-grain-audit/2026-11-19",
        "Accept": "*/*",
    })
    with urllib.request.urlopen(request, timeout=50) as r:
        b = r.read(4_000_000)
        if len(b) >= 4_000_000:
            raise ValueError("Unexpectedly large original data source")
        return b


def canonical(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    if isinstance(v, float) and int(v) == v:
        return str(int(v))
    return str(v).strip()


def source_rows(original):
    from openpyxl import load_workbook
    book = load_workbook(io.BytesIO(original), data_only=True, read_only=True)
    if book.sheetnames != ["database"]:
        raise ValueError(f"Unexpected original worksheets: {book.sheetnames}")
    sh = book["database"]
    original_rows=list(sh.iter_rows(values_only=True))
    book.close()
    if not original_rows:
        raise ValueError("Source workbook was empty")
    header = [canonical(x) for x in original_rows[0]]
    if len(header) != len(set(header)):
        raise ValueError("Duplicate original spreadsheet columns")
    if not set(REQUIRED).issubset(header):
        raise ValueError(f"Original mandatory source columns changed {set(REQUIRED) - set(header)}")
    rows = []
    for raw_idx, original_row in enumerate(original_rows[1:], start=2):
        values=list(original_row[:len(header)])
        row = dict(zip(header,[canonical(x) for x in values]))
        if not any(row.values()):
            continue
        row["_original_excel_row"] = raw_idx
        rows.append(row)
    return header, rows


def source_workbook_audit(original):
    hashval=hashlib.sha256(original).hexdigest()
    if hashval != SOURCE_SHA256:
        raise ValueError(f"Source SHA drift. Expected {SOURCE_SHA256}; observed {hashval}")
    header, rows=source_rows(original)
    source_key_map = {
        "interaction_type": "interaction_type",
        "method": "method",
        "Consumer_taxon": "consumer",
        "Resource_taxon": "resource",
        "interaction strength (IS) reduction": "strength_reduction",
        "IS enhancement": "strength_enhancement",
        "strength sum": "strength_sum",
        "temporal mismatch induced": "temporal_mismatch",
        "spatial mismatch induced": "spatial_mismatch",
        "effect sum": "effect_sum",
    }
    status_cols=list(source_key_map)
    paper_cols=["Authors","Article Title","Source Title","Publication Year","DOI","DOI Link"]
    # Keep original codebook for each true field, including blanks, not a
    # researcher-assigned positive or negative score.
    codebook={}
    for col in status_cols:
        counts=Counter(row.get(col,"") for row in rows)
        codebook[source_key_map[col]]={
            "unique_coding_values":len(counts),
            "missing_count":counts.get("",0),
            "observed_original_value_counts":[{"original_value":k or "<BLANK>", "n":v}
                                                for k,v in sorted(counts.items(),key=lambda x:(x[0]=="",str(x[0])))],
        }
    # DOI full URL and unadorned DOI should resolve to a common publication.
    def normdoi(d):
        d=d.casefold().strip()
        d=re.sub(r"^https?://(?:dx\.)?doi\.org/","",d)
        d=re.sub(r"^doi:\s*","",d)
        return d.strip().rstrip(" .,;")
    dois = [normdoi(row.get("DOI","") or row.get("DOI Link","")) for row in rows]
    title_keys=[(row["Authors"].casefold().strip(),row["Article Title"].casefold().strip(),
                 row["Publication Year"]) for row in rows]
    valid_dois=[x for x in dois if x]
    studies=defaultdict(list)
    for i,row in enumerate(rows):
        # Fall back to text title identity if a DOI was not provided. Never
        # interpret a row number as a unique independent scientific paper.
        key="doi:"+dois[i] if dois[i] else "title:"+repr(title_keys[i])
        studies[key].append(row)
    types=Counter(row.get("interaction_type","") for row in rows)
    method_types=Counter(row.get("method","") for row in rows)
    overlap=Counter()
    for key,rr in studies.items():
        overlap[len(rr)]+=1
    other={
        "original_excel_data_nonempty_rows":len(rows),
        "original_header_count":len(header),
        "source_headers_exact":header,
        "published_review_papers":REVIEW_PUBLISHED_PAPERS,
        "published_review_investigations":REVIEW_PUBLISHED_INTERACTIONS,
        "data_row_count_equal_published_investigations":len(rows)==REVIEW_PUBLISHED_INTERACTIONS,
        "source_distinct_DOIs_nonmissing":len(set(valid_dois)),
        "source_DOI_missing_rows":len(rows)-len(valid_dois),
        "source_distinct_publication_IDs_DOI_then_title":len(studies),
        "source_publications_by_number_of_rows":dict(sorted(overlap.items())),
        "source_rows_with_empty_article_title":sum(not r.get("Article Title","") for r in rows),
        "source_rows_with_empty_author":sum(not r.get("Authors","") for r in rows),
        "source_rows_with_empty_year":sum(not r.get("Publication Year","") for r in rows),
        "rows_missing_core_interaction_type":sum(not r.get("interaction_type","") for r in rows),
        "rows_missing_method":sum(not r.get("method","") for r in rows),
        "review_published_81_is_not_source_82_independent_light_experiments":True,
        "years_original_value_counts":dict(sorted(Counter(r.get("Publication Year","") or "<BLANK>" for r in rows).items())),
        "type_original_value_counts":dict(sorted(types.items())),
        "method_original_value_counts":dict(sorted(method_types.items())),
    }
    # This remains a coding audit. Use no single value from a union of studies
    # without verifying whether it represents 1 (yes), 0 (no), or a direction
    # conditional on another recorded stage.
    combos=Counter(
        (r["interaction_type"],r["interaction strength (IS) reduction"],
         r["IS enhancement"],r["strength sum"],r["temporal mismatch induced"],
         r["spatial mismatch induced"],r["effect sum"])
        for r in rows)
    combo_summary=[
        {"interaction_type":k[0],
         "reduction":k[1],"enhancement":k[2],"strength_sum":k[3],
         "temporal_mismatch":k[4],"spatial_mismatch":k[5],"effect_sum":k[6],"n_rows":n}
        for k,n in sorted(combos.items(),key=lambda x:(x[0][0],x[0][1:]))]
    result={
        "schema":"alan_seymoure2023_true_original_79paper_81interaction_source_coding_audit_v1",
        "source_article_DOI":SOURCE_PAPER,
        "original_figshare_article_ID":SOURCE_ARTICLE,
        "original_figshare_file_ID":SOURCE_FILE_ID,
        "original_source_filename":SOURCE_NAME,
        "original_source_SHA256":hashval,
        "source_observed":other,
        "literal_source_codebook":codebook,
        "literal_interaction_direction_combinations":combo_summary,
        "gate":{
            "original_source_verified":True,
            "source_coding_direction_interpreted":False,
            "unique_lighting_installation_metadata_present":False,
            "publication_ID_not_illumination_experiment_ID":True,
            "original_rows_not_81_independent_effect_sizes":True,
            "causal_ALAN_moderator_estimated":False,
            "meta_effect_estimated":False,
            "source_excel_row_level_records_not_exported":True,
            "scope":"Original-publication-level coding audit; no inferred causal light effect on pollen, seeds, survival or predator attack",
        },
    }
    return result


def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--out", type=Path, required=True)
    args=arg.parse_args()
    meta=json.loads(get(SOURCE_URL))
    files=[f for f in meta.get("files",[]) if f.get("id")==SOURCE_FILE_ID and f.get("name")==SOURCE_NAME]
    if len(files)!=1:
        raise ValueError("Source Figshare article does not match expected original file")
    original=get(files[0]["download_url"])
    result=source_workbook_audit(original)
    result["figshare_metadata_retrieval"]={
        "public_article_ID":meta.get("id"),"file_size_reported":files[0].get("size"),
        "actual_byte_length":len(original),"verified_exact_original_binary":True,
    }
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    s=result["source_observed"]
    print("SEYMOURE_2023_TRUE_SOURCE_SHA",result["original_source_SHA256"])
    print("SEYMOURE_2023_TRUE_ROWS_AND_PUBLICATIONS",
          {k:s[k] for k in ("original_excel_data_nonempty_rows",
                            "source_distinct_publication_IDs_DOI_then_title",
                            "source_DOI_missing_rows","rows_missing_core_interaction_type",
                            "rows_missing_method","data_row_count_equal_published_investigations")})
    print("SEYMOURE_2023_TRUE_CODEBOOK",
          {k:v["observed_original_value_counts"] for k,v in result["literal_source_codebook"].items()
           if k in ("strength_reduction","strength_enhancement","strength_sum",
                    "temporal_mismatch","spatial_mismatch","effect_sum")})
    print("SEYMOURE_2023_TRUE_TYPE_COUNTS",s["type_original_value_counts"])


if __name__=="__main__":
    main()
