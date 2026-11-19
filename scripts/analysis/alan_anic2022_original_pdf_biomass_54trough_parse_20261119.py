#!/usr/bin/env python3
"""Reconstruct author-original Anic 2022 Figshare numeric biomass table from PDF.

Important provenance: one public 24-page PDF contains *actual mesocosm
rows*, not merely a research paper. We deliberately read only the biomass
section before the carbon/nitrogen section, validate an experimental 3×3×6
light×trophic randomization, and export ONLY counts, total=root+shoot QA,
species identities and aggregate treatment cells. Never export plot IDs,
individual biomass measurements, species-by-plot records or full raw tables.

Original Anic source DOI 10.1002/ece3.9157,
Figshare article 20304039/file36265164,
SHA256 baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0.

No multivariate mediation; original 2018 aphid/slug CSV remain source-stopped.
"""
import argparse
from collections import Counter,defaultdict
import hashlib
import io
import json
import math
import re
from pathlib import Path

from scripts.acquisition.alan_anic2022_figshare_original_pdf_table_layout_preflight_20261119 import (
    FIGSHARE_API,FILE_ID,FILE_SHA256,get,
)

BIOMASS_SPECIES=(
    "Achillea millefolium", "Agrostis tenuis",
    "Anthoxanthum odoratum","Centaurea nigra","Cynosurus cristatus",
    "Dactylis glomerata","Festuca ovina",
    "Holcus lanatus","Hypochaeris radicata",
    "Leontodon saxatilis","Lotus corniculatus","Lotus pedunculatus",
    "Plantago lanceolata","Prunella vulgaris",
    "Trifolium dubium","Trifolium pratense",
    "Ranunculus acris","Leucanthemum vulgare",
)
# The paper reports 13 species with sufficient survivors for analysis;
# preserve original 18-species lookup for parser correctness.
SORTED_SPECIES=sorted(BIOMASS_SPECIES,key=len,reverse=True)
EXPECTED_ORIGINAL_PAGE_COUNT=24
BIONUM=re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
OUTPUT_SCHEMA="alan_anic2022_original_54trough_PDF_biomass_parsed_counts_and_cell_means_v1"


def identify_species_prefix(line):
    st=line.strip()
    for sp in SORTED_SPECIES:
        if st.startswith(sp) and (len(st)==len(sp) or st[len(sp)].isspace()):
            return sp,st[len(sp):].split()
    return None,[]


def parse_table_line(line):
    sp,cells=identify_species_prefix(line)
    if sp is None:
        return None
    if len(cells)!=6:
        raise ValueError("Source species row with unexpected mesocosm/light/trophic/biomass columns")
    meso,light,trophic,total,root,shoot=cells
    if not all(BIONUM.fullmatch(v) for v in (total,root,shoot)):
        raise ValueError("Source original three biomass measurement columns nonnumeric")
    vals=[float(total),float(root),float(shoot)]
    if any(not math.isfinite(v) or v<0 or v>1000 for v in vals):
        raise ValueError("Source original dry biomass values impossible")
    if not meso or not light or not trophic:
        raise ValueError("Original mesocosm allocation missing")
    # Prefer code identifiers with no whitespace; never print exact mesocosm.
    if len(meso)>40 or len(light)>30 or len(trophic)>30:
        raise ValueError("Original source treatment or plot code unexpectedly long")
    return {"species":sp,"meso":meso,"light":light,"trophic":trophic,
            "total":vals[0],"root":vals[1],"shoot":vals[2]}


def read_original_biomass_pdf(source,strict=True):
    from pypdf import PdfReader
    if hashlib.sha256(source).hexdigest()!=FILE_SHA256:
        raise ValueError("Authentic original PDF source SHA mismatch")
    if not source.startswith(b"%PDF-"):
        raise ValueError("Original data not a PDF")
    book=PdfReader(io.BytesIO(source))
    if strict and len(book.pages)!=EXPECTED_ORIGINAL_PAGE_COUNT:
        raise ValueError("Original 24 source pages changed")
    parsed=[]
    unparsed_with_species=[]
    for page_index,page in enumerate(book.pages):
        text=page.extract_text() or ""
        # The FIRST 14 pages contain biomass species, mesocosm, L, T,
        # Total(g), root(g), shoot(g); the next 10 pages contain foliar C/N.
        # Do not parse repeated species names in the nutrient table.
        if page_index>=14:
            continue
        for line in text.splitlines():
            sp,parts=identify_species_prefix(line)
            if sp is None:
                continue
            try:
                row=parse_table_line(line)
                if row:parsed.append(row)
            except ValueError as ex:
                unparsed_with_species.append({
                    "page_index":page_index,
                    "species":sp,
                    "token_count_after_species":len(parts),
                    "exception":type(ex).__name__,
                    # Original numeric values NEVER exported.
                })
    results={
        "schema":OUTPUT_SCHEMA,
        "source_original_figshare_sha256":FILE_SHA256,
        "original_total_PDF_pages":len(book.pages),
        "original_biomass_section_pages_read":14,
        "original_species_rows_parsed":len(parsed),
        "original_species_record_types":sorted(set(x["species"] for x in parsed)),
        "original_species_rows_unparsed_metadata_only":unparsed_with_species[:30],
        "original_species_rows_unparsed_total":len(unparsed_with_species),
        "raw_plot_IDs_or_individual_measurements_exported":False,
        "original_2018_Dryad_herbivory_table_not_joined":True,
        "original_2022_PDF_is_not_a_second_randomized_lighting_experiment":True,
        "new_cross_trophic_plant_arthropod_mediation_fitted":False,
    }
    if not parsed:
        results["source_status"]="SOURCE_STOP_NO_BIOMASS_TABLE_ROWS_PARSED"
        results["experimental_light_trophic_full_factorial_verified"]=False
        return results
    ids={}
    conflicts=0
    dup=0
    perSpecies=set()
    cellvalues=defaultdict(list)
    perplot=defaultdict(set)
    total_residual=[]
    lightcounts=Counter()
    trophiccounts=Counter()
    for row in parsed:
        pid=row["meso"]
        allocation=(row["light"],row["trophic"])
        if pid in ids and ids[pid]!=allocation:
            conflicts+=1
        ids.setdefault(pid,allocation)
        key=(row["species"],pid)
        if key in perSpecies:
            dup+=1
        perSpecies.add(key)
        perplot[pid].add(row["species"])
        cellvalues[(row["species"],row["light"],row["trophic"])].append(row["total"])
        total_residual.append(abs(row["total"]-row["root"]-row["shoot"]))
        lightcounts[row["light"]]+=1
        trophiccounts[row["trophic"]]+=1
    assignment=Counter(ids.values())
    # Original IDs are not exported. Light/trophic category names, counts
    # and species×cell aggregates are public methods information only.
    assignment_details=[{"original_light_code":light,"original_trophic_code":trophic,
                         "unique_physical_troughs":num}
                        for (light,trophic),num in sorted(assignment.items())]
    values_out=[]
    for (species,light,trophic),values in sorted(cellvalues.items()):
        values_out.append({
            "species":species,"original_light_code":light,
            "original_trophic_code":trophic,
            "original_mesocosms_with_this_species":len(values),
            "original_biomass_total_g_arithmetic_mean":sum(values)/len(values),
            "original_biomass_total_g_min":min(values),
            "original_biomass_total_g_max":max(values),
            "not_an_original_author_adjusted_effect_or_SE":True,
            "not_a_causal_mediation_estimate":True,
        })
    source_status="AUTHENTIC_ORIGINAL_BIOMASS_SOURCE_PARSED"
    fully_verified=(len(ids)==54 and conflicts==0 and dup==0 and
                    len(set(x["light"] for x in parsed))==3 and
                    len(set(x["trophic"] for x in parsed))==3 and
                    len(assignment)==9 and
                    set(assignment.values())=={6} and
                    not unparsed_with_species and
                    max(total_residual)<0.51)
    if not fully_verified:
        source_status="PARTIAL_OR_FAILED_ORIGINAL_BIOMASS_FACTORIAL_GATE"
    results.update({
        "source_status":source_status,
        "original_unique_physical_mesocosms_in_biomass_table":len(ids),
        "original_light_code_count":len(set(x["light"] for x in parsed)),
        "original_trophic_code_count":len(set(x["trophic"] for x in parsed)),
        "original_factorial_cell_assignments_counts_only":assignment_details,
        "source_original_treatment_conflicts":conflicts,
        "source_original_species_trough_duplicate_rows":dup,
        "original_total_root_shoot_residual_max_g":max(total_residual),
        "original_species_count":len(set(x["species"] for x in parsed)),
        "source_original_species_by_light_trophic_cell_descriptive_means":values_out,
        "experimental_light_trophic_full_factorial_verified":fully_verified,
        "source_plant_biomass_no_original_raw_rows_exported":True,
        "species_conditioned_study_aggregate_not_reproductive_fitness":True,
        "source_2018_aphid_slug_plant_join_verified":False,
    })
    return results


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    meta=json.loads(get(FIGSHARE_API))
    files=[f for f in meta.get("files",[]) if f.get("id")==FILE_ID]
    if len(files)!=1:
        raise ValueError("Exact author Figshare file ID missing")
    original=get(files[0]["download_url"])
    result=read_original_biomass_pdf(original,strict=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("ANIC2022_ORIGINAL_PDF_PLANT_ROWS",result["original_species_rows_parsed"])
    print("ANIC2022_ORIGINAL_PDF_SOURCE_STATUS",result["source_status"])
    print("ANIC2022_ORIGINAL_BIOMASS_SOURCE_SPECIES",result["original_species_record_types"])
    print("ANIC2022_ORIGINAL_BIOMASS_FACTORIAL_ID_COUNTS",
          result.get("original_factorial_cell_assignments_counts_only",[]))
    print("ANIC2022_ORIGINAL_BIOMASS_FAILED_SOURCE_ROW_SHAPES",
          result["original_species_rows_unparsed_metadata_only"][:10])
    print("ANIC2022_54_TROUGH_SOURCE_REANALYSIS_GATE",
          result.get("experimental_light_trophic_full_factorial_verified",False))

if __name__=="__main__":
    main()
