#!/usr/bin/env python3
"""Exposure of 2023 ALAN interaction review to year and methodology coverage.

No assumption that missing "method" means observational, and no inference
that unpublished study designs are randomized experiments.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import re
from pathlib import Path

from scripts.analysis.alan_seymoure2023_true_workbook_coding_audit_20261119 import (
    SOURCE_SHA256, SOURCE_URL, SOURCE_FILE_ID, SOURCE_NAME, get, source_rows
)
from scripts.analysis.alan_seymoure2023_original_direction_study_units_20261119 import (
    _source_scientific_rows, study_key, summarize_source
)


def profile(rows,strict=True):
    scientific,_=_source_scientific_rows(rows,strict=strict)
    base=summarize_source(rows,strict=strict)
    papers=defaultdict(list)
    for r in scientific:
        papers[study_key(r)].append(r)
    cohorts={"pre2020":[],"2020_2022":[],"unusable_year":[]}
    for k,rs in papers.items():
        years={r.get("Publication Year","") for r in rs}
        if len(years)!=1:
            raise ValueError("Conflicting paper publication years")
        year=next(iter(years))
        if re.fullmatch(r"[0-9]{4}",year):
            y=int(year)
            if not 1960<=y<=2022:
                raise ValueError("Paper falls after source 2022 review or before plausible field")
            name="pre2020" if y<2020 else "2020_2022"
        else:name="unusable_year"
        cohorts[name].append((k,rs))
    result={}
    for cohort,papers_list in cohorts.items():
        study_counter=Counter()
        type_stats=defaultdict(lambda:Counter())
        for k,rr in papers_list:
            anyred=int(any(r["interaction strength (IS) reduction"]=="1" for r in rr))
            anyenh=int(any(r["IS enhancement"]=="1" for r in rr))
            anyspace=int(any(r["spatial mismatch induced"]=="1" for r in rr))
            anytime=int(any(r["temporal mismatch induced"]=="1" for r in rr))
            m=int(any(bool(r.get("method","")) for r in rr))
            study_counter["n_studies"]+=1
            study_counter["coded_reduction"]+=anyred
            study_counter["coded_enhancement"]+=anyenh
            study_counter["coded_any_strength"]+=int(anyred or anyenh)
            study_counter["coded_any_spatial"]+=anyspace
            study_counter["coded_any_temporal"]+=anytime
            study_counter["source_method_documented"]+=m
            for r in rr:
                typ=r["interaction_type"]
                d=type_stats[typ]
                d["n_interactions"]+=1
                d["source_method_recorded"]+=bool(r.get("method",""))
                d["coded_increase_only"]+=int(r["IS enhancement"]=="1" and r["interaction strength (IS) reduction"]=="0")
                d["coded_decrease_only"]+=int(r["interaction strength (IS) reduction"]=="1" and r["IS enhancement"]=="0")
                d["coded_mixed"]+=int(r["interaction strength (IS) reduction"]=="1" and r["IS enhancement"]=="1")
                d["coded_neither"]+=int(r["interaction strength (IS) reduction"]=="0" and r["IS enhancement"]=="0")
                d["coded_spatial"]+=int(r["spatial mismatch induced"]=="1")
                d["coded_temporal"]+=int(r["temporal mismatch induced"]=="1")
        result[cohort]={
            "publication_counts":dict(study_counter),
            "interaction_type_counts":{k:dict(v) for k,v in sorted(type_stats.items())}
        }

    # Directional contrast is *composition-dependent*. Exhibit sensitivity
    # across two publication periods without assuming year is randomized.
    year_groups=[]
    for name in ("pre2020","2020_2022"):
        q=result[name]["interaction_type_counts"]
        p=q.get("predator-prey",{})
        h=q.get("herbivory",{})
        def signed(o):
            n=o.get("n_interactions",0)
            return ((o.get("coded_increase_only",0)-o.get("coded_decrease_only",0))/n
                    if n else None)
        delta=signed(p)-signed(h) if signed(p) is not None and signed(h) is not None else None
        year_groups.append({
            "publication_bin":name,
            "predator_interaction_n":p.get("n_interactions",0),
            "herbivory_interaction_n":h.get("n_interactions",0),
            "coded_signed_difference_increase_minus_decrease_per_all_interactions":delta,
            "non_causal_posthoc_year_stratification":True,
            "can_fit_within_period_moderator":False,
        })
    all_papers=len(papers)
    if strict and all_papers!=79:
        raise ValueError("Original 79 publication count drift")
    return {
        "schema":"alan_seymoure2023_original_79paper_publication_year_methods_caution_v1",
        "source_SHA256":SOURCE_SHA256,
        "source_publication_count":all_papers,
        "year_stratified_source_coding":result,
        "year_specific_guild_diagnostics":year_groups,
        "original_method_observed_studies":all_papers-base["publication_level_DOI_aggregated_counts"]["source_studies_all_methods_blank"],
        "original_method_missing_studies":base["publication_level_DOI_aggregated_counts"]["source_studies_all_methods_blank"],
        "source_caveats":{
            "method_blank_is_not_natural_light_control":True,
            "method_present_does_not_prove_random_light_assignment":True,
            "multiple_yearly_observations_not_from_original":True,
            "published_year_is_not_study_year":True,
            "four_binary_interaction_flags_not_comparable_effect_sizes":True,
            "selected_literature_not_probability_sample":True,
            "no_reweighting_to_all_global_species_interactions":True,
            "no_original_study_records_redistributed":True,
        }
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    arg=p.parse_args()
    meta=json.loads(get(SOURCE_URL))
    target=[f for f in meta.get("files",[]) if f.get("id")==SOURCE_FILE_ID and f.get("name")==SOURCE_NAME]
    if len(target)!=1:raise ValueError("Original Figshare metadata changed")
    data=get(target[0]["download_url"])
    if hashlib.sha256(data).hexdigest()!=SOURCE_SHA256:
        raise ValueError("Original source binary hash drift")
    hdr,rows=source_rows(data)
    if len(hdr)!=20:
        raise ValueError("Source variable schema drift")
    ans=profile(rows,strict=True)
    arg.out.parent.mkdir(exist_ok=True,parents=True)
    arg.out.write_text(json.dumps(ans,ensure_ascii=False,indent=2)+"\n")
    print("ALAN_YEAR_STUDY_SOURCE_COHORT",ans["year_stratified_source_coding"])
    print("ALAN_YEAR_REVIEW_METHODOLOGY_GATES",ans["year_specific_guild_diagnostics"])

if __name__=="__main__":
    main()
