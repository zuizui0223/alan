#!/usr/bin/env python3
"""Test whether *plant moisture affinity* can be linked to Page pollination redundancy.

Page 2021: 72 original study-plant-period networks from 33 source studies.
Independent Tichy Europe EIV: MOIST exact plant names, original SHA locked.
No ALAN manipulation in Page. Never interpret EIV as in-situ water exposure.
Disallow original-study-confounded wet-vs-dry contrasts without source overlap.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import pandas as pd
from scripts.analysis.alan_page2021_original_623_moth_fly_sameplant import get_source
from scripts.analysis.alan_page2021_replacement_threshold_20261119 import networks, summary

SCHEMA="alan_page2021_72network_independent_EIV_moisture_support_v1"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
def load_moisture(path):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=EIV_SHA:
        raise RuntimeError("OFFICIAL_EIV_SOURCE_HASH_STOP")
    ei=pd.read_csv(path,encoding="latin1",index_col=0,low_memory=False)
    if "MOIST" not in ei.columns:
        raise RuntimeError("OFFICIAL_EIV_MOIST_COLUMN_STOP")
    ei["MOIST"]=pd.to_numeric(ei["MOIST"],errors="coerce")
    if ei.index.isna().any():
        ei=ei[ei.index.notna()]
    # Explicit duplicate consistency, preserving unmatchable names as missing.
    grouped=ei.groupby(level=0).MOIST.agg(["first","nunique"])
    lookup=grouped.loc[grouped.nunique<=1,"first"] if False else grouped.loc[grouped["nunique"]<=1,"first"]
    return lookup.dropna().to_dict()

def audit(nets, moisture):
    output=[]
    for r in nets:
        clean=r["source_plant"].replace("_"," ").strip()
        x=moisture.get(clean)
        x=float(x) if x is not None and 1<=x<=12 else None
        output.append({
           "study":r["source_study"],"plant":clean,
           "MOIST_exact":x,"wet7":x>7 if x is not None else None,
           "functional_threshold_more_visits":r["requires_more_than_equal_visits"],
           "required_multiplier":r["required_new_visits_per_lost_visit"],
           "top_visitor_group":r["source_top_visitor_group"],
           "source_taxa":r["taxa"],
        })
    matches=[r for r in output if r["MOIST_exact"] is not None]
    wet=[r for r in matches if r["wet7"]]
    dry=[r for r in matches if not r["wet7"]]
    bystudy={}
    for r in matches:
        bystudy.setdefault(r["study"],[]).append(r)
    within={k:v for k,v in bystudy.items()
            if any(z["wet7"] for z in v) and any(not z["wet7"] for z in v)}
    # Do not infer a moderator from between-source-study differences.
    enough=(len(within)>=5 and len({r["plant"] for r in wet})>=10 and
            len({r["plant"] for r in dry})>=10 and len(matches)>=30)
    return {
       "schema":SCHEMA,
       "source_page_all_networks":len(nets),
       "source_page_original_studies":len({r["source_study"] for r in nets}),
       "source_page_original_plants":len({r["source_plant"] for r in nets}),
       "matched_study_networks":len(matches),
       "matched_original_studies":len(bystudy),
       "matched_plant_names":len({r["plant"] for r in matches}),
       "matched_wet_networks":len(wet),
       "matched_dry_networks":len(dry),
       "matched_wet_plant_names":len({r["plant"] for r in wet}),
       "matched_dry_plant_names":len({r["plant"] for r in dry}),
       "source_studies_with_both_wet_and_dry":len(within),
       "source_studies_with_both_wet_and_dry_names":sorted(within),
       "wet_functional_threshold_fraction_descriptive":
           sum(z["functional_threshold_more_visits"] for z in wet)/len(wet)
           if wet else None,
       "dry_functional_threshold_fraction_descriptive":
           sum(z["functional_threshold_more_visits"] for z in dry)/len(dry)
           if dry else None,
       "functional_redundancy_moisture_meta_moderator":
           None,
       "status":("WITHIN_STUDY_MODERATOR_POSSIBLE_NEEDS_HIERARCHICAL_ANALYSIS"
                 if enough else "SOURCE_WITHIN_STUDY_MOISTURE_SUPPORT_STOP"),
       "exclude_between_study_ecology_as_moisture_effect":True,
       "no_altered_pollination_by_ALAN_measured":True,
       "missing_EIV_does_not_mean_dry":True,
       "source_EIV_SHA256":EIV_SHA,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--eiv",required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    report=audit(networks(get_source(),True),load_moisture(args.eiv))
    p=Path(args.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("PAGE_EIV_WET_FUNCTIONAL_REDUNDANCY_SUPPORT "+
          json.dumps(report,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
