#!/usr/bin/env python3
"""Page et al (2021) ORIGINAL 623 SVE/visitation taxon source: true moth/fly overlap.

Data Page/Nicholson original Figshare project 120921, article 16160265,
file id 30155367, byte SHA fixed. A meaningful comparison requires
SAME study, plant, source year and same SVE measurement metric and units.
Visit metrics must also match for functional q*V comparisons.
No source rows are assigned zeros for undocumented moth/fly visits.
No ALAN lighting treatment; no independent trait-based macro fitness effect.
"""
import argparse,csv,hashlib,io,json,math,statistics,urllib.request
from collections import Counter,defaultdict
from pathlib import Path

SHA="e28ebe99e0dfb27b73fbe6638e568bf98fa95b92c920e1d3d981e8eaea27519a"
SCHEMA="alan_page2021_623_original_taxonlevel_moth_fly_comparator_gate_v1"
HEAD={"studyID","study.year","plant.species","pollinator.group","pollinator.taxon",
      "effectiveness.average","effectiveness.metric",
      "effectiveness.average.units","visit.value","visit.metric"}
QKEY=("studyID","plant.species","study.year",
      "effectiveness.metric","effectiveness.average.units")
VKEY=QKEY+("visit.metric",)
PKEY=("studyID","plant.species","study.year")

def check(cond,msg):
    if not bool(cond):raise RuntimeError("PAGE_ORIGINAL_SOURCE_STOP: "+msg)
def get_source():
    url="https://ndownloader.figshare.com/files/30155367"
    req=urllib.request.Request(url,headers={"User-Agent":"Original Figshare single visit effectiveness archive audit","Accept":"text/csv,*/*"})
    with urllib.request.urlopen(req,timeout=70) as r:blob=r.read(2<<20)
    check(hashlib.sha256(blob).hexdigest()==SHA,"Source 623-row CSV sha changed")
    records=list(csv.DictReader(io.StringIO(blob.decode("utf-8-sig"))))
    check(len(records)==623 and HEAD.issubset(records[0]),"Original data rows/field definitions changed")
    return records

def num(s):
    try:
        n=float(str(s).strip())
        return n if math.isfinite(n) else None
    except (TypeError,ValueError):return None

def compare(rows,key,include_numeric):
    groups=defaultdict(list)
    for r in rows:groups[tuple(r[k].strip() for k in key)].append(r)
    outcomes=[]
    for keys,r in groups.items():
        moth=[a for a in r if a["pollinator.group"].strip()=="moth"]
        flies=[a for a in r if a["pollinator.group"].strip()=="fly"]
        if not moth or not flies:continue
        item={"source_study":keys[0],"source_plant":keys[1],
              "source_year":keys[2],"source_pair_definition":list(key),
              "original_moth_taxon_obs":len(moth),
              "original_fly_taxon_obs":len(flies),
              "original_metrics":sorted(set(z["effectiveness.metric"] for z in moth+flies)),
              "original_units":sorted(set(z["effectiveness.average.units"] for z in moth+flies)),
              "original_visitation_metric":sorted(set(z["visit.metric"] for z in moth+flies)),
              "same_source_study_plant_not_independent_light_experiment":True}
        if include_numeric:
            m=[num(z["effectiveness.average"]) for z in moth]
            f=[num(z["effectiveness.average"]) for z in flies]
            m=[z for z in m if z is not None];f=[z for z in f if z is not None]
            item["source_numeric_moth_n"]=len(m)
            item["source_numeric_fly_n"]=len(f)
            item["source_fly_moth_mean_SVE_ratio"]=(
               statistics.fmean(f)/statistics.fmean(m)
               if m and f and statistics.fmean(m)>0 else None)
            item["ratio_is_observed_SVE_measure_not_necessarily_conspecific_stigma_pollen"]=True
        outcomes.append(item)
    return {"original_source_plant_study_comparison_cells":len(groups),
            "original_moth_fly_matched_cells":len(outcomes),
            "original_studies_with_moth_fly":len(set(z["source_study"] for z in outcomes)),
            "original_plants_with_moth_fly":len(set(z["source_plant"] for z in outcomes)),
            "source_comparisons":outcomes,
            "not_a_valid_precision_weighted_meta_effect":True}

def run(rows):
    bygrp=Counter(r["pollinator.group"].strip() for r in rows)
    bymetric=Counter((r["effectiveness.metric"].strip(),
                      r["effectiveness.average.units"].strip()) for r in rows)
    byvisit=Counter(r["visit.metric"].strip() for r in rows)
    return {
      "schema":SCHEMA,
      "source_article_DOI":"10.6084/m9.figshare.16160265.v1",
      "original_2021_article_DOI":"10.1002/ajb2.1764",
      "author_original_CSV_SHA256":SHA,
      "source_taxonlevel_rows":len(rows),
      "source_study_count":len(set(r["studyID"] for r in rows)),
      "source_plant_count":len(set(r["plant.species"] for r in rows)),
      "source_pollinator_group_counts":dict(bygrp),
      "source_effectiveness_metrics_with_units":{
        f"{k[0]} | {k[1]}":v for k,v in bymetric.items()},
      "source_visit_metrics":dict(byvisit),
      "source_SVE_missing":sum(num(r["effectiveness.average"]) is None for r in rows),
      "source_visit_value_missing":sum(num(r["visit.value"]) is None for r in rows),
      "same_plant_study_year_but_unknown_function_metric_comparability":
          compare(rows,PKEY,False),
      "same_study_plant_year_and_original_effect_metric_unit":
          compare(rows,QKEY,True),
      "STRICT_same_study_plant_year_effect_metric_unit_and_visit_metric":
          compare(rows,VKEY,True),
      "moth_or_fly_not_listed_in_study_NOT_natural_absence":True,
      "source_effect_is_not_always_compatible_stigma_pollen":True,
      "no_ALAN_lighting_manipulation":True,
      "source_623_taxon_entries_not_623_independent_research_sites":True,
      "no_crossstudy_fly_moth_ratio_constructed":True
    }
def main():
    a=argparse.ArgumentParser();a.add_argument("--out",required=True)
    q=a.parse_args()
    data=run(get_source());p=Path(q.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(data,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("PAGE2021_TRUE_ORIGINAL_MOTH_FLY_PLANT_SVE_OVERLAP",
      json.dumps({k:v for k,v in data.items() if k not in (
        "source_effectiveness_metrics_with_units","source_visit_metrics")},
      ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
