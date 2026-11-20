#!/usr/bin/env python3
"""Cross-species original vial-level reanalysis of Campos et al. 2026 ALAN feeding.

Original source: https://doi.org/10.5281/zenodo.19022243
Original Excel SHA256 b5c36d9b4c184a3f563388d56a9072fcfe9fa8c988600e9979054402083a7d6d
3 independent species-specific source feeding panels (not 3 studies):
- Potamopyrgus antipodarum 0/1/10 lux, 7 days, vial=10 snails
- Sericostoma vittatum 0/1/10 lux, 7 days, vial=1 larva
- Chironomus riparius 0/1/10 lux, 14 days, vial=10 larvae

Primary descriptive estimand: log(mean feeding at lux / mean feeding at 0).
Resample the whole source vial, not individual organisms or source rows.
Cross-species log ratio contrast at same ALAN lux is dimensionless and
still confounded by exposure duration and normalization differences.
No effect on live wet plants, their reproduction or selection is estimated.
"""
from __future__ import annotations
import argparse, collections, hashlib, json, math, random, statistics, urllib.request
from pathlib import Path
from scripts.acquisition.alan_lu2026_original_workbook_schema import sheet_data

SCHEMA="alan_Campos2026_3source_species_vial_feeding_logresponse_v1"
ORIGINAL_XLSX="b5c36d9b4c184a3f563388d56a9072fcfe9fa8c988600e9979054402083a7d6d"
URL="https://zenodo.org/api/records/19022243/files/data_Campos%20et%20al%202026.xlsx/content"
SEED=20261120
BOOT=4999
SPECIES={
 "feeding rates P. antipodarum":{
  "species":"Potamopyrgus antipodarum","vials_expected_per_lux":10,
  "duration_days":7,"animals_per_vial":10,
  "original_response_unit":"mg conditioned alder leaf consumed / number of surviving organisms / day"},
 "feeding rates S. vittatum":{
  "species":"Sericostoma vittatum","vials_expected_per_lux":20,
  "duration_days":7,"animals_per_vial":1,
  "original_response_unit":"mg conditioned alder leaf consumed / mg organism wet mass / day"},
 "feeding rates C. riparius":{
  "species":"Chironomus riparius","vials_expected_per_lux":15,
  "duration_days":14,"animals_per_vial":10,
  "original_response_unit":"mg conditioned alder leaf consumed / number of surviving organisms / day"},
}

def check(ok,why):
    if not ok:raise ValueError("CAMPOS2026_ORIGINAL_SOURCE_STOP: "+why)

def fetch_original():
    with urllib.request.urlopen(urllib.request.Request(
          URL,headers={"User-Agent":"ALAN-vial-source-comparative-audit"}),timeout=60) as u:
        raw=u.read(2<<20)
    check(hashlib.sha256(raw).hexdigest()==ORIGINAL_XLSX,
          "official Zenodo original spreadsheet SHA changed")
    return sheet_data(raw)

def mean(vals):return statistics.fmean(vals)

def read_group(sheet, spec):
    check(len(sheet)>=3,"no source rows")
    header=[str(x) for x in sheet[0]]
    check(header[0].strip()=="ALAN" and header[1].strip()=="replicate",
          "light/replicate source headers changed")
    check("feeding" in header[2].lower(),"third column not a feeding endpoint")
    group=collections.defaultdict(list)
    ids=collections.defaultdict(set)
    missing=collections.Counter()
    for row in sheet[1:]:
        try:
            al=float(row[0])
            rep=int(float(row[1]))
        except (ValueError,IndexError,TypeError):continue
        check(al in (0.,1.,10.),"unexpected illumination intensity")
        check(rep not in ids[al],"duplicated vial replicate ID under same light")
        ids[al].add(rep)
        try:rate=float(row[2])
        except (ValueError,IndexError,TypeError):
            missing[al]+=1
            continue
        check(math.isfinite(rate) and rate>=0.,"invalid source feeding rate")
        group[al].append(rate)
    require=spec["vials_expected_per_lux"]
    for al in (0.,1.,10.):
        check(len(ids[al])==require,
              f"source exposure group source vial sample changed {spec['species']} lux={al} {len(ids[al])} != {require}")
        check(len(group[al])>=require//2,"too few measured original vials")
        check(mean(group[al])>0,"nonpositive source response")
    return group,{str(k):{"original_vials":len(ids[k]),"observed_feeding_vials":len(group[k]),
                          "missing_feeding_vials":missing[k]}
                 for k in (0.,1.,10.)}

def boot_logresponse(group,rng,reps=BOOT):
    raw_means={int(k):mean(v) for k,v in group.items()}
    points={str(lux):math.log(raw_means[lux]/raw_means[0]) for lux in (1,10)}
    draws={str(lux):[] for lux in (1,10)}
    for _ in range(reps):
        gm={int(lux):mean(rng.choices(vals,k=len(vals))) for lux,vals in group.items()}
        if gm[0]<=0 or gm[1]<=0 or gm[10]<=0:
            raise RuntimeError("ALL_ZERO_BOOTSTRAP_CELL_STOP")
        for lux in (1,10):
            draws[str(lux)].append(math.log(gm[lux]/gm[0]))
    intervals={key:[float(sorted(v)[int(.025*(reps-1))]),
                    float(sorted(v)[int(.975*(reps-1))])] for key,v in draws.items()}
    stats={"mean_feed_by_lux":{str(lux):raw_means[lux] for lux in (0,1,10)},
           "log_ratio_to_zero_by_lux":points,
           "ratio_to_zero_by_lux":{k:math.exp(v) for k,v in points.items()},
           "vial_bootstrap_95pct_log_ratio":intervals,
           "vial_bootstrap_95pct_ratio":{k:[math.exp(x) for x in v] for k,v in intervals.items()},
           "source_bootstrap_reps":reps}
    return stats,draws

def percentile(vec):
    x=sorted(vec);n=len(x)
    return [x[int(.025*(n-1))],x[int(.975*(n-1))]]

def analyze(sheets):
    check(len(sheets)==6,"original Zenodo six sheet change")
    summary={}
    draws={}
    for i,(sheet,spec) in enumerate(SPECIES.items()):
        check(sheet in sheets,"required source feeding sheet missing")
        groups,counts=read_group(sheets[sheet],spec)
        stats,sims=boot_logresponse(groups,random.Random(SEED+1000*i),BOOT)
        summary[spec["species"]]={"source_sheet":sheet,**spec,
             "vial_cohort":counts,**stats}
        draws[spec["species"]]=sims
    cross={}
    pairings=[
      ("Potamopyrgus antipodarum","Chironomus riparius"),
      ("Sericostoma vittatum","Chironomus riparius"),
      ("Potamopyrgus antipodarum","Sericostoma vittatum")
    ]
    for a,b in pairings:
        out={}
        for lux in ("1","10"):
            point=(summary[a]["log_ratio_to_zero_by_lux"][lux]-
                   summary[b]["log_ratio_to_zero_by_lux"][lux])
            differences=[x-y for x,y in zip(draws[a][lux],draws[b][lux])]
            interval=percentile(differences)
            out[lux]={"log_ratio_of_response_ratios":point,
                 "ratio_of_response_ratios":math.exp(point),
                 "vial_bootstrap_95pct_log_ratio":interval,
                 "vial_bootstrap_95pct_ratio":[math.exp(x) for x in interval],
                 "contrast_exposure_time_differs_between_species":
                   summary[a]["duration_days"]!=summary[b]["duration_days"]}
        cross[a+"__minus__"+b]=out
    return {"schema":SCHEMA,"original_zenodo_source_doi":"10.5281/zenodo.19022243",
         "original_journal_doi":"10.3389/fevo.2026.1830790",
         "source_original_XLSX_SHA256":ORIGINAL_XLSX,
         "original_studies":1,"source_species":3,
         "primary_endpoint":"vial-level conditioned dead Alnus glutinosa leaf consumption",
         "not_plant_herbivory_on_living_wetland_vegetation":True,
         "not_an_ALAN_pollen_or_seed_outcome":True,
         "condition_comparison_vial_level_not_individual_animals":True,
         "light_treatments_lux":[0,1,10],
         "original_source_species_feeding":summary,
         "cross_species_descriptive_logratio_contrasts":cross,
         "bootstrap_CI_is_conditional_on_original_vials_not_independent_lamp_installations":True,
         "outcome_processing":"Skip source blank outcomes, not zero. Feeding values including exact zeros are retained; log ratio of strictly positive group arithmetic means.",
         "species_duration_and_normalization_differences_confounded":True,
         "not_an_ALAN_wet_species_main_mechanism":True}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    args=p.parse_args()
    x=analyze(fetch_original())
    file=Path(args.out);file.parent.mkdir(exist_ok=True,parents=True)
    file.write_text(json.dumps(x,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("CAMPOS2026_THREE_SPECIES_VIAL_FEEDING_SOURCE_REANALYSIS "+
          json.dumps({"source":x["source_original_XLSX_SHA256"],
                      "species":x["original_source_species_feeding"],
                      "contrast":x["cross_species_descriptive_logratio_contrasts"]},
                     sort_keys=True),flush=True)
if __name__=="__main__":main()
