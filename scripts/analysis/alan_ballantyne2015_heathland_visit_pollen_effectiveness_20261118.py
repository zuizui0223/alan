#!/usr/bin/env python3
"""True 2015 heathland plant×visitor effective-pollen network, NOT ALAN trial.

Source Ballantyne, Baldock & Willmer 2015 Proc R Soc B 282:20151130.
Original Dryad doi:10.5061/dryad.17pp3, version 9417,
visitation+single-visit conspecific pollen deposition in virgin stigmas.
SHA256 of original single CSV pinned to exact author source.

One Dorset heathland site observed 2013-14, five plant species and 1338
single visits. Source groups 1263 bees, 59 flies, 16 other/ants/butterfly,
NO moth visits. Reconstruct two distinct networks on the SAME 1338 events:
  V_ig = n visits / all observed visits to plant i
  PE_ig = sum control-corrected conspecific pollen / total source
          corrected pollen deposited on plant i (for that plant only).
Their ratio and effective guild diversity quantify functional vs numeric
role without treating pollination visits as independent artificial lights.

Methodological cautions:
* A non-visited plant or insect absent from THIS observation is no proof
  ecological absence.
* No ALAN exposure or light manipulation, no seed set or reproductive output.
* Flower-level pollen counts share dates, pollinator identities and plant
  individuals; date-cluster resampling is only descriptive.
* These are pollen grains arriving on virgin stigmas, not necessarily viable
  fertilizations; normalized plant-specific fractions cannot compare pollen
  production across species.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import io
import json
import math
import random
import statistics
import urllib.request
import zipfile
from pathlib import Path

URL="https://datadryad.org/api/v2/versions/9417/download"
SHA256="6994304ac377540c87b57560ecec50f72e8514494f9edd034eae2b5ed18f61c1"
SCHEMA="alan_ballantyne2015_original_heathland_visit_vs_pollen_effectiveness_v1"
SEED=20261118
B=3999
BEES={
    "Apis_mellifera","Bombus_hortorum","Bombus_jonellus",
    "Bombus_lapidarius","Bombus_pascuorum",
    "Bombus_terrestris/lucorum","halictid","other_solitary_bees"
}
FLIES={"Episyrphus","Eupeodes","large_syrphid","small_syrphid","muscid"}
OTHERS={"ant","butterfly","other"}
PLANTS={
  "Erica_tetralix":"wet heath ecological affinity, not local water table",
  "Erica_cinerea":"relatively drier heath species",
  "Calluna_vulgaris":"heathland shrub species",
  "Ulex_minor":"heathland gorse shrub",
  "Ulex_europaeus":"heathland gorse shrub"
}
EXPECTED={
    "Erica_tetralix":398,
    "Erica_cinerea":351,
    "Calluna_vulgaris":379,
    "Ulex_minor":126,
    "Ulex_europaeus":84
}

def require(cond,why):
    if not bool(cond):raise RuntimeError("ALAN_HEATH_ORIGINAL_SVD_STOP: "+why)

def groups(source):
    if source in BEES:return "bees"
    if source in FLIES:return "flies"
    if source in OTHERS:return "other"
    raise ValueError("Undocumented original visitor group: "+source)

def load():
    req=urllib.request.Request(URL,headers={
        "User-Agent":"Open ecological data check original author stigma pollen",
        "Accept":"application/zip,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=65) as response:
        raw=response.read(8<<20)
    require(zipfile.is_zipfile(io.BytesIO(raw)),
            "original public Dryad not a ZIP; do not accept HTML")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        require(z.testzip() is None,"Original Dryad ZIP CRC mismatch")
        names=[x for x in z.namelist() if x.lower().endswith(".csv")]
        require(len(names)==1,"Dryad pollen CSV not unique")
        buf=z.read(names[0])
    require(hashlib.sha256(buf).hexdigest()==SHA256,
            "Original 2015 Dryad plant visit-pollen bytes changed")
    import csv
    sample=list(csv.DictReader(io.StringIO(buf.decode("latin-1")),delimiter=","))
    require(len(sample)==1338,"Published 1338 raw stigmas not recovered")
    return sample

def normalize(sample):
    records=[]
    for r in sample:
        plant=r["plant_species"].strip()
        taxon=r["visitor_group"].strip()
        require(plant in PLANTS,"Unexpected source heathland plant taxon")
        guild=groups(taxon)
        dt=r["date"].strip()
        v=float(r["conspecific_pollen_deposited"])
        orig=float(r["conspecific_pollen_deposited_(before_control_values_removed)"])
        require(v>=0 and orig>=v and math.isfinite(v) and math.isfinite(orig),
                "Control-subtracted source pollen cannot be nonfinite/negative")
        require(dt!="","Original survey day missing")
        records.append({"plant":plant,"taxon":taxon,"guild":guild,
                        "date":dt,"conspecific_pollen":v,
                        "uncorrected_grains":orig})
    require(Counter(r["plant"] for r in records)==EXPECTED,
            "Source original plant counts changed")
    require(Counter(r["guild"] for r in records)=={
        "bees":1263,"flies":59,"other":16},"Original visitor guild counts changed")
    require(sum(r["conspecific_pollen"]==0 for r in records)==281,
            "Original author 281 control-adjusted non-depositing visits changed")
    return records

def safe_div(num,den):
    return num/den if den>0 else None

def effective_number(weights):
    v=[x for x in weights if x>0]
    total=sum(v)
    return 1/sum((x/total)**2 for x in v) if total>0 else None

def plant_summary(rows,plant):
    v=[r for r in rows if r["plant"]==plant]
    total=len(v)
    total_pollen=sum(r["conspecific_pollen"] for r in v)
    by_group=defaultdict(list)
    by_taxon=defaultdict(list)
    for row in v:
        by_group[row["guild"]].append(row)
        by_taxon[row["taxon"]].append(row)
    result={
       "plant":plant,
       "habitat_affinity_descriptor":PLANTS[plant],
       "original_flower_single_visit_events":total,
       "source_days":len({r["date"] for r in v}),
       "original_positive_single_visit_pollen_count":
          sum(r["conspecific_pollen"]>0 for r in v),
       "original_zero_after_controls_visits":
          sum(r["conspecific_pollen"]==0 for r in v),
       "original_total_conspecific_pollen_grains_after_controls":total_pollen,
       "source_pollen_weighted_vs_visit_weighted_visitor_concentration":{
           "effective_visitor_taxa_number_visitation_Simpson":
               effective_number([len(q) for q in by_taxon.values()]),
           "effective_visitor_taxa_number_stigmatic_pollen_Simpson":
               effective_number([sum(r["conspecific_pollen"] for r in q)
                                 for q in by_taxon.values()]),
           "pollen_weighted_vs_visit_diversity_is_one_observational_community":True
       },
       "groups":{}}
    for guild in ("bees","flies","other"):
        q=by_group[guild]
        n=len(q)
        pollen=sum(x["conspecific_pollen"] for x in q)
        if not n:
            result["groups"][guild]={
                "source_visits":0,"source_stigmatic_pollen_grains":0,
                "source_observed_guild_visitation_not_detected":True,
                "source_absence_in_nature_NOT_inferred":True}
            continue
        result["groups"][guild]={
          "source_visits":n,
          "source_days_with_visits":len({x["date"] for x in q}),
          "source_taxa_with_visits":len({x["taxon"] for x in q}),
          "source_single_visit_positive_deposition":sum(x["conspecific_pollen"]>0 for x in q),
          "fraction_non_depositing_visits":sum(x["conspecific_pollen"]==0 for x in q)/n,
          "source_mean_grains_per_visitor_virgin_stigma":pollen/n,
          "source_total_conspecific_pollen_grains":pollen,
          "visit_share":n/total,
          "pollen_deposition_share":safe_div(pollen,total_pollen),
          "functional_share_divided_by_visit_share":
              safe_div(safe_div(pollen,total_pollen),n/total)
        }
    bee=result["groups"]["bees"]
    fly=result["groups"]["flies"]
    result["fly_vs_bee_same_plant"]={
        "observed_n_fly_single_visits":fly["source_visits"],
        "observed_n_bee_single_visits":bee["source_visits"],
        "mean_fly_per_visit_pollen_relative_to_bee":
            safe_div(fly.get("source_mean_grains_per_visitor_virgin_stigma",0),
                     bee.get("source_mean_grains_per_visitor_virgin_stigma",0)),
        "source_direct_moth_SVD_measurement_absent":True,
        "source_guild_efficiency_mean_not_a_seed_set_or_ALAN_effect":True
    }
    # Per source recorded visitor taxon, not a global guild taxonomy.
    result["original_source_visitor_taxa"]={
        t:{"source_visits":len(q),
           "mean_conspecific_pollen_after_controls":
               sum(z["conspecific_pollen"] for z in q)/len(q),
           "source_days":len({z["date"] for z in q}),
           "source_non_depositing_visits":sum(z["conspecific_pollen"]==0 for z in q)}
        for t,q in sorted(by_taxon.items())}
    return result

def date_cluster_bootstrap(rows,plant):
    v=[r for r in rows if r["plant"]==plant]
    by_date=defaultdict(list)
    for r in v:by_date[r["date"]].append(r)
    dates=list(by_date)
    fdates={r["date"] for r in v if r["guild"]=="flies"}
    bdates={r["date"] for r in v if r["guild"]=="bees"}
    q={"source_independent_calendar_days":len(dates),
       "source_days_with_any_fly":len(fdates),
       "source_days_with_any_bee":len(bdates),
       "original_date_cluster_bootstrap_attempts":B,
       "interpretation":"descriptive calendar day cluster resampling only, not independent light intervention"}
    if len(dates)<5 or len(fdates)<4 or len(bdates)<4:
        q["status"]="source_fly_or_bee_day_replication_insufficient_STOP"
        return q
    rng=random.Random(SEED+sum(ord(x) for x in plant))
    vals=[]
    for i in range(B):
        chosen=rng.choices(dates,k=len(dates))
        bee=[];fly=[]
        for day in chosen:
            for record in by_date[day]:
                if record["guild"]=="bees":bee.append(record["conspecific_pollen"])
                elif record["guild"]=="flies":fly.append(record["conspecific_pollen"])
        if len(fly)>=2 and len(bee)>=2 and sum(bee)>0:
            vals.append(statistics.fmean(fly)/statistics.fmean(bee))
    q["valid_resamples_with_both_group_events"]=len(vals)
    if len(vals)<.95*B:
        q["status"]="bootstrap_sampled_fly_visits_too_rare_STOP"
        return q
    vals.sort()
    q["status"]="source_date_clustered_observational_95CI"
    q["source_date_bootstrap_fly_to_bee_mean_stigmatic_pollen_ratio_95CI"]=[
        vals[int(.025*(len(vals)-1))],vals[int(.975*(len(vals)-1))]]
    return q

def execute():
    r=normalize(load())
    plants={plant:plant_summary(r,plant) for plant in PLANTS}
    for plant in plants:
        plants[plant]["fly_vs_bee_same_plant"]["date_cluster_source_sensitivity"]=date_cluster_bootstrap(r,plant)
    return {
        "schema":SCHEMA,"source_Dryad_DOI":"10.5061/dryad.17pp3",
        "source_article_DOI":"10.1098/rspb.2015.1130",
        "source_exact_original_csv_sha256":SHA256,
        "source_original_plant_x_virgin_stigma_visits":len(r),
        "source_unique_plant_taxa":len(PLANTS),
        "source_original_field_community_count":1,
        "source_has_manipulated_ALAN_treatment":False,
        "source_original_taxonomic_groups":dict(Counter(z["guild"] for z in r)),
        "source_original_visitor_taxa":sorted({z["taxon"] for z in r}),
        "source_original_moth_visits":0,
        "source_original_non_depositing_visit_events":sum(z["conspecific_pollen"]==0 for z in r),
        "original_source_plant_results":plants,
        "no_plant_specific_moth_fly_effectiveness_ratio_estimated":True,
        "no_source_visits_reweighted_as_independent_ecosystem_replications":True,
        "mechanism_scope":"Visit frequency versus conspecific pollen deposited in one heath community, not an ALAN effect on wet-heath pollen or seed output",
        "author_source_study_method":"Virgin-flower stigmas; visitor groups; control-subtracted conspecific pollen counts; per-original-date cluster sensitivity where supported",
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    args=p.parse_args()
    out=execute()
    file=Path(args.out);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("HEATHLAND_ORIGINAL_PLANT_VISITS_VS_EFFECTIVE_POLLEN",json.dumps({
      "source_original_events":out["source_original_plant_x_virgin_stigma_visits"],
      "taxonomic_groups":out["source_original_taxonomic_groups"],
      "source_moths_n":out["source_original_moth_visits"],
      "plants":{p:{
          "days":a["source_days"],
          "visits":a["original_flower_single_visit_events"],
          "nondeposition":a["original_zero_after_controls_visits"],
          "bees":a["groups"]["bees"],
          "flies":a["groups"]["flies"],
          "fly_bee":a["fly_vs_bee_same_plant"],
          "Simpson":a["source_pollen_weighted_vs_visit_weighted_visitor_concentration"]
          } for p,a in out["original_source_plant_results"].items()}
    },sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
