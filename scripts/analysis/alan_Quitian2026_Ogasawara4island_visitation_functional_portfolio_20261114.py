#!/usr/bin/env python3
"""Real Ogasawara 2021 four-island same-season insect pollinator-channel audits.

True ORIGINAL Quitian et al. 2026 study DOI 10.1111/cobi.70304; source author
dataset Zenodo 10.5281/zenodo.19221853 MD5 pinned by source module.
Not ALAN-exposed, not a night-moth survey, not a direct pollen delivery/seed
measurement. Original surveys were DAYTIME 08:30–15:30 in May/Jul/Sep 2021.
Legitimate visit = floral reproductive parts touched, not necessarily effective.

Motivation: get source-measured island community channel *baseline*, not infer
predator compensation from global GloBI or light sensitivity from sampled
islands. Measure honeybee source visit share, non-honeybee order composition,
per-plant effective partner diversity and exact plant×season common support
across anole presence/absence and islands. Test island/anole confounding and
avoid pretending 2,745 source species-visit rows are independent islands.
"""
from __future__ import annotations
import argparse,collections,hashlib,io,json,math
from pathlib import Path
from openpyxl import load_workbook
from scripts.analysis.alan_Quitian2026_original_Ogasawara_4island_network_source_manifest_20261114 import (
 source, SOURCE_MD5, DOI, PAPER
)

COLUMNS=("Island","Invasional Context","Season","Forest_Status","Anole",
         "Plant_sp","Clase","Order","Ord_Fam","Poll_sp","N.Int")
ISLANDS=("A_Chichijima","B_Hahajima","C_Anijima","D_Ototojima")
SEASONS=("A_MAY","B_JULY","C_SEP")

def require(z,msg):
    if not bool(z):raise ValueError("ALAN_OGASAWARA2021_ORIGINAL: "+str(msg))
def entropy_effective(vis):
    n=sum(vis.values())
    if n<=0:return None
    a=[x/n for x in vis.values() if x>0]
    return float(math.exp(-sum(p*math.log(p) for p in a)))
def normalize(s):
    return " ".join(str(s or "").strip().split())
def api(s):
    t=s.strip().lower()
    return t=="apis mellifera" or t.startswith("apis mellifera ") or t=="a. mellifera"

def read():
    b=source()
    workbook=load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    require(workbook.sheetnames==["Data interactions"],"Original XLSX sheet name changed")
    sheet=workbook.active
    rows=sheet.iter_rows(values_only=True)
    columns=tuple(normalize(x) for x in next(rows))
    require(columns==COLUMNS,"Original source sheet exact 11-column schema changed")
    out=[]
    for vals in rows:
        if all(v is None for v in vals):continue
        x={k:v for k,v in zip(COLUMNS,vals)}
        for name in COLUMNS[:-1]:
            require(x[name] is not None,"Missing required original document field "+name)
            x[name]=normalize(x[name])
        try:
            count=float(x["N.Int"])
            require(math.isfinite(count) and count>=0 and int(count)==count,
                    "Original noninteger or negative interactions")
            x["N.Int"]=int(count)
        except (ValueError,TypeError) as e:
            raise ValueError("Original N.Int source field nonnumeric") from e
        out.append(x)
    workbook.close()
    require(len(out)==2745,"Original 2745 source plant-pollinator rows changed")
    require(set(x["Island"] for x in out)==set(ISLANDS),"Original source island list changed")
    require(set(x["Season"] for x in out)==set(SEASONS),"Original source survey seasons changed")
    require(set(x["Anole"] for x in out)=={"Presence","Absence"},
            "Original source predator field changed")
    return out,hashlib.sha256(b).hexdigest()

def group(rows,keys):
    out=collections.defaultdict(list)
    for r in rows:
        out[tuple(r[k] for k in keys)].append(r)
    return out

def model(rows):
    total=sum(r["N.Int"] for r in rows)
    plants={r["Plant_sp"] for r in rows}
    pollinators={r["Poll_sp"] for r in rows}
    actual={"source_interaction_rows":len(rows),
            "source_genuine_recorded_legitimate_visits":total,
            "source_positive_plant_names":len(plants),
            "source_positive_pollinator_names_and_morphotypes":len(pollinators),
            "source_apis_name_variants":sorted({r["Poll_sp"] for r in rows if api(r["Poll_sp"])}),
            "original_pollinator_class_levels":sorted({r["Clase"] for r in rows}),
            "original_order_levels":sorted({r["Order"] for r in rows}),
            "original_family_levels":len({r["Ord_Fam"] for r in rows}),
            "source_unique_contexts":len({(r["Island"],r["Invasional Context"]) for r in rows}),
            "source_distinct_anole_status":sorted({r["Anole"] for r in rows}),
            "survey_year":2021,
            "source_plant_pollinator_visits_NOT_fruit_seed_efficacy":True}
    require(total>=3000 and total<=5500,"Original total visits unexpectedly drifted")
    require(48<=len(plants)<=75 and 100<=len(pollinators)<=135,
            "Original plant or pollinator species source support changed")
    # Use original invaded site CONTEXT within each physical island.
    context=group(rows,["Island","Invasional Context"])
    ctable=[]
    for (island,ctx),rs in sorted(context.items()):
        status={r["Anole"] for r in rs}
        forest={r["Forest_Status"] for r in rs}
        require(len(status)==1 and len(forest)<=2,
                "Original anole status not unique to site context")
        visits=sum(r["N.Int"] for r in rs)
        honey=sum(r["N.Int"] for r in rs if api(r["Poll_sp"]))
        nonbee=sum(r["N.Int"] for r in rs if not api(r["Poll_sp"]))
        ctable.append({"anonymous_original_island":island,
                       "original_invasion_context":ctx,
                       "anole_status":next(iter(status)),
                       "original_forest_status_levels":sorted(forest),
                       "legitimate_visits":visits,
                       "honeybee_Apis_mellifera_visits":honey,
                       "honeybee_fraction_of_source_visits":honey/visits if visits else None,
                       "other_positive_visits_not_assumed_native":nonbee,
                       "distinct_plant_names":len({r["Plant_sp"] for r in rs}),
                       "distinct_pollinator_source_names":len({r["Poll_sp"] for r in rs}),
                       "original_surveys_are_DAYTIME_only":True})
    seas=group(rows,["Island","Season"])
    island_season=[]
    for (island,season),rs in sorted(seas.items()):
        n=sum(x["N.Int"] for x in rs)
        honey=sum(x["N.Int"] for x in rs if api(x["Poll_sp"]))
        island_season.append({"original_island":island,"original_season":season,
            "original_interaction_rows":len(rs),
            "total_legitimate_visits":n,
            "Apis_mellifera_visits":honey,
            "honeybee_source_visit_fraction":honey/n if n else None,
            "distinct_original_plants":len({r["Plant_sp"] for r in rs}),
            "distinct_original_visitors":len({r["Poll_sp"] for r in rs}),
            "source_order_visits":dict(sorted(collections.Counter({
              x["Order"]:sum(y["N.Int"] for y in rs if y["Order"]==x["Order"])
              for x in rs}).items())),
            "not_direct_fitness_or_effective_pollen":True})
    plant_season=group(rows,["Island","Invasional Context","Season","Plant_sp"])
    portfolio=[]
    for k,rs in sorted(plant_season.items()):
        visits=collections.Counter()
        for x in rs:visits[x["Poll_sp"]]+=x["N.Int"]
        n=sum(visits.values())
        if n<=0:continue
        honey=sum(x["N.Int"] for x in rs if api(x["Poll_sp"]))
        portfolio.append({"island":k[0],"context":k[1],"season":k[2],
            "plant":k[3],"recorded_visits":n,
            "Apis_visit_share":honey/n,
            "pollinator_partner_number_with_positive_counts":len(visits),
            "effective_pollinator_partner_number_expShannon":entropy_effective(visits),
            "original_floral_abundance_NOT_IN_THIS_XLSX":True})
    # Never export plant-specific identity, site/year rows in artifact.
    by_ctx=group(portfolio,["island","context"])
    for item in ctable:
        q=by_ctx[(item["anonymous_original_island"],
                  item["original_invasion_context"])]
        item["source_plant_season_opportunities_count"]=len(q)
        item["mean_effective_pollinator_partners_per_documented_plant_season"]=(
            sum(x["effective_pollinator_partner_number_expShannon"] for x in q)/len(q) if q else None)
        item["median_documented_plant_season_Apis_visit_share"]=(
            sorted(x["Apis_visit_share"] for x in q)[len(q)//2] if q else None)
    anole_by_island={}
    for island in ISLANDS:
        cats=sorted(set(x["anole_status"] for x in ctable if x["anonymous_original_island"]==island))
        anole_by_island[island]={"observed_anole_states":cats,
                  "both_presence_absence_observed_WITHIN_island":len(cats)==2}
    # Independent ecological island count is FOUR, even if contexts/seasons/rows are larger.
    n_cross=sum(x["both_presence_absence_observed_WITHIN_island"] for x in anole_by_island.values())
    plant_matched=[]
    # Same named plant and same calendar season observed on both side(s) of Anole state;
    # these records are NOT source-independent experimental predator assignments.
    pool=collections.defaultdict(lambda:collections.defaultdict(set))
    for k,rs in plant_season.items():
        isld,ctx,sea,plant=k
        state=next(r["Anole"] for r in rs)
        pool[(sea,plant)][state].add(isld)
    for k,d in pool.items():
        if set(d)=={"Presence","Absence"}:
            plant_matched.append(k)
    conclusion={
      "schema":"ALAN_Ogasawara2021_four_islands_real_legitimate_plant_pollinator_network_v1",
      "status":"ORIGINAL_ZENODO_REAL_2021_ISLAND_CHANNEL_PORTFOLIO_WITHOUT_ALAN",
      "source":{
        "paper_DOI":PAPER,"original_dataset_DOI":DOI,
        "original_Zenodo_XLSX_MD5":SOURCE_MD5,
        "original_XLSX_SHA256":None,
        "survey_months":["May","July","September"],
        "original_year":2021,
        "survey_hours_local":"08:30–15:30",
        "one_archipelago_four_islands_NOT_4_globally_independent_regions":True,
        "original_legitimate_visit_not_actual_heterospecific_or_conspecific_pollen":True
      },
      "original_source_support":actual,
      "source_island_contexts":ctable,
      "source_island_seasons":island_season,
      "island_anole_field_design":anole_by_island,
      "islands_with_actual_within_island_anole_state_variation":n_cross,
      "plant_x_calendar_season_with_positive_both_predator_state_source_support":len(plant_matched),
      "independent_islands":4,
      "total_original_contexts":len(ctable),
      "plant_season_original_documented_records_no_raw_rows_exported":len(portfolio),
      "interpretation":{
        "native_invasive_status_of_nonApis_not_inferred_without_taxonomy":True,
        "original_source_visited_plant_not_verified_effective_pollen_or_seed":True,
        "anole_presence_not_direct_moth_or_honeybee_predation_observation":True,
        "anole_presence_and_invasion_status_may_be_confounded_by_island_habitat":True,
        "daytime_only_original_census_makes_night_moth_shift_UNOBSERVED":True,
        "no_source_field_nightlight_exposure":True,
        "no_moth_or_other_pollinator_count_missing_means_absent":True,
        "no_light_by_predator_seed_mediation_test":True,
        "effective_partner_number_is_VISIT_NETWORK_niche_not_pollen_function":True,
        "no_2019_2025_photo_ALAN_same_island_year_join_from_2021_field_study":True,
        "source_only_independent_archipelago_and_sites_not_global_generalization":True
      }
    }
    return conclusion

def main(out):
    out=Path(out)
    require(not out.exists(),"Refuse overwriting author-source results")
    data,sha=read()
    result=model(data)
    result["source"]["original_XLSX_SHA256"]=sha
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ALAN_ISLAND_GENUINE_OGASAWARA_POLLINATOR_CHANNELS",
          json.dumps(result,ensure_ascii=False),flush=True)
    return result
if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    main(a.out)
