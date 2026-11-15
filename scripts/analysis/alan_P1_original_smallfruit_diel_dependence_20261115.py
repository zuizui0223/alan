#!/usr/bin/env python3
"""Recover ORIGINAL diel-exclusion fruit outcomes, not ALAN treatment effects.

Blareau & Requier, J Applied Entomology, online 2025/issue 2026:
10.1111/jen.70019 ; source Figshare 10.6084/m9.figshare.29930195.v1.
Author files are downloaded and SHA256 pinned before aggregation.

Only group counts and non-identifying estimates are saved. No source rows or
plant identifiers are committed. All SE/CI for site/plant-level effects remain
unknown because files omit original individual plant ID.
"""
from __future__ import annotations

import collections
import csv
import hashlib
import io
import json
import math
import pathlib
import re
import urllib.request

FILES = {
    "fruit_set": {
        "url": "https://ndownloader.figshare.com/files/57253715",
        "sha256": "c5c4a169b345cccafa6d555cfa1726ae5f6606bb83312c2602bee3626cc2dd0f",
        "columns": ["Plant_sp","Treatment","Fruit_set_success","Fruit_set_failure","Year"],
        "expected_rows": 80,
    },
    "fruit_quality": {
        "url": "https://ndownloader.figshare.com/files/57253712",
        "sha256": "b988588bf14b73ca2f619845e09eedc92ef464b77e64417431fbb2b725334fd5",
        "columns": ["Plant_sp","Treatment","Mass.g","Height.mm","Width.mm","Fruit_vol.mm3","Seed_set","Year"],
        "expected_rows": 223,
    },
}

def http_bytes(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent":"ALAN-P1-Figshare-original-replay/2026-11-15 (https://github.com/zuizui0223/alan)"})
    with urllib.request.urlopen(req,timeout=35) as response:
        data=response.read(2_000_000)
    if not data:
        raise RuntimeError("Empty source payload: "+url)
    return data


def read_original(name: str, getter=http_bytes) -> list[dict]:
    config=FILES[name]
    data=getter(config["url"])
    digest=hashlib.sha256(data).hexdigest()
    if digest!=config["sha256"]:
        raise ValueError("Original source SHA mismatch; cannot run: "+name)
    reader=csv.DictReader(io.StringIO(data.decode("utf-8-sig")),delimiter=";")
    if reader.fieldnames!=config["columns"]:
        raise ValueError("Source schema unexpected: "+name+" "+repr(reader.fieldnames))
    rows=list(reader)
    if len(rows)!=config["expected_rows"]:
        raise ValueError("Source row count mismatch: "+name+" "+str(len(rows)))
    return rows


def int_count(x: str) -> int:
    value=str(x).strip()
    if not re.fullmatch(r"\d+(?:\.0+)?",value):
        raise ValueError("Noninteger success/failure without original audit: "+repr(value))
    return int(float(value))


def float_or_none(x: str) -> float | None:
    value=str(x).strip()
    if value.lower() in ("","na","nan","null","-"):
        return None
    try:
        ans=float(value)
    except ValueError:
        return None
    return ans if math.isfinite(ans) else None


def normalize_treatment(value: str) -> str:
    norm=re.sub(r"[^a-z0-9+]","",str(value).strip().lower().replace("&","+"))
    if norm in ("closed","bagged","no","none","excluded","exclusion","fullclosed"): return "closed"
    if norm in ("day","dayonly","diurnal"): return "day_only"
    if norm in ("night","nightonly","nocturnal"): return "night_only"
    if norm in ("day+night","daynight","open","both","dayandnight","full","day+nightonly","control"): return "day_night_open"
    return "UNRESOLVED("+str(value)+")"


def reduce_source(fruit_rows: list[dict], quality_rows: list[dict]) -> dict:
    groups=collections.defaultdict(lambda:{"inflorescences":0,"fruit_success":0,"fruit_failure":0})
    all_treatments=set()
    for r in fruit_rows:
        t=normalize_treatment(r["Treatment"])
        year=str(r["Year"]).strip()
        taxon=str(r["Plant_sp"]).strip()
        if not (year and taxon):
            raise ValueError("Missing taxon or year; no silent rows")
        key=(taxon,year,t)
        s=int_count(r["Fruit_set_success"])
        f=int_count(r["Fruit_set_failure"])
        all_treatments.add(t)
        if s+f <= 0:
            raise ValueError("Empty original inflorescence row")
        groups[key]["inflorescences"]+=1
        groups[key]["fruit_success"]+=s
        groups[key]["fruit_failure"]+=f
    quality=collections.defaultdict(lambda:{"measured_fruit_rows":0,"nonmissing_seed_rows":0,"seed_count_sum_measured":0,"unknown_seed_rows":0})
    for r in quality_rows:
        key=(str(r["Plant_sp"]).strip(),str(r["Year"]).strip(),normalize_treatment(r["Treatment"]))
        quality[key]["measured_fruit_rows"]+=1
        value=float_or_none(r["Seed_set"])
        if value is None:
            quality[key]["unknown_seed_rows"]+=1
        elif value<0:
            raise ValueError("Negative seed count")
        else:
            quality[key]["nonmissing_seed_rows"]+=1
            quality[key]["seed_count_sum_measured"]+=value
    group_rows=[]
    for key in sorted(set(groups)|set(quality)):
        raw=groups.get(key,{"inflorescences":0,"fruit_success":0,"fruit_failure":0})
        q=quality.get(key,{"measured_fruit_rows":0,"nonmissing_seed_rows":0,"seed_count_sum_measured":0,"unknown_seed_rows":0})
        initial=raw["fruit_success"]+raw["fruit_failure"]
        row={
          "taxon_original":key[0],"year":key[1],"treatment":key[2],
          **raw, "initial_flowers":initial,
          "fruit_set_probability_raw":None if initial==0 else raw["fruit_success"]/initial,
          **q, "quality_rows_to_fruit_success_ratio_NOT_a_coverage_fraction":None if raw["fruit_success"]==0 else q["measured_fruit_rows"]/raw["fruit_success"],
          "quality_vs_fruit_count_status":("MATCH" if q["measured_fruit_rows"]==raw["fruit_success"] else ("QUALITY_ROWS_EXCEED_FRUIT_SET_SUCCESS" if q["measured_fruit_rows"]>raw["fruit_success"] else "QUALITY_ROWS_FEWER_THAN_FRUIT_SET_SUCCESS")),
          "safe_to_join_quality_to_initial_flowers":False,
          "seed_count_per_measured_fruit":None if not q["nonmissing_seed_rows"] else q["seed_count_sum_measured"]/q["nonmissing_seed_rows"]
        }
        row["se_and_ci_cluster_adjusted"]=None
        group_rows.append(row)
    index={(x["taxon_original"],x["year"],x["treatment"]):x for x in group_rows}
    estimands=[]
    pairs=sorted(set((x["taxon_original"],x["year"]) for x in group_rows))
    for taxon,year in pairs:
        day=index.get((taxon,year,"day_only"))
        open_=index.get((taxon,year,"day_night_open"))
        pday=None if not day else day["fruit_set_probability_raw"]
        popen=None if not open_ else open_["fruit_set_probability_raw"]
        V=None if not popen or pday is None else 1-pday/popen
        estimands.append({
            "taxon_original":taxon,"year":year,
            "p_day_only":pday, "p_day_night_open":popen,
            "V_fruit_night_access":V,
            "source_comparison":"night-access marginal functional dependence; NO ALAN manipulation",
            "SE_or_CI":None,
            "nocturnal_insect_ever_observed":False,
            "cannot_assign_moth_effectiveness":True
        })
    species_code_mapping={"BC":"Ribes nigrum","RB":"Rubus idaeus","RD":"Ribes rubrum"}
    expected_flower_counts={"BC":375,"RB":54,"RD":421}
    observed_flower_counts={code:sum(r["initial_flowers"] for r in group_rows if r["taxon_original"]==code) for code in species_code_mapping}
    if len(fruit_rows)==FILES["fruit_set"]["expected_rows"] and observed_flower_counts != expected_flower_counts:
        raise AssertionError("Original paper flower totals disagree: "+repr(observed_flower_counts))
    flagged=[{"taxon_original":x["taxon_original"],"year":x["year"],"treatment":x["treatment"],"fruit_set_success":x["fruit_success"],"quality_rows":x["measured_fruit_rows"],"status":x["quality_vs_fruit_count_status"]} for x in group_rows if x["quality_vs_fruit_count_status"]!="MATCH"]
    return {
      "schema":"ALAN_P1_B1_2025_original_fruit_crop_diel_gate_v2",
      "species_original_code_mapping_confirmed_by_published_flower_totals":species_code_mapping,
      "source_paper_species_flower_totals_verified":observed_flower_counts,
      "source_quality_join_mismatch_group_count":len(flagged),
      "source_quality_join_mismatch_groups":flagged,
      "source_article":"10.1111/jen.70019","source_dataset":"10.6084/m9.figshare.29930195.v1",
      "source_sha256":{k:x["sha256"] for k,x in FILES.items()},
      "original_rows":{k:x["expected_rows"] for k,x in FILES.items()},
      "n_inflorescences":sum(x["inflorescences"] for x in group_rows),
      "n_initial_flowers":sum(x["initial_flowers"] for x in group_rows),
      "n_original_fruits_success":sum(x["fruit_success"] for x in group_rows),
      "n_measured_fruit_quality_rows":sum(x["measured_fruit_rows"] for x in group_rows),
      "treatments_raw_normalized":sorted(all_treatments),
      "sample_counts_by_taxon_year_treatment":group_rows,
      "new_night_access_V_by_taxon_year":estimands,
      "source_contains_individual_plant_id":False,
      "conditional_seed_coverage_issue":"Quality file is conditioned on formed fruit but multiple taxon-year-treatment counts disagree with fruit-set success, some quality rows exceed fruit successes. Never call 223/258 coverage probability nor combine quality seed sums into total per originally assigned flower without source row links.",
      "ALAN_risk_ratio_estimated":False,
      "night_pollinator_guild_confirmed":False,
      "meta_reproductive_beta_V":None,
      "interpretation":"Cross-stage only, noninferential until plant-level correlation, missing fruit outcomes and site-year data reconciled. Negative V values are retained."
    }


def main():
    fr=read_original("fruit_set")
    q=read_original("fruit_quality")
    out=reduce_source(fr,q)
    path=pathlib.Path("results/alan_P1_B1_smallfruit_genuine_diel_V_20261115.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
      "source_rows":out["original_rows"],"inflorescences":out["n_inflorescences"],
      "quality_group_mismatches":out["source_quality_join_mismatch_group_count"],
      "initial_flowers":out["n_initial_flowers"],
      "formed_fruits":out["n_original_fruits_success"],
      "quality_rows":out["n_measured_fruit_quality_rows"],
      "treatments":out["treatments_raw_normalized"],
      "V":out["new_night_access_V_by_taxon_year"],
      "ALAN_effect":None
    },ensure_ascii=False))


if __name__=="__main__":
    main()
