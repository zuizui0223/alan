#!/usr/bin/env python3
"""Giavi2021 ALAN daytime visitors (6 paired Swiss meadows) x Kew F trait.

Original fixed 21 original plant-taxon source effects in alan/main, source
Git blob SHA1 b54cadf945edf2eac129d1b45ce9d5edcbda78cf.
Original RBG Kew BIFloraExplorer 3227 source taxa F_PLANTATT and F_Doring
source Git blob SHA1 e6a1da929af81f370d2a1b79177b489e3e573a06.
No unresolved genus (sp), aggregate (aggr) or sensu lato labels accepted.
The 21 source plant rows are not 21 independent light treatments.
Compute exact plant moisture-label rank permutational association in all
scored named species and 2+ shared dark-light meadow pairs sensitivity.
These are not direct reproductive seed outcomes or mediation of macro wet bias.
"""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,math
from pathlib import Path
from scipy.stats import spearmanr
from scripts.analysis.alan_anic2022_13species_kew_plantatt_20261120 import (
 RAW_GIT_BLOB,pull_kew_binary,trait_lookup
)
SCHEMA="alan_Giavi2021_6meadow_planttaxa_Kew_PLANTATT_F_original_source_support_v1"
GIAVI_FILE="results/all_species_forest/wet_environment/giavi_diurnal.csv"
GIAVI_BLOB_SHA1="b54cadf945edf2eac129d1b45ce9d5edcbda78cf"
def require(cond,msg):
    if not cond: raise RuntimeError("GIAVI_Kew_WETNESS_SOURCE_STOP: "+msg)
def source():
    b=Path(GIAVI_FILE).read_bytes()
    sha=hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()
    require(sha==GIAVI_BLOB_SHA1,"Giavi source blob SHA changed")
    a=list(csv.DictReader(b.decode("utf-8-sig").splitlines()))
    require(len(a)==21,"21 source taxa changed")
    require(len({r["label"] for r in a})==21,"duplicate plant taxon")
    return a
def source_name(label):
    if label.endswith("_sp") or label.endswith("_sp.") or "_aggr" in label or "_s.l." in label:
        return None
    a=label.split("_")
    return " ".join(a) if len(a)==2 and all(a) else None
def add_doring_only_records(raw,trait_map):
    import io
    source=csv.DictReader(io.StringIO(raw.decode("latin1")))
    for row in source:
        original_name=row["taxon_name_binom"].strip()
        f=row["F_Doring"].strip()
        if f in ("","NA"):continue
        try:value=float(f)
        except ValueError:continue
        if not (math.isfinite(value) and 1<=value<=12):continue
        if original_name not in trait_map:
            # Original PLANTATT F is unavailable; keep it NULL.
            trait_map[original_name]={
                "F_PLANTATT":None,"F_Doring":f,
                "Kew_taxonomic_status":row["taxonomic_status"].strip(),
                "Kew_aggregate_label":row["aggregate"].strip()}
    return trait_map

def matches(rows,traits):
    out=[]
    for r in rows:
        name=source_name(r["label"])
        p=traits.get(name) if name else None
        def val(x):
            if x is None or x in ("","NA"):return None
            return float(x)
        out.append({
            "original_taxon":r["label"],
            "resolved_original_binomial":name,
            "shared_lit_dark_meadow_pairs":int(r["shared_meadow_pairs"]),
            "source_lit_minus_dark_day_visitation_log1p":float(r["estimate"]),
            "F_PLANTATT":val(p["F_PLANTATT"]) if p else None,
            "F_Doring":val(p["F_Doring"]) if p else None,
            "unresolved_original_plant_is_not_called_dry":name is None,
            "not_an_independent_light_experiment":True,
        })
    return out
def rho(rows,key):
    a=[r[key] for r in rows]
    b=[r["source_lit_minus_dark_day_visitation_log1p"] for r in rows]
    if len(a)<3 or len(set(a))<=1:return None
    return float(spearmanr(a,b).statistic)
def exact(rows,key):
    if len(rows)<3:return {"status":"NO_REPLICATED_SOURCE_SUPPORT","n":len(rows)}
    f=[r[key] for r in rows]
    y=[r["source_lit_minus_dark_day_visitation_log1p"] for r in rows]
    obs=rho(rows,key)
    require(obs is not None,"source trait lacks variance")
    # Capped to small 8-9-species exact distinct assignments;
    # tests a mathematical plant-label null, NOT original 6 random lamps.
    pm=set(itertools.permutations(f))
    require(0<len(pm)<=200000,"source null permutations unsupported")
    tail=sum(abs(spearmanr(q,y).statistic)>=abs(obs)-1e-12 for q in pm)
    return {
       "status":"EXACT_ORDINAL_SPECIES_MOISTURE_LABEL_NULL_EXPLORATORY",
       "species":len(rows),
       "physical_meadow_light_pairs_total":6,
       "observed_Spearman_rho":obs,
       "distinct_source_label_arrangements":len(pm),
       "two_sided_exact_label_permutation_p":tail/len(pm),
       "wet_F_gt7":sum(z>7 for z in f),
       "drier_F_le7":sum(z<=7 for z in f),
       "no_causal_independent_species_light_effect":True
    }
def analyze(records,traits):
    allrows=matches(records,traits)
    output={}
    for key in ("F_PLANTATT","F_Doring"):
        f=[r for r in allrows if r[key] is not None]
        strict=[r for r in f if r["shared_lit_dark_meadow_pairs"]>=2]
        output[key]={
            "source_names_with_F":len(f),
            "source_names_with_2plus_meadow_pairs":len(strict),
            "all_names_matched_source_exact":exact(f,key),
            "strict_2plus_meadow_pairs":exact(strict,key),
            "source_scored_original_species_and_pair_count":[
                {"source_name":z["resolved_original_binomial"],"F":z[key],
                 "original_shared_meadow_pairs":z["shared_lit_dark_meadow_pairs"]}
                for z in f],
            "leave_one_species_out_rho":{
                z["resolved_original_binomial"]:rho([v for v in f if v is not z],key)
                for z in f}
        }
    return {
        "schema":SCHEMA,
        "Giavi_original_data_doi":"10.5281/zenodo.4540407",
        "Giavi_paper_doi":"10.1038/s41467-021-22011-8",
        "Giavi_source_21_taxon_git_blob_sha1":GIAVI_BLOB_SHA1,
        "Kew_BIFloraExplorer_original_git_blob_sha1":RAW_GIT_BLOB,
        "source_original_six_physical_meadow_pairs":6,
        "source_original_plant_taxon_labels":len(allrows),
        "source_true_named_species":sum(x["resolved_original_binomial"] is not None for x in allrows),
        "source_original_plant_label_and_f_score":allrows,
        "two_independent_UK_F_score_sensitivities":output,
        "source_original_meadow_by_plant_covariance_not_recovered":True,
        "source_visitation_not_compatible_pollen_or_seeds":True,
        "moist_species_bias_causal_mechanism_not_estimated":True,
        "cannot_promote_to_general_wetness_moderator_with_fewer_4wet4dry":True,
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    b=pull_kew_binary();z,n=trait_lookup(b);z=add_doring_only_records(b,z)
    result=analyze(source(),z)
    result["original_Kew_plant_rows"]=n
    path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GIAVI2021_SIXMEADOW_Kew_F_MOIST_DAYVISITATION_REAL_SOURCE_AUDIT "+
      json.dumps({"source_true_named_species":result["source_true_named_species"],
                  "source_F_sensitivities":result["two_independent_UK_F_score_sensitivities"]},
                 ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
