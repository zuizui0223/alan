#!/usr/bin/env python3
"""Macgregor et al. (2016) original EIDC: matched lit/unlit moth diagnostics.

Source DOI 10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708.
This audits the separate measures of over-lamp flight, captured moths
(transect/light trap), and explicitly recorded body-pollen. It DOES NOT
infer unsampled moth absences or effective pollen deposition on flowers.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import t as student_t, ttest_1samp, spearmanr

URL="https://data-package.ceh.ac.uk/data/31cc5cec-d33b-4dd6-a932-061ff947e708.zip"
SOURCE_FILES={
    "sites":("data/F1_sites.csv","0584befa8a5068be620ba977f179d55238222ddecab492076ff8864d2238c6e6",22),
    "moths":("data/F2_moths.csv","6a460e3314a51a609523553b04248416aa4494afe411f70d63ab702991cadc11",1599),
    "flight":("data/F3_overhead.csv","866e4257b8aed7652ee5b5738ed17e5023ad1aae8af2f4ef0d887e3bd4455b79",90),
    "pollen":("data/F4_pollen.csv","e205e9d7320ec12e97e4d8bcbd1e31ed03fbde92d03dcf4688f520e595b9735b",361)}
SCHEMA="ALAN_true_Macgregor2016_moth_flight_capture_pollen_matched_20site_v1"

def require(x,reason):
    if not x: raise ValueError("ALAN_EIDC_TRUE_MOTH: "+str(reason))

def signed_tables():
    r=requests.get(URL,timeout=(35,110),headers={"User-Agent":"ALAN-independently-source-audited-ecology/1.0"})
    r.raise_for_status()
    require(zipfile.is_zipfile(io.BytesIO(r.content)),"Original EIDC archive unavailable")
    frames={}
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        require(z.testzip() is None,"Original EIDC zip CRC failed")
        for key,(name,digest,n) in SOURCE_FILES.items():
            content=z.read(name)
            require(hashlib.sha256(content).hexdigest()==digest,"Original EIDC CSV changed: "+key)
            d=pd.read_csv(io.BytesIO(content))
            require(len(d)==n,"Original EIDC row count changed: "+key)
            frames[key]=d
    return frames

def paired_site_effect(pair):
    """Per-site mean log ratio, one independent matched site among replicates."""
    x=pair.copy()
    require(len(x)>=3 and {"Site","Lit","Unlit"}.issubset(x.columns),"No actual pair")
    x["log_ratio"]=np.log(x["Lit"]+.5)-np.log(x["Unlit"]+.5)
    per_site=x.groupby("Site",sort=True).log_ratio.mean()
    G=len(per_site)
    require(G>=3,"Too few true matched field sites")
    effect=float(per_site.mean())
    sd=float(per_site.std(ddof=1))
    se=sd/np.sqrt(G)
    critical=float(student_t.ppf(.975,G-1))
    p=float(ttest_1samp(per_site,0).pvalue) if sd>0 else None
    return {
        "independent_field_sites":G,"original_matched_sampling_events":len(x),
        "lit_more_events":int((x.Lit>x.Unlit).sum()),
        "unlit_more_events":int((x.Lit<x.Unlit).sum()),
        "exact_ties":int((x.Lit==x.Unlit).sum()),
        "per_site_mean_log_lit_vs_unlit_ratio_with_0p5":effect,
        "per_site_equal_exp_log_ratio":float(np.exp(effect)),
        "t_reference_95CI_exp_log_ratio":[float(np.exp(effect-critical*se)),
                                         float(np.exp(effect+critical*se))],
        "unadjusted_two_sided_t_p":p,
        "CI_unit":"site not each within-site observed-night repeat",
        "caution":"exploratory matched observational; log ratio has a fixed 0.5 pseudocount"
    }

def run(out):
    frames=signed_tables()
    sites,moths,flight,pollen=(frames[k] for k in ("sites","moths","flight","pollen"))
    require(moths.Reference.notna().all() and moths.Reference.nunique()==1599,
            "Individual recorded captured moth ID not unique")
    require(set(moths.Treatment.dropna().astype(str))=={"Lit","Unlit"},
            "New treatment categories")
    require(set(moths.Method.dropna().astype(str))=={"Transect","Trap"},
            "New sampling methods")
    for column in ["FirstCount","SecondCount","ThirdCount"]:
        flight[column]=pd.to_numeric(flight[column],errors="raise")
        require(flight[column].notna().all() and flight[column].ge(0).all(),
                "No invented overhead flight-count zeros")
    flight["observed_passes"]=flight[["FirstCount","SecondCount","ThirdCount"]].sum(axis=1)
    require(int(flight.observed_passes.sum())==434 and
            set(flight.Treatment)=={"Lit","Unlit"},"Unexpected source overhead counts")
    require(not flight.duplicated(["Site","Date","Treatment"]).any(),
            "Overhead flight sample duplicates a true matched session")
    paired=flight.pivot(index=["Site","Date"],columns="Treatment",
                        values="observed_passes").reset_index()
    require(len(paired)==45 and paired[["Lit","Unlit"]].notna().all().all(),
            "True original flight pairs were not all measured")
    require(int(paired.Lit.sum())==286 and int(paired.Unlit.sum())==148,
            "Source overhead counts changed")
    overhead={
        "n_original_overhead_rows":len(flight),
        "n_actual_matched_field_site_by_observation_date_pairs":len(paired),
        "number_source_field_sites":int(paired.Site.nunique()),
        "original_pooled_lit_passes":int(paired.Lit.sum()),
        "original_pooled_unlit_passes":int(paired.Unlit.sum()),
        "raw_ratio_lit_to_unlit":float(paired.Lit.sum()/paired.Unlit.sum()),
        "pair_inference":paired_site_effect(paired)
    }
    # Do not put an unobserved zero-night in the model: F2 is a moth INDIVIDUAL
    # file, not a full negative-collection-effort ledger.
    captures={}
    for method,d in moths.groupby("Method",sort=True):
        counted=d.groupby(["Site","Date","Treatment"]).size().rename("moths").reset_index()
        p=counted.pivot(index=["Site","Date"],columns="Treatment",values="moths").reset_index()
        both=p.loc[p[["Lit","Unlit"]].notna().all(axis=1)].copy()
        captures[method]={
            "n_individual_original_moth_rows":len(d),
            "n_distinct_source_moth_species":int(d.BinomialName.dropna().nunique()),
            "raw_source_capture_counts":{"Lit":int(d.Treatment.eq("Lit").sum()),
                                         "Unlit":int(d.Treatment.eq("Unlit").sum())},
            "n_site_date_events_with_nonzero_both_treatments":len(both),
            "matched_nonzero_pairs_are_not_complete_effort_zeros":True,
            "paired_positive_only_diagnostic":
                paired_site_effect(both) if both.Site.nunique()>=3 else None,
            "n_source_sampling_event_keys_present":len(counted)
        }
    # F4 has explicitly-positive recorded pollen rows, may list a moth twice
    # for two records; do NOT assume a missing F4 record represents a negative
    # pollen test. Preserve unmatched raw reference(s) as source QA.
    pollen=pollen.copy()
    pollen["ref"]=pollen.MothReference.astype(str).str.strip()
    moths["ref"]=moths.Reference.astype(str).str.strip()
    numeric=pollen.drop(columns=["MothReference","BinomialName","ref"]).apply(
        pd.to_numeric,errors="coerce")
    require(numeric.fillna(0).ge(0).all().all() and
            numeric.fillna(0).sum(axis=1).gt(0).all(),"Original pollen type counts not valid")
    pollen["n_original_recorded_grains"]=numeric.fillna(0).sum(axis=1).to_numpy()
    sums=pollen.groupby("ref",observed=True).agg(
        pollen_row_records=("ref","size"),
        grains=("n_original_recorded_grains","sum")).reset_index()
    met=sums.merge(moths[["ref","Site","Date","Method","Treatment"]],
                   how="left",on="ref",validate="one_to_one",indicator=True)
    linked=met.loc[met["_merge"].eq("both")].copy()
    require(len(linked)>0,"No independently identified moth-to-pollen source keys")
    pollen_by_method_treatment=[]
    for (method,treatment),g in linked.groupby(["Method","Treatment"],sort=True):
        pollen_by_method_treatment.append({
            "method":str(method),"treatment":str(treatment),
            "actual_matched_individual_pollen_moth_IDs":len(g),
            "source_recorded_pollen_grains":int(g.grains.sum()),
            "pollen_source_rows":int(g.pollen_row_records.sum())})
    pollen_result={
        "n_original_pollen_records":len(pollen),"n_unique_original_pollen_moth_refs":len(sums),
        "n_independent_pollen_moth_refs_linked_to_real_captured_moths":len(linked),
        "n_pollen_references_unlinked_to_captured_individuals":int(met["_merge"].ne("both").sum()),
        "total_original_recorded_grains_in_F4":int(pollen.n_original_recorded_grains.sum()),
        "method_treatment_sums_of_original_positive_pollen_records":pollen_by_method_treatment,
        "no_negative_pollen_inferred_for_F2_only_moths":True,
        "not_effective_stigma_deposition":True,
    }
    # Diagnostic: within-site dose of measured lux vs actual flight differences.
    nsite=sites.loc[sites.Site_Number.notna()].copy()
    nsite.Site_Number=nsite.Site_Number.astype(int)
    require(nsite.Site_Number.nunique()==20,"Original site metadata different")
    luxfield="Lit_site_Light_Meter_reading/lux"
    site_effect=paired.assign(logratio=np.log(paired.Lit+.5)-np.log(paired.Unlit+.5)
               ).groupby("Site").logratio.mean().rename("flight_log_ratio").reset_index()
    nsite[luxfield]=pd.to_numeric(nsite[luxfield],errors="coerce")
    dose=site_effect.merge(nsite[["Site_Number",luxfield]],
                           left_on="Site",right_on="Site_Number",validate="one_to_one")
    dose=dose.loc[dose[luxfield].notna()]
    if len(dose)>=8 and dose[luxfield].nunique()>=3:
        rho,p=spearmanr(dose[luxfield],dose.flight_log_ratio)
        lamp_dose={"n_site_pairs_with_true_metered_lux":len(dose),
                   "spearman_rho_lamp_lux_to_flight_ratio":float(rho),
                   "unadjusted_spearman_two_sided_p":float(p)}
    else:
        lamp_dose={"status":"not_enough_original_lamp_meter_values",
                   "n_site_pairs_with_true_metered_lux":len(dose)}
    result={"schema":SCHEMA,"publication_DOI":"10.1111/gcb.13371",
        "original_EIDC_dataset_DOI":"10.5285/31cc5cec-d33b-4dd6-a932-061ff947e708",
        "original_inner_CSV_SHA256":{k:v[1] for k,v in SOURCE_FILES.items()},
        "source_verified_record_counts":{k:len(v) for k,v in frames.items()},
        "overhead_flight":overhead,"original_moth_individual_capture":captures,
        "original_Pollen_grain_positive_record_integrity":pollen_result,
        "exploratory_lamp_lux_dose_for_flight":lamp_dose,
        "association_does_not_identify_plant_seed_set_or_local_moth_demography":True,
        "source_nonobserved_capture_events_not_imputed_as_zero":True,
        "independent_field_sites_not_pollen_moth_count_for_CI":True,
        "not_global_nor_island_biogeography":True}
    out=Path(out)
    require(not out.exists(),"Never overwrite a source-audited real result")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ALAN_REAL_MATCHED_LIT_DARK_MOTH_FLIGHT_AND_POLLEN",json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.out)
