#!/usr/bin/env python3
"""Source-authentic riparian ALAN aquatic-insect -> Tetragnatha -> prey audit.

Parkinson, Lawson & Tiegs (2020) PLOS ONE 15:e0240138.
Original repository LizzParkinson/ALAN_Data at commit
6985359fd29dfe321c740691bf97ea9d564db769 (source authors, not simulated).
ALAN was randomized among 20 lake mesocosms, 10 lit & 10 control.
Spider counts and web prey were repeated 3 times on the SAME mesocosms.
Pan traps classify aquatic-origin vs terrestrial-origin invertebrate counts
from ONE overnight sample. They do not prove true adult emergence, moth
pollination, local floral visits or plant fecundity.

Inference unit is physical mesocosm, not 60 visits, 1,992 individual insects,
multiple spider body-mass samples or taxonomic families. All nine observed
outcomes are from ONE experiment and must not be used as nine independent
"studies" in a meta-analysis. Label permutations (20 choose 10) are unblocked
with respect to the original fish 2x2 factor (unrecoverable fish identity);
interpret as descriptive randomization-reference sensitivities, not exact
conditional-design p-values. Randomization assumption does not identify the
role of predators as mediators; prey & predator are treatment outcomes.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

SOURCE_REPO="LizzParkinson/ALAN_Data"
SOURCE_COMMIT="6985359fd29dfe321c740691bf97ea9d564db769"
SOURCE_BLOB_SHA1={
    "Spider and Prey Insect Abundance.csv":"63965aaedae5254ea0064b9346a7fc2a8f7b610a",
    "Aquatic and Terrestrial Invertebrate Inputs.csv":"4e2b7065d67e0604ca121ea26c96b5db6d64de94",
    "Spider Body Mass.csv":"8ff2eefe6ff365e2b0f6eed038ed7b6db74247fc"
}
DOI="10.1371/journal.pone.0240138"
SCHEMA="alan_true_original_parkinson2020_riparian_predators_aquatic_insect_stage_sep_v1"
DATES=("7/8/2016","7/23/2016","8/4/2016")
PRESET_PRIMARY=("aquatic_pan","terrestrial_pan","spider_mean_3dates","prey_in_web_mean_3dates")
EXPLORATORY=("caenidae_pan","chironomidae_pan","aquatic_fraction_pan","prey_per_spider_ratio")
B=5999
SEED=20261118

def need(cond,why):
    if not bool(cond):raise ValueError("ALAN_RIPARIAN_SOURCE_STOP: "+str(why))

def git_blob_sha1(p):
    raw=Path(p).read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def source_rows(path,name):
    p=Path(path)
    need(p.is_file() and git_blob_sha1(p)==SOURCE_BLOB_SHA1[name],
         "Original author data Git blob changed: "+name)
    with p.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def grouped_source(spider_path,insects_path,mass_path):
    spiders=source_rows(spider_path,"Spider and Prey Insect Abundance.csv")
    pan=source_rows(insects_path,"Aquatic and Terrestrial Invertebrate Inputs.csv")
    mass=source_rows(mass_path,"Spider Body Mass.csv")
    need(len(spiders)==60 and len(pan)==20 and len(mass)==60,
         "20 original mesocosm x 3 visits, 20 pan traps, or 60 spider masses missing")
    by={}
    for r in spiders:
        ident=int(r["Mesocosm ID"])
        by.setdefault(ident,[]).append(r)
    need(set(by)==set(range(1,21)) and all(len(v)==3 for v in by.values()),
         "original repeated mesocosm identity damaged")
    need(len({(int(r["Mesocosm ID"]),r["Sampling Date"]) for r in spiders})==60,
         "duplicate spider sample visit/mesocosm")
    need(all(set(r["Sampling Date"] for r in v)==set(DATES) for v in by.values()),
         "incorrect three source sampling dates")
    pan_by={int(r["Mesocosm ID"]):r for r in pan}
    need(set(pan_by)==set(by),"invertebrate source mesocosm mismatch")
    labels={}
    records=[]
    for ident in sorted(by):
        pr=pan_by[ident]
        treatments={r["Treatment"] for r in by[ident]}|{pr["Treatment"]}
        need(len(treatments)==1 and treatments.issubset({"ALAN","no ALAN"}),
             "treatment switched within original mesocosm")
        lit=(next(iter(treatments))=="ALAN")
        labels[ident]=lit
        count=lambda k:float(pr[k])
        aq=count("# individual aquatic invertberates")
        land=count("# of individual terrestrial invertebrates")
        total=count("# individual invertebrates")
        caen=count("# Caenidae (A)")
        chir=count("# Chironomidae (C)")
        noct=count("# Noctuidae (m) ")
        need(all(x>=0 and float(x).is_integer() for x in (aq,land,total,caen,chir,noct)),
             "non-count data in original pan trap")
        need(aq+land==total,"aquatic/terrestrial breakdown fails original total")
        need(caen+chir<=aq,"aquatic taxa exceed aquatic total")
        need(noct<=land,"terrestrial moth count exceeds terrestrial total")
        ss=np.array([float(r["# Individual tetragnathids"]) for r in by[ident]],float)
        prey=np.array([float(r["# individual insects in webs"]) for r in by[ident]],float)
        need((ss>=0).all() and (prey>=0).all() and
             np.all(ss==np.floor(ss)) and np.all(prey==np.floor(prey)),
             "invalid repeated spider or web prey counts")
        sorted_visits=sorted(by[ident],key=lambda x:DATES.index(x["Sampling Date"]))
        records.append({
            "mesocosm":ident,"alan":lit,
            "aquatic_pan":aq,"terrestrial_pan":land,
            "caenidae_pan":caen,"chironomidae_pan":chir,"noctuidae_pan":noct,
            "total_pan":total,
            "aquatic_fraction_pan":aq/total if total>0 else float("nan"),
            "spider_mean_3dates":float(ss.mean()),
            "prey_in_web_mean_3dates":float(prey.mean()),
            "spider_3dates":[int(r["# Individual tetragnathids"]) for r in sorted_visits],
            "web_prey_3dates":[int(r["# individual insects in webs"]) for r in sorted_visits]
        })
    need(sum(labels.values())==10 and len(labels)==20,
         "original 10 ALAN/10 non-ALAN treatments not preserved")
    need(all(r["Treatment"] in ("ALAN","no ALAN") and
             r["Tetragnathid sex"] in ("F","M") for r in mass),
         "original spider sex/treatment coding changed")
    need(sum(r["total_pan"] for r in records)==1992,
         "author source total different; never silently insert paper's 1990")
    need(sum(sum(r["spider_3dates"]) for r in records)==645,
         "author original 431+214 spider observations changed")
    return records,{"spider_csv_rows":len(spiders),"pan_trap_csv_rows":len(pan),
                    "individual_body_mass_csv_rows":len(mass),
                    "independent_original_mesocosms":len(records),
                    "lit_mesocosms":10,"control_mesocosms":10,
                    "raw_pan_invertebrate_total":1992,
                    "published_article_total_stated_1990":1990,
                    "raw_sum_differs_from_published_by":2,
                    "spider_repeated_observed_count_total":645,
                    "spider_body_mass_has_no_mesocosm_ID_so_not_analysed_as_independent_unit":True}

def exact_unblocked_test(values,treatment,members=None):
    values=np.array(values,float)
    treatment=np.array(treatment,bool)
    need(len(values)==20 and treatment.sum()==10 and np.isfinite(values).all(),
         "randomization source must be 20 measured mesocosms, 10 per arm")
    idx=np.fromiter((i for pair in itertools.combinations(range(20),10) for i in pair),
                    dtype=np.int16)
    permutations=math.comb(20,10)
    ids=idx.reshape(permutations,10)
    sums=values[ids].sum(axis=1)
    obs_l=float(values[treatment].mean())
    obs_c=float(values[~treatment].mean())
    total=float(values.sum())
    obs=obs_l-obs_c
    sims=(2*sums-total)/10.0
    p=float(np.count_nonzero(np.abs(sims)>=abs(obs)-1e-12)/permutations)
    return {"observed_raw_mean_difference":obs,
            "descriptive_unblocked_fixed_10_per_arm_permutation_p":p,
            "enumerated_unblocked_label_sets":permutations,
            "not_exact_for_original_unobserved_fish_blocking":True}

def bootstrap_mean_ratio(a,b,rng,n=B):
    a,b=np.array(a,float),np.array(b,float)
    need(len(a)==10 and len(b)==10,"bootstrap must resample independent mesocosms")
    la=a[rng.integers(len(a),size=(n,len(a)))].mean(axis=1)
    lb=b[rng.integers(len(b),size=(n,len(b)))].mean(axis=1)
    valid=lb>0
    v=la[valid]/lb[valid]
    return {
       "percentile_95CI":[float(np.quantile(v,.025)),float(np.quantile(v,.975))]
           if len(v)>=.95*n else None,
       "bootstrap_n":n,
       "effective_bootstrap_comparisons":int(len(v)),
       "confidence_intervals_are_mesocosm_resampling_not_independent_study_meta":True
    }

def bh(ps):
    inds=np.argsort(ps)
    out=[float("nan")]*len(ps)
    bound=1.0
    for rank in range(len(ps)-1,-1,-1):
        j=int(inds[rank])
        bound=min(bound,float(ps[j])*len(ps)/(rank+1))
        out[j]=bound
    return out

def describe_effect(records,k,rng,primary=False):
    values=np.array([r[k] for r in records],float)
    labels=np.array([r["alan"] for r in records],bool)
    a,b=values[labels],values[~labels]
    first=float(a.mean());second=float(b.mean())
    result={"endpoint":k,"original_independent_units_per_light_arm":10,
       "lit_mean":first,"dark_mean":second,
       "lit_minus_dark":first-second,
       "lit_dark_mean_ratio":(first/second) if second>0 else None,
       "percent_change_from_control":100*(first/second-1) if second>0 else None,
       "kind":"primary_stage" if primary else "secondary_descriptive",
       "statistic_is_one_study_not_independent_meta_effect":True}
    result.update(exact_unblocked_test(values,labels))
    if second>0:
        result.update(bootstrap_mean_ratio(a,b,rng))
    return result

def within_treatment_corr(r,field):
    light=np.array([x["alan"] for x in r],bool)
    x=np.array([u[field] for u in r],float)
    y=np.array([u["spider_3dates"][-1] for u in r],float)
    def corr(a,b):
        if np.std(a)<=1e-9 or np.std(b)<=1e-9:return None
        return float(np.corrcoef(a,b)[0,1])
    raw=corr(x,y)
    x=x.copy();y=y.copy()
    for g in (False,True):
        x[light==g]-=np.mean(x[light==g])
        y[light==g]-=np.mean(y[light==g])
    return {"aquatic_field":field,
      "same_mesocosm_pan_trap_Aug4_vs_spider_count_Aug4_pearson_r":raw,
      "treatment_centered_correlational_r_not_mediation":corr(x,y),
      "independent_light_plot_count":20}

def calculate(spider,insects,mass):
    source,frame=grouped_source(spider,insects,mass)
    rng=np.random.default_rng(SEED)
    effects={}
    for k in PRESET_PRIMARY+("caenidae_pan","chironomidae_pan",
                              "aquatic_fraction_pan",):
        effects[k]=describe_effect(source,k,rng,k in PRESET_PRIMARY)
    p=[effects[k]["descriptive_unblocked_fixed_10_per_arm_permutation_p"]
       for k in PRESET_PRIMARY]
    adjusted=bh(p)
    for k,q in zip(PRESET_PRIMARY,adjusted):effects[k]["four_primary_stage_endpoints_BH_q"]=q
    # Noctuidae moth sample is only TWO individuals across all 20 plots; a
    # nocturnal plant-pollinator network or visitation response is unmeasured.
    need(sum(r["noctuidae_pan"] for r in source)==2,
         "original Noctuidae capture is sparse; source changed")
    lit=[r for r in source if r["alan"]]
    ctrl=[r for r in source if not r["alan"]]
    spider_l=sum(r["spider_mean_3dates"] for r in lit)
    spider_c=sum(r["spider_mean_3dates"] for r in ctrl)
    web_l=sum(r["prey_in_web_mean_3dates"] for r in lit)
    web_c=sum(r["prey_in_web_mean_3dates"] for r in ctrl)
    per_spider=(web_l/spider_l)/(web_c/spider_c)
    need(np.isfinite(per_spider) and per_spider>0,"prey per spider undefined")
    bs=[]
    for i in range(B):
        il=rng.integers(10,size=10)
        ic=rng.integers(10,size=10)
        sa=sum(lit[j]["spider_mean_3dates"] for j in il)
        sb=sum(ctrl[j]["spider_mean_3dates"] for j in ic)
        pa=sum(lit[j]["prey_in_web_mean_3dates"] for j in il)
        pb=sum(ctrl[j]["prey_in_web_mean_3dates"] for j in ic)
        if sa>0 and sb>0 and pb>0:
            bs.append((pa/sa)/(pb/sb))
    effects["per_capita_web_prey_not_pollination"]={
      "observed_ratio_of_pooled_mean_prey_per_spider":per_spider,
      "mesocosm_bootstrap_95CI":[float(np.quantile(bs,.025)),
                                float(np.quantile(bs,.975))],
      "not_an_attack_rate_or_independently_observed_fitness":True}
    effects["Noctuidae_terrestrial_moth_insufficient"]={
      "source_terrestrial_pan_trap_noctuid_moths_all_plots":2,
      "lit_moth_records":int(sum(r["noctuidae_pan"] for r in lit)),
      "control_moth_records":int(sum(r["noctuidae_pan"] for r in ctrl)),
      "moth_pollination_ALAN_effect":"NOT ESTIMABLE"
    }
    date_fx={}
    for ix,d in enumerate(DATES):
        v=np.array([r["spider_3dates"][ix] for r in source],float)
        w=np.array([r["web_prey_3dates"][ix] for r in source],float)
        labels=np.array([r["alan"] for r in source])
        date_fx[d]={"spider_lit_mean":float(v[labels].mean()),
                    "spider_control_mean":float(v[~labels].mean()),
                    "webprey_lit_mean":float(w[labels].mean()),
                    "webprey_control_mean":float(w[~labels].mean()),
                    "one_mesocosm_reobserved_not_new_study":True}
    out={
      "schema":SCHEMA,
      "source_original_doi":DOI,
      "source_repository":SOURCE_REPO,
      "source_fixed_commit":SOURCE_COMMIT,
      "source_blob_sha1":SOURCE_BLOB_SHA1,
      "source_units":frame,
      "predefined_four_outcomes":list(PRESET_PRIMARY),
      "sampling_dates":list(DATES),
      "source_light_design":"20 one-m2 lake mesocosms randomly assigned original lighting, repeated biweekly spiders; pan trap August 3-4 only",
      "original_fish_factor_originally_randomized_but_lost_not_a_covariate":True,
      "stage_effects":effects,
      "date_specific_descriptive":date_fx,
      "simultaneous_resource_and_predator_co_response":{
          "aquatic_pan_vs_Aug4_spider":within_treatment_corr(source,"aquatic_pan"),
          "caenidae_pan_vs_Aug4_spider":within_treatment_corr(source,"caenidae_pan")
      },
      "limitations":[
        "20 independent lighting mesocosms, NOT 60 independent survey dates",
        "One experiment: no meta-analysis of different response stages as independent studies",
        "August 4 pan-trap captures aquatic-origin adults; does NOT identify their rate of emergence vs attraction from surrounding water",
        "Spiders and pan-trap prey are outcomes jointly changed by light, not demonstrated causal mediation",
        "Body mass samples cannot be assigned to original 20 mesocosms, so no individual mass meta coefficient",
        "Noctuidae captures only 2 across all mesocosms; no flower visitation or seed/fruit output measured",
        "Original paper calls 1990 pan-trap individuals; raw deposited tables sum to 1992, audit discrepancy retained",
        "The author's text states +51% spiders, whereas raw 431/214 count ratio gives +101%; numerator/denominator definitions are not the same",
        "No observed gecko predation, no plant interactions and no mesocosm-specific plant species composition",
        "Unblocked descriptive randomization label tests ignore original abandoned fish-factor block assignment",
        "All results experimental for spider/insect outcomes ONLY, not for wetland plant fitness"
      ],
      "moth_pollen_to_plant_reproduction_identified":False,
      "independent_cross_study_meta_wetness_x_light_fitted":False
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--spider-source",required=True)
    p.add_argument("--aquatic-source",required=True)
    p.add_argument("--mass-source",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    v=calculate(a.spider_source,a.aquatic_source,a.mass_source)
    file=Path(a.out);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(v,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("PARKINSON_AUTHENTIC_RIPARIAN_ALAN_PREDATOR_RESOURCES",json.dumps(
      {"experimental_units":v["source_units"]["independent_original_mesocosms"],
       "four_outcomes":{k:{h:y for h,y in v["stage_effects"][k].items()
           if h in ("lit_mean","dark_mean","lit_dark_mean_ratio",
                    "descriptive_unblocked_fixed_10_per_arm_permutation_p",
                    "four_primary_stage_endpoints_BH_q","percentile_95CI")}
           for k in PRESET_PRIMARY},
       "per_spider":v["stage_effects"]["per_capita_web_prey_not_pollination"],
       "noctuid_moths":v["stage_effects"]["Noctuidae_terrestrial_moth_insufficient"],
       "date_effects":v["date_specific_descriptive"]},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
