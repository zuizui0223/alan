#!/usr/bin/env python3
"""Exploratory original 13-meadow-species MOIST x lighting code diagnostic.

Original Anic et al 2022 data: SHA-authentic 24-page Figshare numeric PDF,
665 source species×physical mesocosm rows, 54 truly independent 3x3x6
light×trophic plots. Author light codes L1/L2/L3 and trophic codes T1/T2/T3
are NOT mapped to unlit/amber/white/P/H/Pr until an author code key exists.
Reconstruct actual original species×cell n/means/SD through existing parser,
never transform absent source biomass rows to zero.
Independent European EIV Tichy MOIST: original source SHA verified.
Missing/unmapped species are not classified as dry or wet.

Predeclared 3 contrasts: (L2/L1,L3/L1,L3/L2), effect per species is
log(equal-trophic-cell mean biomass ratio). Primary *exploratory*
association is Spearman across source species, corrected across all
three comparisons by the max-absolute-rho trait-label permutation.
Original 54 plots are shared across species: no species-as-independent
lighting-treatment or cross-species causal inference.
"""
from __future__ import annotations
import argparse,collections,json,math,random,statistics
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
from scripts.analysis.alan_anic2022_original_pdf_biomass_54trough_parse_20261119 import (
  FILE_SHA256,FILE_ID,FIGSHARE_API,get,read_original_biomass_pdf
)
SCHEMA="alan_Anic2022_13plant_exact_EIV_MOIST_x_unmapped_light_codes_v1"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
SEED=20261120
REPS=4999
L=("L1","L2","L3")
T=("T1","T2","T3")
PAIRS=(("L1","L2"),("L1","L3"),("L2","L3"))

def require(ok,msg):
    if not ok:raise RuntimeError("ANIC13_MOISTURE_SOURCE_STOP: "+msg)

def source_eiv(path):
    from hashlib import sha256
    b=Path(path).read_bytes()
    require(sha256(b).hexdigest()==EIV_SHA,"original Europe EIV byte hash changed")
    t=pd.read_csv(path,encoding="latin1",index_col=0,low_memory=False)
    require("MOIST" in t.columns,"true EIV original MOIST column missing")
    t["MOIST"]=pd.to_numeric(t["MOIST"],errors="coerce")
    a=t.groupby(level=0).MOIST.agg(["first","nunique"])
    return a.loc[a["nunique"]<=1,"first"].dropna().to_dict()

def per_species_effects(d,eiv):
    require(d["experimental_light_trophic_full_factorial_verified"],
            "original physical 54-trough balanced treatment grid not verified")
    require(d["original_species_rows_parsed"]==665 and
            d["original_unique_physical_mesocosms_in_biomass_table"]==54,
            "original raw plant biomass table source count changed")
    groups=collections.defaultdict(dict)
    for q in d["source_original_species_by_light_trophic_cell_descriptive_means"]:
        groups[q["species"]][(q["original_light_code"],q["original_trophic_code"])]=q
    out=[]
    for name,bycell in sorted(groups.items()):
        moisture=eiv.get(name)
        m=(float(moisture) if moisture is not None and
           math.isfinite(float(moisture)) and 1<=float(moisture)<=12 else None)
        row={"species":name,"MOIST_exact":m,"wet7":m>7 if m is not None else None,
             "source_cell_counts":len(bycell),"source_original_EIV_exact_match":m is not None}
        means={}
        complete={}
        for l in L:
            c=[bycell.get((l,t)) for t in T]
            good=(all(x is not None and x["original_mesocosms_with_this_species"]>=2
                      and x["original_biomass_total_g_arithmetic_mean"]>0
                      for x in c))
            complete[l]=good
            means[l]=(statistics.fmean(x["original_biomass_total_g_arithmetic_mean"] for x in c)
                      if good else None)
        row["light_code_aggregate_means"]=means
        row["source_eligible_codes"]=[l for l in L if complete[l]]
        row["logratio_light_codes"]={
          f"{b}/{a}":math.log(means[b]/means[a]) if complete[a] and complete[b] else None
          for a,b in PAIRS
        }
        out.append(row)
    require(len(out)==13,"source 13 surviving taxa changed")
    return out

def corr(x,y):
    if len(x)<4 or len(set(x))<3 or len(set(y))<3:return None
    stat=spearmanr(x,y).statistic
    return float(stat) if math.isfinite(stat) else None

def phenotype_modulation(rows,reps=REPS):
    matched=[z for z in rows if z["MOIST_exact"] is not None]
    contrasts={}
    for a,b in PAIRS:
        col=f"{b}/{a}"
        r=[x for x in matched if x["logratio_light_codes"][col] is not None]
        wet=[x for x in r if x["wet7"]]
        dry=[x for x in r if not x["wet7"]]
        rho=corr([x["MOIST_exact"] for x in r],[x["logratio_light_codes"][col] for x in r])
        contrasts[col]={"source_species_n":len(r),
          "wet_species_n":len(wet),"drier_species_n":len(dry),
          "MOIST_min":min((x["MOIST_exact"] for x in r),default=None),
          "MOIST_max":max((x["MOIST_exact"] for x in r),default=None),
          "Spearman_MOIST_vs_log_Biomass_ratio_exploratory":rho,
          "median_logratio_wet7":statistics.median(x["logratio_light_codes"][col] for x in wet) if wet else None,
          "median_logratio_drier_le7":statistics.median(x["logratio_light_codes"][col] for x in dry) if dry else None,
          "endpoint":"total shoot+root biomass, not reproduction",
          "light_codes_unmapped_to_original_spectral_arms":True}
    # Check true data support BEFORE interpreting any nominal association.
    # Require sufficient plants across wet/dry spectrum, not just an apparent rho.
    general_eligible=(len(matched)>=10 and
         sum(x["wet7"] for x in matched)>=4 and
         sum(not x["wet7"] for x in matched)>=4 and
         max(x["MOIST_exact"] for x in matched)-min(x["MOIST_exact"] for x in matched)>=3)
    if not general_eligible:
        return {"status":"SOURCE_MOISTURE_MODERATOR_COMMON_SUPPORT_STOP",
            "matched_original_plants":len(matched),"study_original_plants":len(rows),
            "distinct_matched_wet_plant_names":sum(z["wet7"] for z in matched),
            "distinct_matched_drier_plant_names":sum(not z["wet7"] for z in matched),
            "descriptive_only":contrasts,
            "multiple_contrast_confirmatory_inference":None}
    # Complete-case across all three comparisons; same taxon set avoids
    # outcome-selective changing denominator across three null tests.
    all3=[x for x in matched if all(y is not None for y in x["logratio_light_codes"].values())]
    if len(all3)<10 or sum(z["wet7"] for z in all3)<4 or sum(not z["wet7"] for z in all3)<4:
        return {"status":"SOURCE_ALL_THREE_LIGHT_PAIR_COMMON_SUPPORT_STOP",
            "matched_original_plants":len(matched),"complete_three_pairs":len(all3),
            "descriptive_only":contrasts,"multiple_contrast_confirmatory_inference":None}
    values=[z["MOIST_exact"] for z in all3]
    effect={c:[z["logratio_light_codes"][c] for z in all3] for c in contrasts}
    observed={c:corr(values,y) for c,y in effect.items()}
    require(all(x is not None for x in observed.values()),"trait or source effects constant")
    rng=random.Random(SEED)
    stats=[]
    stats_pair={c:[] for c in effect}
    for _ in range(reps):
        shuffled=values.copy();rng.shuffle(shuffled)
        rs={c:corr(shuffled,y) for c,y in effect.items()}
        require(all(v is not None for v in rs.values()),"permutation invalid")
        stats.append(max(abs(v) for v in rs.values()))
        for c,v in rs.items():stats_pair[c].append(abs(v))
    adj={c:(1+sum(v>=abs(rho)-1e-12 for v in stats))/(reps+1)
        for c,rho in observed.items()}
    raw={c:(1+sum(v>=abs(rho)-1e-12 for v in stats_pair[c]))/(REPS+1)
        for c,rho in observed.items()}
    return {"status":"EXPLORATORY_SOURCE_TRAIT_VS_UNMAPPED_LIGHT_CODE_ASSOCIATIONS",
      "matched_original_plants":len(matched),"complete_three_pairs":len(all3),
      "wet_source_plants":sum(z["wet7"] for z in all3),
      "drier_source_plants":sum(not z["wet7"] for z in all3),
      "source_taxa_are_NOT_independent_lighting_installations":True,
      "within_same_54_trough_species_covariance_unmodeled":True,
      "partial_or_phylogenetic_moisture_effect_not_estimated":True,
      "null_species_label_exchangeability_assumption_exploratory":True,
      "permutation_replicates":reps,
      "source_species_balanced_Spearman":observed,
      "nominal_species_label_permutation_p":raw,
      "max3_absolute_rho_adjusted_permutation_p":adj,
      "descriptive_by_original_light_code":contrasts}

def execute(eiv_path):
    import json
    source=json.loads(get(FIGSHARE_API))
    meta=[x for x in source.get("files",[]) if x.get("id")==FILE_ID]
    require(len(meta)==1,"original Figshare file id missing")
    raw=get(meta[0]["download_url"])
    d=read_original_biomass_pdf(raw,strict=True)
    eiv=source_eiv(eiv_path)
    rows=per_species_effects(d,eiv)
    moderation=phenotype_modulation(rows)
    return {"schema":SCHEMA,
        "original_pdf_source_sha256":FILE_SHA256,
        "source_eiv_original_SHA256":EIV_SHA,
        "original_source_2016_biomass_665_rows":True,
        "original_physical_troughs":54,
        "original_source_studies":1,
        "original_species":len(rows),
        "original_species_exact_eiv_matches":sum(x["source_original_EIV_exact_match"] for x in rows),
        "original_species_measurements":rows,
        "effect_modifier_inference_gate":moderation,
        "wet_species_macro_effect_not_estimated":True,
        "source_light_codes_L1_L2_L3_NOT_mapped_to_dark_amber_white":True,
        "no_original_independent_true_seed_pollination_outcomes":True,
        "no_original_raw_plant_or_trough_records_committed":True}
def main():
    p=argparse.ArgumentParser();p.add_argument("--eiv",required=True);p.add_argument("--out",required=True)
    a=p.parse_args()
    r=execute(a.eiv)
    q=Path(a.out);q.parent.mkdir(exist_ok=True,parents=True)
    q.write_text(json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("ANIC2022_REAL13_MOISTURE_COVERAGE_AND_LIGHT_SOURCE_GATES "+json.dumps({
      "species":r["original_species"],"exact_eiv_matches":r["original_species_exact_eiv_matches"],
      "moderation":r["effect_modifier_inference_gate"],
      "source_species_moisture_light_code_ratios":[{
          "species":x["species"],"MOIST":x["MOIST_exact"],"codes":x["logratio_light_codes"]
        } for x in r["original_species_measurements"]]
    },sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
